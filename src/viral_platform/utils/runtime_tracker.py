"""Runtime instrumentation shared by every module's "run" callback.

Each module wraps its run-button callback body in ``tracker.track_callback()``
and marks sub-regions with ``tracker.phase("analysis"|"processing"|"visualization")``.
Results are stored in a process-wide registry and printed to the CLI once the
export button's callback finishes (see ``print_runtime_report``).

Phase definitions:
    analysis       -- run button pressed -> raw analysis/computation finished
                       (does not include figure/table building).
    processing     -- converting raw analysis results into figures/tables.
    visualization  -- wrapping those figures/tables into the Dash components
                       that get rendered in the app.
    total_callback -- wall-clock time of the entire callback function body.
    end_to_end     -- user clicks run -> results are returned for rendering.
                       For single-callback modules this equals total_callback;
                       for multi-callback flows (e.g. upload) it can span
                       multiple chained callbacks.

Each callback invocation is also sampled in the background for process CPU%%,
process RAM, and (if an NVIDIA GPU + pynvml are available) GPU utilization/
memory, reporting the peak and average observed while the callback ran.
"""

import logging
import os
import threading
import time
from collections import OrderedDict
from contextlib import contextmanager

logger = logging.getLogger(__name__)

try:
    import psutil
except ImportError:  # pragma: no cover - psutil is a declared dependency
    psutil = None

_GPU_HANDLES = None
# psutil.cpu_percent() is normalized to a single core (100% == 1 core fully busy),
# so a process using N cores can report up to N*100%. Divide by the logical core
# count to report the more intuitive "% of total machine capacity" (0-100 scale).
_LOGICAL_CPU_COUNT = (psutil.cpu_count(logical=True) or 1) if psutil else 1


def _get_gpu_handles():
    """Lazily initialize NVML and return device handles, or None if unavailable."""
    global _GPU_HANDLES
    if _GPU_HANDLES is not None:
        return _GPU_HANDLES or None
    try:
        import pynvml

        pynvml.nvmlInit()
        device_count = pynvml.nvmlDeviceGetCount()
        _GPU_HANDLES = [pynvml.nvmlDeviceGetHandleByIndex(i) for i in range(device_count)]
    except Exception:
        _GPU_HANDLES = []
    return _GPU_HANDLES or None


class _ResourceSampler:
    """Background sampler recording peak/average CPU, RAM, and GPU usage."""

    def __init__(self, interval=0.1):
        self.interval = interval
        self._process = psutil.Process(os.getpid()) if psutil else None
        self._thread = None
        self._stop_event = threading.Event()
        self.cpu_samples = []
        self.ram_samples_mb = []
        self.gpu_util_samples = []
        self.gpu_mem_samples_mb = []

    def start(self):
        if self._process is None:
            return
        self._process.cpu_percent(None)  # prime the internal counter, first read is meaningless
        self._stop_event.clear()
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def _run(self):
        gpu_handles = _get_gpu_handles()
        while not self._stop_event.wait(self.interval):
            # Waiting before sampling (rather than after) keeps each cpu_percent()
            # window close to `self.interval`; sub-tick windows produce wildly
            # inflated readings on Windows due to limited clock-tick resolution.
            try:
                raw_cpu_percent = self._process.cpu_percent(None)
                self.cpu_samples.append(raw_cpu_percent / _LOGICAL_CPU_COUNT)
                self.ram_samples_mb.append(self._process.memory_info().rss / (1024 * 1024))
                if gpu_handles:
                    import pynvml

                    util_total = 0
                    mem_total = 0.0
                    for handle in gpu_handles:
                        util_total += pynvml.nvmlDeviceGetUtilizationRates(handle).gpu
                        mem_total += pynvml.nvmlDeviceGetMemoryInfo(handle).used / (1024 * 1024)
                    self.gpu_util_samples.append(util_total / len(gpu_handles))
                    self.gpu_mem_samples_mb.append(mem_total)
            except Exception:
                pass

    def stop(self):
        """Stop sampling and return a dict of peak/average metrics."""
        if self._process is None:
            return _EMPTY_RESOURCE_METRICS.copy()
        self._stop_event.set()
        if self._thread is not None:
            self._thread.join(timeout=1)
        return {
            "cpu_peak_percent": max(self.cpu_samples, default=0.0),
            "cpu_avg_percent": _mean(self.cpu_samples),
            "ram_peak_mb": max(self.ram_samples_mb, default=0.0),
            "ram_avg_mb": _mean(self.ram_samples_mb),
            "gpu_peak_percent": max(self.gpu_util_samples, default=None) if self.gpu_util_samples else None,
            "gpu_avg_percent": _mean(self.gpu_util_samples) if self.gpu_util_samples else None,
            "gpu_peak_mem_mb": max(self.gpu_mem_samples_mb, default=None) if self.gpu_mem_samples_mb else None,
            "gpu_avg_mem_mb": _mean(self.gpu_mem_samples_mb) if self.gpu_mem_samples_mb else None,
        }


_EMPTY_RESOURCE_METRICS = {
    "cpu_peak_percent": None,
    "cpu_avg_percent": None,
    "ram_peak_mb": None,
    "ram_avg_mb": None,
    "gpu_peak_percent": None,
    "gpu_avg_percent": None,
    "gpu_peak_mem_mb": None,
    "gpu_avg_mem_mb": None,
}


def _mean(values):
    return (sum(values) / len(values)) if values else 0.0


_lock = threading.Lock()
_trackers = {}
_results = OrderedDict()

_PHASE_KEYS = ("analysis", "processing", "visualization")


class _Invocation:
    """State for one callback invocation; kept per-thread so concurrent runs don't mix."""

    def __init__(self):
        self.phase_totals = {}
        self.start = time.perf_counter()
        self.sampler = _ResourceSampler()
        self.sampler.start()


class ModuleRuntimeTracker:
    """Tracks phase and total runtimes for a single module's run callback."""

    def __init__(self, module_label):
        self.module_label = module_label
        self._local = threading.local()

    @contextmanager
    def track_callback(self):
        """Wrap an entire callback body to record total callback runtime."""
        self.begin()
        try:
            yield self
        finally:
            self.finish()

    def begin(self):
        """Start timing a callback invocation (use for functions with multiple return points)."""
        self._local.invocation = _Invocation()

    def finish(self):
        """Stop timing the current callback invocation and store the results."""
        invocation = getattr(self._local, "invocation", None)
        self._local.invocation = None
        if invocation is None:
            return
        now = time.perf_counter()
        resource_metrics = invocation.sampler.stop()
        # Runs that never entered a phase are early returns (e.g. initial page load); don't overwrite real results.
        if not invocation.phase_totals:
            return
        elapsed = now - invocation.start
        self._record(
            phase_totals=invocation.phase_totals,
            total_elapsed=elapsed,
            end_to_end_elapsed=elapsed,
            resource_metrics=resource_metrics,
        )

    def track(self, func):
        """Decorator that times a callback regardless of which return path it takes."""
        import functools

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            self.begin()
            try:
                return func(*args, **kwargs)
            finally:
                self.finish()

        return wrapper

    @contextmanager
    def phase(self, name):
        if name not in _PHASE_KEYS:
            raise ValueError(f"Unknown runtime phase '{name}', expected one of {_PHASE_KEYS}")
        start = time.perf_counter()
        try:
            yield
        finally:
            elapsed = time.perf_counter() - start
            invocation = getattr(self._local, "invocation", None)
            if invocation is not None:
                invocation.phase_totals[name] = invocation.phase_totals.get(name, 0.0) + elapsed

    def _record(self, phase_totals, total_elapsed, end_to_end_elapsed, resource_metrics):
        with _lock:
            _results[self.module_label] = {
                "analysis_runtime_s": round(phase_totals.get("analysis", 0.0), 4),
                "processing_runtime_s": round(phase_totals.get("processing", 0.0), 4),
                "visualization_runtime_s": round(phase_totals.get("visualization", 0.0), 4),
                "total_callback_runtime_s": round(total_elapsed, 4),
                "end_to_end_runtime_s": round(end_to_end_elapsed, 4),
                **_round_resource_metrics(resource_metrics),
            }


def _round_resource_metrics(metrics):
    return {key: (round(value, 2) if value is not None else None) for key, value in metrics.items()}


def get_tracker(module_label):
    """Return the process-wide tracker instance for a given module label."""
    with _lock:
        tracker = _trackers.get(module_label)
        if tracker is None:
            tracker = ModuleRuntimeTracker(module_label)
            _trackers[module_label] = tracker
        return tracker


def get_all_results():
    """Return a snapshot of every module's recorded runtimes."""
    with _lock:
        return OrderedDict(_results)


def format_runtime_report():
    """Build a human-readable CLI report of every recorded module's runtimes."""
    with _lock:
        snapshot = OrderedDict(_results)

    if not snapshot:
        return "No module runtimes have been recorded yet."

    lines = ["", "=" * 78, "MODULE RUNTIME REPORT", "=" * 78]
    if psutil is None:
        lines.append("(psutil not installed - CPU/RAM tracking unavailable)")
    if not _get_gpu_handles():
        lines.append("(no NVIDIA GPU / pynvml detected - GPU tracking unavailable)")

    for module_label, metrics in snapshot.items():
        lines.append(f"\n[{module_label}]")
        lines.append(f"  Analysis runtime:       {metrics['analysis_runtime_s']:.4f} s")
        lines.append(f"  Processing runtime:     {metrics['processing_runtime_s']:.4f} s")
        lines.append(f"  Visualization runtime:  {metrics['visualization_runtime_s']:.4f} s")
        lines.append(f"  Total callback runtime: {metrics['total_callback_runtime_s']:.4f} s")
        lines.append(f"  End-to-end runtime:     {metrics['end_to_end_runtime_s']:.4f} s")
        lines.append(f"  CPU usage (peak/avg):   {_fmt(metrics['cpu_peak_percent'])} / {_fmt(metrics['cpu_avg_percent'])} % of total machine capacity")
        lines.append(f"  RAM usage (peak/avg):   {_fmt(metrics['ram_peak_mb'])} / {_fmt(metrics['ram_avg_mb'])} MB")
        lines.append(f"  GPU usage (peak/avg):   {_fmt(metrics['gpu_peak_percent'])} / {_fmt(metrics['gpu_avg_percent'])} %")
        lines.append(f"  GPU memory (peak/avg):  {_fmt(metrics['gpu_peak_mem_mb'])} / {_fmt(metrics['gpu_avg_mem_mb'])} MB")
    lines.append("")
    lines.append("=" * 78)
    return "\n".join(lines)


def _fmt(value):
    return f"{value:.2f}" if value is not None else "N/A"



def print_runtime_report():
    """Print the runtime report for every module to the CLI (stdout + logger)."""
    report = format_runtime_report()
    print(report)
    logger.info(report)
    return report
