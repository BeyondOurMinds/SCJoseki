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
"""

import logging
import threading
import time
from collections import OrderedDict
from contextlib import contextmanager

logger = logging.getLogger(__name__)

_lock = threading.Lock()
_trackers = {}
_results = OrderedDict()

_PHASE_KEYS = ("analysis", "processing", "visualization")


class ModuleRuntimeTracker:
    """Tracks phase and total runtimes for a single module's run callback."""

    def __init__(self, module_label):
        self.module_label = module_label
        self._phase_totals = {}
        self._callback_start = None
        self._end_to_end_start = None

    @contextmanager
    def track_callback(self):
        """Wrap an entire callback body to record total callback runtime."""
        self.begin()
        try:
            yield self
        finally:
            self.finish()

    def start_end_to_end(self):
        """Mark the moment the user-facing action began (e.g. an upstream callback)."""
        self._end_to_end_start = time.perf_counter()

    def begin(self):
        """Start timing a callback invocation (use for functions with multiple return points)."""
        self._phase_totals = {}
        self._callback_start = time.perf_counter()

    def finish(self):
        """Stop timing the current callback invocation and store the results."""
        now = time.perf_counter()
        end_to_end_start = self._end_to_end_start or self._callback_start
        self._end_to_end_start = None
        self._record(
            total_elapsed=now - self._callback_start,
            end_to_end_elapsed=now - end_to_end_start,
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
            self._phase_totals[name] = self._phase_totals.get(name, 0.0) + elapsed

    def _record(self, total_elapsed, end_to_end_elapsed):
        with _lock:
            _results[self.module_label] = {
                "analysis_runtime_s": round(self._phase_totals.get("analysis", 0.0), 4),
                "processing_runtime_s": round(self._phase_totals.get("processing", 0.0), 4),
                "visualization_runtime_s": round(self._phase_totals.get("visualization", 0.0), 4),
                "total_callback_runtime_s": round(total_elapsed, 4),
                "end_to_end_runtime_s": round(end_to_end_elapsed, 4),
            }


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
    for module_label, metrics in snapshot.items():
        lines.append(f"\n[{module_label}]")
        lines.append(f"  Analysis runtime:       {metrics['analysis_runtime_s']:.4f} s")
        lines.append(f"  Processing runtime:     {metrics['processing_runtime_s']:.4f} s")
        lines.append(f"  Visualization runtime:  {metrics['visualization_runtime_s']:.4f} s")
        lines.append(f"  Total callback runtime: {metrics['total_callback_runtime_s']:.4f} s")
        lines.append(f"  End-to-end runtime:     {metrics['end_to_end_runtime_s']:.4f} s")
    lines.append("")
    lines.append("=" * 78)
    return "\n".join(lines)


def print_runtime_report():
    """Print the runtime report for every module to the CLI (stdout + logger)."""
    report = format_runtime_report()
    print(report)
    logger.info(report)
    return report
