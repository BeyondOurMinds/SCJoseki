import logging
import os
import re
import shutil
import tempfile
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path

try:
    import tkinter as tk
    from tkinter import filedialog
except ModuleNotFoundError:  # pragma: no cover - headless/container environments
    tk = None
    filedialog = None

import anndata as ad
import pandas as pd
import plotly.graph_objects as go

from viral_platform.state.dataset_store import (
    get_state_store,
    get_working_dataset,
)
from viral_platform.utils.logging_config import get_captured_logs_text


logger = logging.getLogger(__name__)


# Cached Dash component traversal

def _walk_cached_nodes(node):
    """Recursively find Plotly figures and Dash tables in cached components."""
    if node is None:
        return

    if isinstance(node, (list, tuple)):
        for item in node:
            yield from _walk_cached_nodes(item)
        return

    if hasattr(node, "to_plotly_json"):
        try:
            node = node.to_plotly_json()
        except Exception:
            return

    if isinstance(node, dict):
        # Dash component
        if (
            "type" in node
            and "props" in node
            and isinstance(node.get("props"), dict)
        ):
            props = node.get("props") or {}

            # Plotly figure
            figure_payload = props.get("figure")

            if isinstance(figure_payload, go.Figure):
                yield (
                    "figure",
                    {
                        "figure": figure_payload,
                        "component_id": props.get("id"),
                    },
                )

            elif isinstance(figure_payload, dict):
                yield (
                    "figure",
                    {
                        "figure": figure_payload,
                        "component_id": props.get("id"),
                    },
                )

            # Dash DataTable
            table_columns = props.get("columns")
            table_data = props.get("data")

            if (
                isinstance(table_columns, list)
                and isinstance(table_data, list)
            ):
                yield (
                    "table",
                    {
                        "columns": table_columns,
                        "data": table_data,
                        "component_id": props.get("id"),
                    },
                )

            # Recursively inspect other properties
            for prop_name, prop_value in props.items():
                if prop_name in {"figure", "columns", "data"}:
                    continue

                yield from _walk_cached_nodes(prop_value)

            return

        # Raw Plotly figure dictionary
        if (
            isinstance(node.get("data"), list)
            and isinstance(node.get("layout"), dict)
        ):
            yield ("figure", node)
            return

        # Actual Plotly Figure
        if isinstance(node, go.Figure):
            yield ("figure", node)
            return

        # Generic dictionary
        for value in node.values():
            yield from _walk_cached_nodes(value)


def _safe_component_name(text, fallback):
    """Convert a cache key into a safe filename component."""
    cleaned = re.sub(
        r"[^A-Za-z0-9_.-]+",
        "_",
        text or "",
    )

    cleaned = cleaned.strip("._")

    return cleaned or fallback


def _cache_key_export_prefix(cache_key, item_type):
    """Map cache keys to readable export filename prefixes."""
    explicit = {
        "de-table-container": "de_table",
        "volcano-plot-container": "volcano_plot",
        "de-heatmap-container": "de_heatmap",
    }
    prefix = explicit.get(cache_key)
    if prefix:
        return prefix
    default_prefix = _safe_component_name(cache_key, item_type)
    return default_prefix.replace("-", "_")


def _extract_celltype_label(component_id):
    """Extract a usable cell-type label from a Dash component id if available."""
    if isinstance(component_id, dict):
        celltype = component_id.get("celltype")
        if celltype is not None and str(celltype).strip():
            return _safe_component_name(str(celltype), "celltype")
    return None


def _extract_tables_and_figures(results_cache):
    """Extract tables and Plotly figures from cached Dash components."""
    tables = []
    figures = []

    for cache_key, component in (results_cache or {}).items():
        table_index = 1
        figure_index = 1
        used_table_names = set()
        used_figure_names = set()

        for item_type, payload in _walk_cached_nodes(component):

            if item_type == "table":
                columns = payload.get("columns") or []
                data = payload.get("data") or []

                if not isinstance(data, list):
                    continue

                base_name = _cache_key_export_prefix(cache_key, "table")
                celltype_label = _extract_celltype_label(
                    payload.get("component_id")
                )
                if celltype_label:
                    table_id = f"{base_name}_{celltype_label}"
                else:
                    table_id = f"{base_name}_table_{table_index}"

                if table_id in used_table_names:
                    table_id = f"{table_id}_{table_index}"
                used_table_names.add(table_id)

                table_index += 1

                df = pd.DataFrame(data)

                column_ids = [
                    column.get("id")
                    for column in columns
                    if (
                        isinstance(column, dict)
                        and column.get("id")
                    )
                ]

                if (
                    column_ids
                    and all(
                        column in df.columns
                        for column in column_ids
                    )
                ):
                    df = df[column_ids]

                tables.append(
                    (
                        table_id,
                        df,
                    )
                )

            elif item_type == "figure":
                base_name = _cache_key_export_prefix(cache_key, "figure")
                component_id = payload.get("component_id") if isinstance(payload, dict) else None
                figure_payload = payload.get("figure") if isinstance(payload, dict) else payload
                celltype_label = _extract_celltype_label(component_id)
                if celltype_label:
                    figure_id = f"{base_name}_{celltype_label}"
                else:
                    figure_id = f"{base_name}_figure_{figure_index}"

                if figure_id in used_figure_names:
                    figure_id = f"{figure_id}_{figure_index}"
                used_figure_names.add(figure_id)

                figure_index += 1

                figures.append(
                    (
                        figure_id,
                        figure_payload,
                    )
                )

    return tables, figures


# Plotly export
def _prepare_figure_for_export(figure_like):
    """Prepare a Plotly figure for readable SVG export."""
    if isinstance(figure_like, go.Figure):
        fig = go.Figure(figure_like)
    else:
        fig = go.Figure(figure_like)

    legend = fig.layout.legend

    if (
        legend
        and getattr(legend, "orientation", None) == "h"
    ):
        current_margin = fig.layout.margin or {}

        bottom_margin = max(
            int(
                getattr(
                    current_margin,
                    "b",
                    0,
                )
                or 0
            ),
            220,
        )

        fig.update_layout(
            margin=dict(
                l=int(
                    getattr(
                        current_margin,
                        "l",
                        0,
                    )
                    or 0
                ),
                r=int(
                    getattr(
                        current_margin,
                        "r",
                        0,
                    )
                    or 0
                ),
                t=int(
                    getattr(
                        current_margin,
                        "t",
                        60,
                    )
                    or 60
                ),
                b=bottom_margin,
            ),
            width=max(
                int(
                    getattr(
                        fig.layout,
                        "width",
                        0,
                    )
                    or 0
                ),
                1600,
            ),
            height=max(
                int(
                    getattr(
                        fig.layout,
                        "height",
                        0,
                    )
                    or 0
                ),
                900,
            ),
        )

    return fig


# H5AD
def _write_h5ad_file(adata, path):
    """Write AnnData directly to the requested H5AD path."""
    setting_name = "allow_write_nullable_strings"

    if hasattr(ad.settings, setting_name):
        previous_setting = getattr(
            ad.settings,
            setting_name,
        )

        try:
            setattr(
                ad.settings,
                setting_name,
                True,
            )

            adata.write_h5ad(str(path))

        finally:
            setattr(
                ad.settings,
                setting_name,
                previous_setting,
            )

    else:
        adata.write_h5ad(str(path))


# Log export
def _build_log_text():
    """Build the log file included in the ZIP export."""
    adata = get_working_dataset()

    captured_logs = get_captured_logs_text().strip()

    lines = [
        "SCJoseki export",
        f"created_utc: {datetime.now(timezone.utc).isoformat()}",
        f"cells: {int(adata.n_obs) if adata is not None else 0}",
        f"genes: {int(adata.n_vars) if adata is not None else 0}",
        "",
        "captured_logger_output:",
        (
            captured_logs
            if captured_logs
            else "<no logger messages captured>"
        ),
    ]

    return "\n".join(lines)


def _resolve_export_directory():
    """Return a writable directory for exported files.

    In Docker/headless Linux environments, callers can mount a host
    volume at /exports or set SCJOSEKI_EXPORT_DIR to any writable path.
    On desktop platforms, return None so that a native save dialog
    can be displayed.
    """
    override_dir = os.getenv("SCJOSEKI_EXPORT_DIR")

    if override_dir:
        export_dir = Path(override_dir).expanduser()
        export_dir.mkdir(parents=True, exist_ok=True)
        return export_dir

    # Linux/headless environments without a display
    if (
        os.name != "nt"
        and not os.environ.get("DISPLAY")
    ):
        export_dir = Path("/exports")
        export_dir.mkdir(parents=True, exist_ok=True)
        return export_dir

    # Desktop environments: use native save dialog
    return None


# Native save dialog when a GUI is available; otherwise fall back to a filesystem path.
def _choose_save_path(
    title,
    initial_filename,
    filetypes,
):
    """
    Open a native save dialog when Tk is available.

    In headless/container runs, fall back to a writable export directory so the
    export flow remains usable without a local GUI display.
    """
    export_dir = _resolve_export_directory()
    if export_dir is not None:
        return export_dir / initial_filename

    try:
        root = tk.Tk()
    except Exception:
        fallback_dir = Path("/exports")
        fallback_dir.mkdir(parents=True, exist_ok=True)
        return fallback_dir / initial_filename

    root.withdraw()
    root.attributes("-topmost", True)

    try:
        selected_path = filedialog.asksaveasfilename(
            title=title,
            initialfile=initial_filename,
            filetypes=filetypes,
            defaultextension=filetypes[0][1],
        )
    except Exception:
        selected_path = ""
    finally:
        try:
            root.destroy()
        except Exception:
            pass

    if not selected_path:
        return None

    return Path(selected_path)


# ZIP creation
def _create_zip_without_h5ad(adata):
    """
    Create the results ZIP without the H5AD.

    The returned ZIP is stored in a temporary file rather than
    returned as a giant in-memory Dash payload.
    """
    start_time = time.perf_counter()

    temp_dir = tempfile.mkdtemp(
        prefix="scjoseki_export_"
    )

    zip_path = Path(temp_dir) / (
        "scjoseki_export_"
        f"{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        ".zip"
    )

    file_payloads = {}
    notes = []

    # Metadata
    metadata_start = time.perf_counter()

    metadata_csv = (
        adata.obs
        .copy()
        .to_csv(index=True)
    )

    file_payloads[
        "metadata/metadata.csv"
    ] = metadata_csv.encode("utf-8")

    logger.info(
        "Metadata export completed in %.2f seconds (%d bytes).",
        time.perf_counter() - metadata_start,
        len(file_payloads["metadata/metadata.csv"]),
    )

    # Tables and figures
    results_cache = get_state_store().get(
        "results_cache",
        {},
    )

    tables, figures = _extract_tables_and_figures(
        results_cache
    )

    logger.info(
        "Export contains %d tables and %d figures.",
        len(tables),
        len(figures),
    )

    # Tables
    table_start = time.perf_counter()

    for table_name, table_df in tables:
        file_payloads[
            f"tables/{table_name}.csv"
        ] = table_df.to_csv(
            index=False
        ).encode("utf-8")

    logger.info(
        "Table export completed in %.2f seconds.",
        time.perf_counter() - table_start,
    )

    # Figures
    exported_figure_count = 0

    for figure_name, figure_dict in figures:
        try:
            logger.info(
                "Starting SVG export: %s",
                figure_name,
            )

            start_time = time.perf_counter()

            export_fig = _prepare_figure_for_export(
                figure_dict
            )

            logger.info(
                "Prepared figure '%s' in %.2f seconds.",
                figure_name,
                time.perf_counter() - start_time,
            )

            if getattr(
                export_fig,
                "data",
                None,
            ) is None:
                raise ValueError(
                    "Figure has no data."
                )

            svg_start = time.perf_counter()

            svg_bytes = export_fig.to_image(
                format="svg"
            )

            logger.info(
                "Rendered '%s' to SVG in %.2f seconds (%d bytes).",
                figure_name,
                time.perf_counter() - svg_start,
                len(svg_bytes),
            )

            file_payloads[
                f"figures/{figure_name}.svg"
            ] = svg_bytes

            exported_figure_count += 1

        except Exception as exc:
            logger.exception(
                "Could not export figure '%s' as SVG.",
                figure_name,
            )

            notes.append(
                f"Could not export figure "
                f"'{figure_name}' as SVG: {exc}"
            )

    notes.append(
        f"Exported {len(tables)} table(s)."
    )

    notes.append(
        f"Exported {exported_figure_count} figure(s) as SVG."
    )

    # Log
    log_text = _build_log_text()

    file_payloads[
        "logs/session_export.log"
    ] = log_text.encode("utf-8")

    # Write ZIP to disk
    zip_start = time.perf_counter()

    with zipfile.ZipFile(
        zip_path,
        mode="w",
        compression=zipfile.ZIP_DEFLATED,
    ) as archive:

        for relative_path, payload in file_payloads.items():
            archive.writestr(
                relative_path,
                payload,
                compress_type=zipfile.ZIP_DEFLATED,
            )

    zip_size = zip_path.stat().st_size

    logger.info(
        "ZIP creation completed in %.2f seconds.",
        time.perf_counter() - zip_start,
    )

    logger.info(
        "ZIP size: %.2f MB.",
        zip_size / (1024 * 1024),
    )

    logger.info(
        "ZIP preparation completed in %.2f seconds.",
        time.perf_counter() - start_time,
    )

    return zip_path, notes, temp_dir


# Complete export workflow
def create_export_bundle():
    """
    Create and save the complete SCJoseki export.

    Workflow:

    1. Create ZIP containing:
       - metadata CSV
       - tables
       - SVG figures
       - log

    2. Prompt the user for a ZIP save location.

    3. Copy the ZIP to that location.

    4. Prompt the user for an H5AD save location.

    5. Write the current AnnData directly to that location.

    Returns
    -------
    dict
        Export result information.
    """
    total_start = time.perf_counter()

    adata = get_working_dataset()

    if adata is None:
        raise ValueError(
            "No working dataset is loaded. "
            "Upload data before exporting."
        )

    zip_path = None
    temp_dir = None

    try:
        # Step 1: Create ZIP without H5AD
        logger.info(
            "Starting SCJoseki export."
        )

        zip_path, notes, temp_dir = (
            _create_zip_without_h5ad(adata)
        )

        # Step 2: Ask user where to save ZIP
        zip_filename = zip_path.name

        selected_zip_path = _choose_save_path(
            title="Save SCJoseki results ZIP",
            initial_filename=zip_filename,
            filetypes=[
                (
                    "SCJoseki ZIP archive",
                    "*.zip",
                ),
                (
                    "All files",
                    "*.*",
                ),
            ],
        )

        if selected_zip_path is None:
            return {
                "status": "cancelled",
                "message": (
                    "Export cancelled before the ZIP "
                    "was saved."
                ),
                "zip_path": None,
                "h5ad_path": None,
                "notes": notes,
            }

        # Step 3: Copy ZIP to selected location
        logger.info(
            "Saving ZIP to: %s",
            selected_zip_path,
        )

        shutil.copyfile(
            zip_path,
            selected_zip_path,
        )

        logger.info(
            "ZIP saved successfully."
        )

        # Step 4: Ask user where to save H5AD
        selected_h5ad_path = _choose_save_path(
            title="Save SCJoseki H5AD dataset",
            initial_filename="session.h5ad",
            filetypes=[
                (
                    "AnnData H5AD file",
                    "*.h5ad",
                ),
                (
                    "All files",
                    "*.*",
                ),
            ],
        )

        if selected_h5ad_path is None:
            return {
                "status": "partial",
                "message": (
                    "The results ZIP was saved, "
                    "but the H5AD save was cancelled."
                ),
                "zip_path": selected_zip_path,
                "h5ad_path": None,
                "notes": notes,
            }

        # Step 5: Write H5AD directly to selected location
        logger.info(
            "Saving H5AD to: %s",
            selected_h5ad_path,
        )

        h5_start = time.perf_counter()

        _write_h5ad_file(
            adata,
            selected_h5ad_path,
        )

        h5_size = selected_h5ad_path.stat().st_size

        logger.info(
            "H5AD saved successfully in %.2f seconds (%.2f MB).",
            time.perf_counter() - h5_start,
            h5_size / (1024 * 1024),
        )

        logger.info(
            "Total export workflow completed in %.2f seconds.",
            time.perf_counter() - total_start,
        )

        return {
            "status": "success",
            "message": (
                "Export completed successfully."
            ),
            "zip_path": selected_zip_path,
            "h5ad_path": selected_h5ad_path,
            "notes": notes,
        }

    finally:
        # The temporary ZIP is only needed until the user chooses its final destination.
        if temp_dir is not None:
            try:
                shutil.rmtree(
                    temp_dir,
                    ignore_errors=True,
                )
            except Exception:
                logger.exception(
                    "Could not clean up temporary export directory."
                )