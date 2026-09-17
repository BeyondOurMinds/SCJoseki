import os
import uuid
import logging
from pathlib import Path

from dash import Input, Output, State, html, no_update

from viral_platform.app import UPLOAD_FOLDER
from viral_platform.io.loaders import load_file_from_path
from viral_platform.utils.runtime_tracker import get_tracker

logger = logging.getLogger(__name__)

_upload_tracker = get_tracker("File Upload")


def _cleanup_uploaded_cache_file(file_path):
    """Remove a processed upload file and empty cache folders left behind."""
    path = Path(file_path)
    root = Path(UPLOAD_FOLDER).resolve()

    try:
        if path.exists() and path.is_file():
            path.unlink()
            logger.info("Removed processed upload cache file: %s", path)
    except OSError as exc:
        logger.warning("Could not remove upload cache file %s: %s", path, exc)
        return

    # Remove now-empty parent directories until we reach the upload cache root.
    current = path.parent
    while current != root:
        try:
            current.rmdir()
            logger.info("Removed empty upload cache directory: %s", current)
        except OSError:
            break
        current = current.parent


def register_upload_callbacks(app):
    @app.callback(
        Output("output-data-upload", "children"),
        Output("active-dataset-version", "data"),
        Input("file-uploader", "isCompleted"),
        State("file-uploader", "fileNames"),
        State("file-uploader", "upload_id"),
        prevent_initial_call=False,
    )
    def process_uploaded_file(is_completed, file_names, upload_id):
        if not is_completed:
            logger.info("Upload callback received non-complete state.")
            return no_update, no_update

        _upload_tracker.begin()

        root_folder = Path(UPLOAD_FOLDER)
        if upload_id:
            root_folder = root_folder / str(upload_id)

        file_paths = []
        if file_names:
            file_paths = [str(root_folder / filename) for filename in file_names]

        if not file_paths:
            logger.info("Upload callback triggered with no files.")
            _upload_tracker.finish()
            return html.P("No file uploaded yet."), no_update

        file_path = file_paths[0]
        filename = os.path.basename(file_path)
        logger.info("Chunked upload completed: %s", file_path)

        try:
            with _upload_tracker.phase("analysis"):
                adata = load_file_from_path(file_path)
        except Exception as exc:
            logger.exception("Failed to process uploaded file: %s", filename)
            _upload_tracker.finish()
            return (
                html.Div([
                    html.H5(f"Uploaded file: {filename}"),
                    html.P(f"Failed to read file: {exc}"),
                ]),
                no_update,
            )

        _cleanup_uploaded_cache_file(file_path)

        logger.info(
            "Successfully processed uploaded file %s (%s cells, %s genes).",
            filename,
            adata.n_obs,
            adata.n_vars,
        )

        with _upload_tracker.phase("processing"):
            sample_count = adata.uns.get("sample_count")
            dataset_message = f"Loaded dataset with {adata.n_obs} cells and {adata.n_vars} genes."
            if sample_count:
                dataset_message = (
                    f"Loaded dataset with {adata.n_obs} cells and {adata.n_vars} genes "
                    f"from {sample_count} sample(s)."
                )

        with _upload_tracker.phase("visualization"):
            output = html.Div([
                html.H5(f"Uploaded file: {filename}"),
                html.P(dataset_message),
            ])

        _upload_tracker.finish()

        return (
            output,
            str(uuid.uuid4()),
        )

    @app.callback(
        Output("file-input", "children"),
        Input("output-data-upload", "children"),
        prevent_initial_call=False,
    )
    def update_file_status(_upload_output):
        if not _upload_output:
            return "No file uploaded"

        logger.info("Upload status updated in UI.")
        return "File processed. See upload result below."