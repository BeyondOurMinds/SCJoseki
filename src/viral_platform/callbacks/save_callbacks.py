import logging

from dash import Input, Output, html
from dash.exceptions import PreventUpdate

from viral_platform.io.session_saves import (
    create_export_bundle,
)
from viral_platform.utils.runtime_tracker import get_tracker, print_runtime_report


logger = logging.getLogger(__name__)

_export_tracker = get_tracker("Export")


def register_save_callbacks(app):

    @app.callback(
        Output(
            "header-export-status",
            "children",
        ),
        Input(
            "header-export-btn",
            "n_clicks",
        ),
        prevent_initial_call=True,
    )
    def _export_handler(n_clicks):

        if not n_clicks:
            raise PreventUpdate

        _export_tracker.begin()

        try:
            with _export_tracker.phase("analysis"):
                result = create_export_bundle()

        except Exception as exc:
            logger.exception(
                "Full export bundle failed."
            )

            _export_tracker.finish()
            print_runtime_report()

            return html.P(
                f"Export failed: {exc}",
                style={
                    "color": "#b02a37"
                },
            )

        with _export_tracker.phase("processing"):
            status = result.get(
                "status"
            )

        if status == "cancelled":
            output = html.P(
                result.get(
                    "message",
                    "Export cancelled.",
                ),
                style={
                    "color": "#856404"
                },
            )

        elif status == "partial":
            with _export_tracker.phase("visualization"):
                details = [
                    html.P(
                        result.get(
                            "message",
                            "Export partially completed.",
                        ),
                        style={
                            "color": "#856404"
                        },
                    )
                ]

                zip_path = result.get(
                    "zip_path"
                )

                if zip_path:
                    details.append(
                        html.P(
                            f"ZIP saved to: {zip_path}",
                            style={
                                "color": "#856404"
                            },
                        )
                    )

                output = html.Div(details)

        else:
            # ---------------------------------------------------------
            # Successful export
            # ---------------------------------------------------------
            with _export_tracker.phase("visualization"):
                details = [
                    html.P(
                        "Export completed successfully.",
                        style={
                            "color": "#146c43"
                        },
                    )
                ]

                zip_path = result.get(
                    "zip_path"
                )

                h5ad_path = result.get(
                    "h5ad_path"
                )

                if zip_path:
                    details.append(
                        html.P(
                            f"ZIP: {zip_path}",
                            style={
                                "color": "#146c43"
                            },
                        )
                    )

                if h5ad_path:
                    details.append(
                        html.P(
                            f"H5AD: {h5ad_path}",
                            style={
                                "color": "#146c43"
                            },
                        )
                    )

                for note in result.get(
                    "notes",
                    [],
                ):
                    details.append(
                        html.P(
                            note
                        )
                    )

                output = html.Div(details)

        _export_tracker.finish()
        print_runtime_report()

        return output