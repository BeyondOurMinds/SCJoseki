import re

from dash import Input, Output, State, ctx, dcc, html, no_update
import dash_bootstrap_components as dbc

from viral_platform.analysis.viral_gene_detection import (
	find_custom_viral_genes,
	find_viral_genes,
	list_viral_gene_sets,
	normalize_gene_name,
)
from viral_platform.state.dataset_store import cache_results, get_dataset, get_state_store, get_working_dataset, update_state_store
from viral_platform.utils.runtime_tracker import get_tracker

_viral_detection_tracker = get_tracker("Viral Gene Detection")


def help_icon():
	"""Create the shared help icon used in viral detection UI rows.

	Use:
	- Standardizes helper icon style in result/action panels.

	Interacts with:
	- create_viral_gene_detection_results.

	Inputs:
	- None.

	Outputs:
	- dash.html.Span component.
	"""
	return html.Span(
		"?",
		style={
			"display": "inline-flex",
			"alignItems": "center",
			"justifyContent": "center",
			"width": "18px",
			"height": "18px",
			"borderRadius": "50%",
			"border": "1.5px solid #6c757d",
			"fontSize": "11px",
			"color": "#6c757d",
			"cursor": "pointer",
			"marginLeft": "6px",
			"verticalAlign": "middle",
		},
	)


def virus_card(icon_char, name, badge_text, detail):
	"""Render one per-virus summary card.

	Use:
	- Displays per-virus match counts and context text.

	Interacts with:
	- create_viral_gene_detection_results.

	Inputs:
	- icon_char (str), name (str), badge_text (str), detail (Dash child).

	Outputs:
	- dash.html.Div card component.
	"""
	return html.Div(
		[
			html.Div(
				[
					html.Span(icon_char, style={"fontSize": "28px", "marginRight": "12px"}),
					html.Div(
						[
							html.Div(
								[
									html.A(name, href="#", style={"fontWeight": "600", "color": "#0d6efd", "textDecoration": "none"}),
									dbc.Badge(badge_text, color="secondary", style={"marginLeft": "8px", "fontSize": "11px"}),
								],
								style={"display": "flex", "alignItems": "center"},
							),
							html.P(detail, style={"margin": "2px 0 0", "fontSize": "13px", "color": "#6c757d"}),
						],
					),
				],
				style={"display": "flex", "alignItems": "center"},
			),
		],
		style={
			"padding": "12px 16px",
			"border": "1px solid #dee2e6",
			"borderRadius": "6px",
			"backgroundColor": "#fff",
			"height": "100%",
		},
	)


def _split_input_genes(raw_text):
	"""Split add-gene input text into cleaned gene tokens.

	Use:
	- Handles comma/newline user input for append operations.

	Interacts with:
	- curate_viral_gene_list.

	Inputs:
	- raw_text (str|None).

	Outputs:
	- list[str]: trimmed gene tokens.
	"""
	if not raw_text:
		return []
	return [token.strip() for token in re.split(r"[,\n\r]+", raw_text) if token.strip()]


def _build_viral_gene_entries(genes, matched_features_by_gene):
	"""Build dropdown/list entries for detected viral genes.

	Use:
	- Formats editable labels and includes matched-as only for alias-style matches.

	Interacts with:
	- create_viral_gene_detection_results, curate_viral_gene_list.

	Inputs:
	- genes (iterable[str])
	- matched_features_by_gene (dict[str, iterable[str]])

	Outputs:
	- list[dict]: dropdown entries with label/value.
	"""
	entries = []
	for gene in sorted(set(genes)):
		matched = sorted(matched_features_by_gene.get(gene, []))
		alias_matches = [
			feature
			for feature in matched
			if normalize_gene_name(feature).upper() != gene.upper()
		]
		label = f"{gene} (matched as {', '.join(alias_matches)})" if alias_matches else gene
		entries.append({"label": label, "value": gene})
	return entries


def _build_gene_options(genes):
	"""Build standardized dropdown options for normalized gene names.

	Use:
	- Keeps host-virus interaction dropdown options in a consistent label/value format.

	Interacts with:
	- run_viral_gene_detection, curate_viral_gene_list.

	Inputs:
	- genes (iterable[str]).

	Outputs:
	- list[dict]: options as {"label": gene, "value": gene}.
	"""
	return [{"label": gene, "value": gene} for gene in sorted(set(genes))]


def _format_detected_grouped_list(set_to_genes, empty_message):
	"""Format grouped detected-gene text as '(set): gene1, gene2' chunks.

	Use:
	- Produces display-only grouped list text for detection result accordions.

	Interacts with:
	- create_viral_gene_detection_results, curate_viral_gene_list.

	Inputs:
	- set_to_genes (dict[str, iterable[str]]), empty_message (str).

	Outputs:
	- str: grouped list text or empty_message.
	"""
	parts = []
	for set_name in sorted((set_to_genes or {}).keys()):
		genes = sorted(set(set_to_genes.get(set_name, [])))
		if genes:
			parts.append(f"({set_name}): {', '.join(genes)}")
	return " ".join(parts) if parts else empty_message


def _normalize_automatic_detection_payload(detected, selected_virus):
	"""Return a per-virus mapping regardless of finder return shape."""
	if not isinstance(detected, dict):
		return {}

	# Single-virus payload: {"features": [...], "genes": [...], ...}
	if "features" in detected and "genes" in detected:
		virus_key = selected_virus or "Selected Virus"
		return {virus_key: detected}

	# All-virus payload: {"Virus A": {...}, "Virus B": {...}}
	return detected


def create_viral_gene_detection_results(pass_fail, color, gene_count_per_virus, detected_features, gene_entries, unique_count, not_found=None, detected_gene_sets=None):
	"""Render full viral detection results UI with curation controls.

	Use:
	- Builds summary cards, detected-list editor, append controls, and alert output.

	Interacts with:
	- virus_card, _build_viral_gene_entries, curate_viral_gene_list callback.

	Inputs:
	- pass_fail/color: status display values.
	- gene_count_per_virus/detected_features: aggregate metrics.
	- gene_entries: dropdown/list entries.
	- unique_count (int), not_found (list[str]|None).

	Outputs:
	- dash.html.Div containing interactive viral detection results.
	"""
	all_genes_count = unique_count if unique_count else 0
	grouped_list_text = _format_detected_grouped_list(
		detected_gene_sets,
		"No viral genes currently selected.",
	)
	result = html.Div(
		[
			dbc.Row(
				[
					dbc.Col(
						html.Div(
							[
								html.Span(
									f"{pass_fail}",
									style={
										"display": "inline-flex",
										"alignItems": "center",
										"justifyContent": "center",
										"width": "22px",
										"height": "22px",
										"borderRadius": "50%",
										"backgroundColor": f"{color}",
										"color": "#fff",
										"fontWeight": "700",
										"fontSize": "13px",
										"marginRight": "8px",
									},
								),
								html.Div([html.Span("Detection Results", style={"fontWeight": "600", "fontSize": "20px"})]),
							],
							style={"display": "flex", "alignItems": "center"},
						)
					)
				]
			),
			html.P("Detected Viral Genes in Dataset", style={"fontWeight": "600", "marginBottom": "10px"}),
			html.Div(
				[
					dbc.Row(
						[
							*[
								dbc.Col(
									virus_card(
										"\U0001f537",
										virus,
										f"{count} genes",
										html.Span(
											[
												f"Matched {count} of {all_genes_count} genes found ({((count / all_genes_count) * 100.0 if all_genes_count else 0.0):.1f}%)",
												html.Br(),
												f"{len(detected_features)} features detected",
											]
										),
									),
									xs=12,
									md=6,
									lg=4,
								)
								for virus, count in gene_count_per_virus.items()
							],
							dbc.Col(
								html.Div(
									[
										html.P(
											"Total unique viral genes detected",
											style={"fontSize": "12px", "color": "#6c757d", "margin": "0 0 2px"},
										),
										html.Span(
											f"{unique_count}",
											id="detected-viral-unique-count",
											style={"fontSize": "32px", "fontWeight": "700", "color": "#198754"},
										),
									],
									style={
										"padding": "12px 16px",
										"textAlign": "center",
										"border": "1px solid #dee2e6",
										"borderRadius": "6px",
										"backgroundColor": "#fff",
										"height": "100%",
									},
								),
								xs=12,
								md=6,
								lg=4,
							),
						],
						className="g-2",
					),
				],
				style={
					"border": "1px solid #dee2e6",
					"borderRadius": "6px",
					"backgroundColor": "#fff",
					"padding": "8px",
					"marginBottom": "12px",
				},
			),
			html.Div(
				[
					dbc.Accordion(
						[
							dbc.AccordionItem(
								[
									dcc.Dropdown(
										id="detected-viral-genes-dropdown",
										options=gene_entries,
										value=[],
										multi=True,
										placeholder="Select detected viral genes to remove",
									),
									dbc.Button(
										"Remove Selected",
										id="remove-detected-viral-genes-button",
										color="danger",
										size="sm",
										style={"marginTop": "8px", "marginBottom": "8px"},
									),
									html.Div(
										id="current-detected-viral-gene-list",
										children=html.P(
											grouped_list_text,
											style={"fontSize": "13px", "color": "#6c757d", "whiteSpace": "pre-wrap"},
										),
									),
								],
								title="Detected Genes" + (f", {len(not_found)} not found" if not_found else ""),
							)
						]
					)
				],
				style={"backgroundColor": "#fff", "border": "1px solid #dee2e6", "borderRadius": "8px", "padding": "16px", "marginBottom": "16px"},
			),
			dbc.Row(
				[
					dbc.Col(
						[
							html.Div(
								[
									html.Span("Add/Append Viral Genes", style={"fontWeight": "600"}),
									html.Span(" (Optional)", style={"color": "#6c757d"}),
									help_icon(),
								],
								style={"display": "flex", "alignItems": "center", "gap": "4px", "marginBottom": "4px"},
							),
							html.P("Add additional viral genes that may not be in our database.", style={"fontSize": "12px", "color": "#6c757d", "margin": "0"}),
						],
						width=3,
					),
					dbc.Col(dbc.Input(id="append-viral-genes-input", placeholder="Enter gene names separated by commas or new lines", type="text")),
					dbc.Col(
						dbc.Button("+ Add Genes", id="append-viral-genes-button", color="primary", outline=True, size="sm", style={"whiteSpace": "nowrap"}),
						width="auto",
						style={"display": "flex", "alignItems": "center"},
					),
				],
				align="center",
				style={"backgroundColor": "#fff", "border": "1px solid #dee2e6", "borderRadius": "8px", "padding": "12px 16px", "marginBottom": "16px"},
			),
			dbc.Alert(id="append-viral-genes-alert", children="", color="danger", is_open=False, dismissable=True),
		]
	)
	cache_results(**{"viral-gene-detection-results-container": result})
	return result


def register_vd_callbacks(app):
	"""Register all Dash callbacks for viral detection and curation flows.

	Use:
	- Connects panel controls to automatic/custom detection and post-detection edits.

	Interacts with:
	- find_viral_genes/find_custom_viral_genes, dataset store getters/updaters,
	  normalize_gene_name, and UI renderer helpers.

	Inputs:
	- app (dash.Dash): application instance.

	Outputs:
	- None (side effect: callback registration).
	"""
	@app.callback(
		Output("virus-select-dropdown", "options"),
		Input("virus-select-dropdown", "id"),
	)
	def populate_viral_set_options(_):
		"""Populate viral-set dropdown options from files available on disk.

		Inputs: ignored trigger value.
		Outputs: list[dict] dropdown options.
		Interacts with: list_viral_gene_sets.
		"""
		options = [{"label": "Auto-detect (search all known viruses)", "value": "__auto__"}]
		options.extend({"label": set_name, "value": set_name} for set_name in list_viral_gene_sets())
		return options

	@app.callback(
		Output("custom-detection-card", "style"),
		Output("custom-detection-card-title", "style"),
		Output("custom-gene-list-container", "hidden"),
		Input("detection-method-radio", "value"),
	)
	def update_custom_card_style(selected_method):
		"""Toggle custom detection card emphasis and custom list visibility.

		Inputs: selected_method (str).
		Outputs: custom card style, title style, textarea hidden flag.
		Interacts with: detection method radio UI.
		"""
		if selected_method == "custom":
			return {
				"border": "2px solid #0d6efd",
				"borderRadius": "8px",
				"padding": "12px 16px",
				"backgroundColor": "#f0f6ff",
				"cursor": "pointer",
			}, {"fontWeight": "600", "color": "#0d6efd"}, False
		return {
			"border": "2px solid #dee2e6",
			"borderRadius": "8px",
			"padding": "12px 16px",
			"backgroundColor": "#ffffff",
			"cursor": "pointer",
		}, {"fontWeight": "600", "color": "#212529"}, True

	@app.callback(
		Output("automatic-detection-card", "style"),
		Output("automatic-detection-card-title", "style"),
		Input("detection-method-radio", "value"),
	)
	def update_automatic_card_style(selected_method):
		"""Toggle automatic detection card emphasis.

		Inputs: selected_method (str).
		Outputs: automatic card style and title style.
		Interacts with: detection method radio UI.
		"""
		if selected_method == "automatic":
			return {
				"border": "2px solid #0d6efd",
				"borderRadius": "8px",
				"padding": "12px 16px",
				"backgroundColor": "#f0f6ff",
				"cursor": "pointer",
			}, {"fontWeight": "600", "color": "#0d6efd"}
		return {
			"border": "2px solid #dee2e6",
			"borderRadius": "8px",
			"padding": "12px 16px",
			"backgroundColor": "#ffffff",
			"cursor": "pointer",
		}, {"fontWeight": "600", "color": "#212529"}

	@app.callback(
		Output("viral-gene-detection-loading-signal", "children"),
		Output("viral-gene-detection-results-container", "children"),
		Input("run-viral-gene-detection-button", "n_clicks"),
		State("detection-method-radio", "value"),
		State("custom-gene-list-input", "value"),
		State("virus-select-dropdown", "value"),
	)
	@_viral_detection_tracker.track
	def run_viral_gene_detection(n_clicks, selected_method, custom_gene_list, selected_virus):
		"""Execute viral detection and return rendered results.

		Use:
		- Runs automatic or custom detection, aggregates metrics, persists state.

		Interacts with:
		- find_viral_genes/find_custom_viral_genes, update_state_store,
		  create_viral_gene_detection_results.

		Inputs:
		- n_clicks, selected_method, custom_gene_list, selected_virus.

		Outputs:
		- loading text and detection-results component.
		"""
		if n_clicks is None or n_clicks == 0:
			return no_update, no_update

		if selected_method == "automatic":
			with _viral_detection_tracker.phase("analysis"):
				detected = find_viral_genes(selected_virus)
			detected_by_virus = _normalize_automatic_detection_payload(detected, selected_virus)
			gene_count_per_virus = {}
			detected_gene_sets = {}
			all_genes = []
			detected_features = []
			matched_features_by_gene = {}
			for virus_name, value in detected_by_virus.items():
				if not isinstance(value, dict):
					continue

				features = value.get("features", [])
				genes = value.get("genes", [])
				if features == [] and genes == []:
					continue
				gene_count_per_virus[virus_name] = len(genes)
				detected_gene_sets[virus_name] = sorted(set(genes))
				all_genes.extend(genes)
				detected_features.extend(features)
				for gene, features in value.get("matched_features_by_gene", {}).items():
					matched_features_by_gene.setdefault(gene, set()).update(features)

			all_genes = sorted(set(all_genes))
			detected_features = sorted(set(detected_features))
			unique_count = len(all_genes)
			gene_entries = _build_viral_gene_entries(all_genes, matched_features_by_gene)
			update_state_store(
				viral_detection={
					"viral_genes": ", ".join(all_genes),
					"viral_features": ", ".join(detected_features),
					"detected_gene_sets": {
						set_name: sorted(set(genes))
						for set_name, genes in sorted(detected_gene_sets.items())
					},
					"matched_features_by_gene": {
						gene: sorted(features)
						for gene, features in sorted(matched_features_by_gene.items())
					},
				}
			)

			if unique_count > 0:
				with _viral_detection_tracker.phase("visualization"):
					result = create_viral_gene_detection_results("\u2713", "#198754", gene_count_per_virus, detected_features, gene_entries, unique_count, detected_gene_sets=detected_gene_sets)
				return (
					"Automatic detection completed",
					result,
				)
			with _viral_detection_tracker.phase("visualization"):
				result = create_viral_gene_detection_results("\u2717", "#dc3545", gene_count_per_virus, detected_features, [], 0, detected_gene_sets=detected_gene_sets)
			return (
				"Automatic detection completed, but no viral genes were detected.",
				result,
			)

		detected = find_custom_viral_genes(custom_gene_list)
		all_genes = sorted(set(detected["genes"]))
		detected_features = sorted(set(detected["features"]))
		matched_features_by_gene = {
			gene: set(features)
			for gene, features in detected.get("matched_features_by_gene", {}).items()
		}
		unique_count = len(all_genes)
		not_found = detected["not_found"] if detected["not_found"] else None
		detected_gene_sets = {"Custom List": all_genes}
		gene_entries = _build_viral_gene_entries(all_genes, matched_features_by_gene)
		update_state_store(
			viral_detection={
				"viral_genes": ", ".join(all_genes),
				"viral_features": ", ".join(detected_features),
				"detected_gene_sets": {"Custom List": all_genes},
				"matched_features_by_gene": {
					gene: sorted(features)
					for gene, features in sorted(matched_features_by_gene.items())
				},
			}
		)
		gene_count_per_virus = {"Custom List": unique_count}
		if unique_count > 0:
			with _viral_detection_tracker.phase("visualization"):
				result = create_viral_gene_detection_results("\u2713", "#198754", gene_count_per_virus, detected_features, gene_entries, unique_count, not_found, detected_gene_sets)
			return (
				"Custom detection completed",
				result,
			)
		with _viral_detection_tracker.phase("visualization"):
			result = create_viral_gene_detection_results("\u2717", "#dc3545", gene_count_per_virus, detected_features, [], 0, not_found, detected_gene_sets)
		return (
			"Custom detection completed, but no viral genes were detected.",
			result,
		)

	@app.callback(
		Output("detected-viral-genes-dropdown", "options"),
		Output("detected-viral-genes-dropdown", "value"),
		Output("current-detected-viral-gene-list", "children"),
		Output("detected-viral-unique-count", "children"),
		Output("append-viral-genes-alert", "children"),
		Output("append-viral-genes-alert", "is_open"),
		Input("remove-detected-viral-genes-button", "n_clicks"),
		Input("append-viral-genes-button", "n_clicks"),
		State("detected-viral-genes-dropdown", "value"),
		State("append-viral-genes-input", "value"),
		prevent_initial_call=True,
	)
	def curate_viral_gene_list(remove_clicks, append_clicks, selected_genes, appended_text):
		"""Apply interactive remove/add edits to detected viral genes.

		Use:
		- Removes selected detected genes or appends validated genes from user input.
		- Keeps viral genes/features/gene-feature mapping synchronized in state.

		Interacts with:
		- get_state_store/update_state_store, get_dataset/get_working_dataset,
		  normalize_gene_name, _split_input_genes, _build_viral_gene_entries.

		Inputs:
		- remove/add click counts, selected genes, appended text.

		Outputs:
		- updated dropdown options/value, list body, count, alert text/open state.
		"""
		history = get_state_store()
		viral_state = history.get("viral_detection", {})

		genes = set([g.strip() for g in viral_state.get("viral_genes", "").split(",") if g.strip()])
		features = set([f.strip() for f in viral_state.get("viral_features", "").split(",") if f.strip()])
		matched_features_by_gene = {
			gene: set(gene_features)
			for gene, gene_features in viral_state.get("matched_features_by_gene", {}).items()
		}
		detected_gene_sets = {
			set_name: set(set_genes)
			for set_name, set_genes in viral_state.get("detected_gene_sets", {}).items()
		}

		if not detected_gene_sets and genes:
			detected_gene_sets = {"Detected": set(genes)}

		alert_text = ""
		alert_open = False
		triggered = ctx.triggered_id

		if triggered == "remove-detected-viral-genes-button":
			for gene in (selected_genes or []):
				genes.discard(gene)
				for feature in matched_features_by_gene.get(gene, set()):
					features.discard(feature)
				matched_features_by_gene.pop(gene, None)
				for set_name in list(detected_gene_sets.keys()):
					detected_gene_sets[set_name].discard(gene)

		elif triggered == "append-viral-genes-button":
			adata = get_dataset() or get_working_dataset()
			if adata is None:
				alert_text = "No dataset available. Please upload a dataset before appending genes."
				alert_open = True
			else:
				query_genes = _split_input_genes(appended_text)
				not_found = []
				for query in query_genes:
					matches = []
					for feature in adata.var_names:
						matched_gene = normalize_gene_name(feature, [query])
						if matched_gene == query:
							matches.append(feature)

					if matches:
						gene_name = query
						genes.add(gene_name)
						matched_features_by_gene.setdefault(gene_name, set()).update(matches)
						features.update(matches)
						detected_gene_sets.setdefault("Manual Append", set()).add(gene_name)
					else:
						not_found.append(query)

				if not_found:
					alert_text = (
						"Could not match these genes in the dataset: "
						+ ", ".join(not_found)
						+ ". Please check spelling and case."
					)
					alert_open = True

		detected_gene_sets = {
			gene: set(gene_features)
			for gene, gene_features in detected_gene_sets.items()
		}
		detected_gene_sets = {
			set_name: sorted(set_genes)
			for set_name, set_genes in sorted(detected_gene_sets.items())
			if set_genes
		}

		update_state_store(
			viral_detection={
				"viral_genes": ", ".join(sorted(genes)),
				"viral_features": ", ".join(sorted(features)),
				"detected_gene_sets": detected_gene_sets,
				"matched_features_by_gene": {
					gene: sorted(gene_features)
					for gene, gene_features in sorted(matched_features_by_gene.items())
				},
			}
		)

		entries = _build_viral_gene_entries(genes, matched_features_by_gene)
		list_text = _format_detected_grouped_list(detected_gene_sets, "No viral genes currently selected.")

		return (
			entries,
			[],
			html.P(list_text, style={"fontSize": "13px", "color": "#6c757d", "whiteSpace": "pre-wrap"}),
			str(len(genes)),
			alert_text,
			alert_open,
		)
