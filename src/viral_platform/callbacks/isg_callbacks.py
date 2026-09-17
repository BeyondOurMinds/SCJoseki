import re

from dash import Input, Output, State, ctx, dcc, html, no_update
import dash_bootstrap_components as dbc
import scanpy as sc

from viral_platform.analysis.isg_analysis import (
	find_custom_isg_genes,
	find_isg_genes,
	list_isg_sets,
	load_isg_set,
	normalize_gene_name,
)
from viral_platform.state.dataset_store import (
	cache_results,
	get_dataset,
	get_state_store,
	get_working_dataset,
	set_working_dataset,
	sync_state_with_dataset,
	update_state_store,
)
from viral_platform.utils.runtime_tracker import get_tracker

_isg_detection_tracker = get_tracker("ISG: Detection")
_isg_summary_tracker = get_tracker("ISG: Summary Stats")


def help_icon():
	"""Create the shared help icon element used across ISG panel rows.

	Use:
	- Provides consistent inline helper icon styling in result/action sections.

	Interacts with:
	- create_isg_detection_results.

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


def isg_set_card(icon_char, name, badge_text, detail):
	"""Render one ISG set summary card in the detection results header.

	Use:
	- Displays per-set match counts and descriptive text.

	Interacts with:
	- _build_isg_set_cards, create_isg_detection_results.

	Inputs:
	- icon_char (str), name (str), badge_text (str), detail (Dash child).

	Outputs:
	- dash.html.Div card element.
	"""
	return html.Div(
		[
			html.Div(
				[
					html.Span(
						icon_char,
						style={
							"fontSize": "28px",
							"marginRight": "12px",
						},
					),
					html.Div(
						[
							html.Div(
								[
									html.A(
										name,
										href="#",
										style={
											"fontWeight": "600",
											"color": "#0d6efd",
											"textDecoration": "none",
										},
									),
									dbc.Badge(
										badge_text,
										color="secondary",
										style={"marginLeft": "8px", "fontSize": "11px"},
									),
								],
								style={"display": "flex", "alignItems": "center"},
							),
							html.P(
								detail,
								style={"margin": "2px 0 0", "fontSize": "13px", "color": "#6c757d"},
							),
						]
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


def _build_isg_set_cards(gene_count_per_set, total_genes_per_set, detected_features):
	"""Build result cards for all detected ISG sets.

	Use:
	- Converts aggregate counts into display cards and percentages.

	Interacts with:
	- isg_set_card, create_isg_detection_results.

	Inputs:
	- gene_count_per_set (dict[str, int])
	- total_genes_per_set (dict[str, int])
	- detected_features (list[str])

	Outputs:
	- list[dash.html.Div]: rendered set cards.
	"""
	cards = []
	for set_name, count in gene_count_per_set.items():
		total_count = total_genes_per_set.get(set_name, 0)
		pct = (count / total_count * 100.0) if total_count else 0.0
		cards.append(
			isg_set_card(
				"\U0001f7e6",
				set_name,
				f"{count} genes",
				html.Span(
					[
						f"Matched {count} of {total_count} detected ISGs ({pct:.1f}%)",
						html.Br(),
						f"{len(detected_features)} features detected",
					]
				),
			)
		)
	return cards


def _split_input_genes(raw_text):
	"""Split add-gene text input into cleaned gene tokens.

	Use:
	- Normalizes comma/newline user input before append validation.

	Interacts with:
	- curate_isg_gene_list.

	Inputs:
	- raw_text (str|None): user-entered text.

	Outputs:
	- list[str]: trimmed non-empty gene tokens.
	"""
	if not raw_text:
		return []
	tokens = re.split(r"[,\n\r]+", raw_text)
	return [token.strip() for token in tokens if token.strip()]


def _build_isg_gene_entries(genes, matched_features_by_gene):
	"""Create dropdown entries for detected ISGs with optional alias context.

	Use:
	- Formats curation dropdown/list labels and includes matched-as only for alias cases.

	Interacts with:
	- create_isg_detection_results, curate_isg_gene_list.

	Inputs:
	- genes (iterable[str])
	- matched_features_by_gene (dict[str, iterable[str]])

	Outputs:
	- list[dict]: dropdown options with label/value fields.
	"""
	entries = []
	for gene in sorted(set(genes)):
		matched = sorted(matched_features_by_gene.get(gene, []))
		alias_matches = [
			feature
			for feature in matched
			if normalize_gene_name(feature) != gene
		]
		label = f"{gene} (matched as {', '.join(alias_matches)})" if alias_matches else gene
		entries.append({"label": label, "value": gene})
	return entries


def _format_detected_grouped_list(set_to_genes, empty_message):
	"""Format grouped detected-gene text as '(set): gene1, gene2' chunks.

	Use:
	- Produces display-only grouped list text for detection result accordions.

	Interacts with:
	- create_isg_detection_results, curate_isg_gene_list.

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


def _normalize_automatic_detection_payload(detected, selected_set):
	"""Return a per-set mapping regardless of the finder return shape."""
	if not isinstance(detected, dict):
		return {}

	if "features" in detected and "genes" in detected:
		set_key = selected_set or "Selected ISG Set"
		return {set_key: detected}

	return detected


def create_isg_detection_results(
	pass_fail,
	color,
	gene_count_per_set,
	total_genes_per_set,
	detected_features,
	detected_gene_entries,
	unique_count,
	not_found=None,
	detected_gene_sets=None,
):
	"""Render the full ISG detection results section.

	Use:
	- Produces summary cards, editable detected list, add/remove controls, and alerts.

	Interacts with:
	- _build_isg_set_cards, _build_isg_gene_entries, curate_isg_gene_list callback.

	Inputs:
	- pass_fail, color: status display values.
	- gene_count_per_set/total_genes_per_set/detected_features: aggregate metrics.
	- detected_gene_entries: dropdown/list items.
	- unique_count (int), not_found (list[str]|None).

	Outputs:
	- dash.html.Div containing the interactive results UI.
	"""
	grouped_list_text = _format_detected_grouped_list(
		detected_gene_sets,
		"No ISGs currently selected.",
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
			html.P("Detected ISGs in Dataset", style={"fontWeight": "600", "marginBottom": "10px"}),
			html.Div(
				[
					dbc.Row(
						[
							*[
								dbc.Col(card, xs=12, md=6, lg=4)
								for card in _build_isg_set_cards(gene_count_per_set, total_genes_per_set, detected_features)
							],
							dbc.Col(
								html.Div(
									[
										html.P(
											"Total unique ISGs detected",
											style={"fontSize": "12px", "color": "#6c757d", "margin": "0 0 2px"},
										),
										html.Span(
											f"{unique_count}",
											id="detected-isg-unique-count",
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
										id="detected-isg-genes-dropdown",
										options=detected_gene_entries,
										value=[],
										multi=True,
										placeholder="Select detected ISGs to remove",
									),
									dbc.Button(
										"Remove Selected",
										id="remove-detected-isg-genes-button",
										color="danger",
										size="sm",
										style={"marginTop": "8px", "marginBottom": "8px"},
									),
									html.Div(
										id="current-detected-isg-gene-list",
										children=html.P(
											grouped_list_text,
											style={"fontSize": "13px", "color": "#6c757d", "whiteSpace": "pre-wrap"},
										),
									),
								],
								title="Detected ISGs" + (f", {len(not_found)} not found" if not_found else ""),
							)
						]
					)
				],
				style={
					"backgroundColor": "#fff",
					"border": "1px solid #dee2e6",
					"borderRadius": "8px",
					"padding": "16px",
					"marginBottom": "16px",
				},
			),
			dbc.Row(
				[
					dbc.Col(
						[
							html.Div(
								[
									html.Span("Add/Append ISGs", style={"fontWeight": "600"}),
									html.Span(" (Optional)", style={"color": "#6c757d"}),
									help_icon(),
								],
								style={
									"display": "flex",
									"alignItems": "center",
									"gap": "4px",
									"marginBottom": "4px",
								},
							),
							html.P(
								"Add additional ISGs that may not be in the curated set.",
								style={"fontSize": "12px", "color": "#6c757d", "margin": "0"},
							),
						],
						width=3,
					),
					dbc.Col(
						dbc.Input(
							id="append-isg-genes-input",
							placeholder="Enter gene names separated by commas or new lines",
							type="text",
						)
					),
					dbc.Col(
						dbc.Button(
							"+ Add Genes",
							id="append-isg-genes-button",
							color="primary",
							outline=True,
							size="sm",
							style={"whiteSpace": "nowrap"},
						),
						width="auto",
						style={"display": "flex", "alignItems": "center"},
					),
				],
				align="center",
				style={
					"backgroundColor": "#fff",
					"border": "1px solid #dee2e6",
					"borderRadius": "8px",
					"padding": "12px 16px",
					"marginBottom": "16px",
				},
			),
			dbc.Alert(
				id="append-isg-genes-alert",
				children="",
				color="danger",
				is_open=False,
				dismissable=True,
			),
		]
	)
	cache_results(**{"isg-detection-results-container": result})
	return result


def isg_summary_results(adata):
	"""Render summary-statistics table for computed ISG scores.

	Use:
	- Displays high-level distribution metrics for adata.obs['ISG_score'].

	Interacts with:
	- run_isg_summary_stats callback.

	Inputs:
	- adata (AnnData): dataset with ISG_score column in obs.

	Outputs:
	- dash_bootstrap_components.Card with summary table.
	"""
	isg_score = adata.obs["ISG_score"]
	positive_cells = (isg_score > 0).sum()

	return dbc.Card(
		dbc.CardBody(
			[
				html.H5("ISG Summary Statistics"),
				dbc.Table(
					[
						html.Tbody(
							[
								html.Tr([html.Td("Cells with Positive ISG Score"), html.Td(f"{positive_cells}/{adata.n_obs}")]),
								html.Tr([html.Td("Maximum ISG Score"), html.Td(f"{isg_score.max():.4f}")]),
								html.Tr([html.Td("Average ISG Score"), html.Td(f"{isg_score.mean():.4f}")]),
								html.Tr([html.Td("Median ISG Score"), html.Td(f"{isg_score.median():.4f}")]),
								html.Tr([html.Td("Std. Dev. ISG Score"), html.Td(f"{isg_score.std():.4f}")]),
							]
						)
					]
				),
			]
		)
	)


def register_isg_callbacks(app):
	"""Register all Dash callbacks for ISG detection, curation, and summary analysis.

	Use:
	- Wires UI controls to analysis functions and state-store synchronization.

	Interacts with:
	- find_isg_genes, find_custom_isg_genes, load_isg_set, score_genes,
	  and dataset_store state helpers.

	Inputs:
	- app (dash.Dash): application instance.

	Outputs:
	- None (side effect: callback registration).
	"""
	@app.callback(
		Output("isg-set-select-dropdown", "options"),
		Input("isg-set-select-dropdown", "id"),
	)
	def populate_isg_set_options(_):
		"""Populate ISG set dropdown options from files available on disk.

		Inputs: ignored trigger value.
		Outputs: list[dict] dropdown options.
		Interacts with: list_isg_sets.
		"""
		options = [{"label": "Auto-detect (search all ISG sets)", "value": "__auto__"}]
		options.extend({"label": set_name, "value": set_name} for set_name in list_isg_sets())
		return options

	@app.callback(
		Output("custom-isg-detection-card", "style"),
		Output("custom-isg-detection-card-title", "style"),
		Output("custom-isg-gene-list-container", "hidden"),
		Input("isg-detection-method-radio", "value"),
	)
	def update_custom_isg_card_style(selected_method):
		"""Toggle visual emphasis and visibility for custom detection controls.

		Inputs: selected_method (str).
		Outputs: custom card style, title style, textarea hidden flag.
		Interacts with: ISG detection method radio UI.
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
		Output("automatic-isg-detection-card", "style"),
		Output("automatic-isg-detection-card-title", "style"),
		Input("isg-detection-method-radio", "value"),
	)
	def update_automatic_isg_card_style(selected_method):
		"""Toggle visual emphasis for automatic detection card.

		Inputs: selected_method (str).
		Outputs: automatic card style and title style.
		Interacts with: ISG detection method radio UI.
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
		Output("isg-detection-loading-signal", "children"),
		Output("isg-detection-results-container", "children"),
		Input("run-isg-detection-button", "n_clicks"),
		State("isg-detection-method-radio", "value"),
		State("custom-isg-gene-list-input", "value"),
		State("isg-set-select-dropdown", "value"),
	)
	@_isg_detection_tracker.track
	def run_isg_detection(n_clicks, selected_method, custom_gene_list, selected_set):
		"""Execute ISG detection and render results.

		Use:
		- Runs automatic or custom detection, aggregates counts, and updates shared state.

		Interacts with:
		- find_isg_genes/find_custom_isg_genes, update_state_store,
		  create_isg_detection_results.

		Inputs:
		- n_clicks, selected_method, custom_gene_list, selected_set.

		Outputs:
		- loading text, results component.
		"""
		if n_clicks is None or n_clicks == 0:
			return no_update, no_update

		if selected_method == "automatic":
			with _isg_detection_tracker.phase("analysis"):
				detected = find_isg_genes(selected_set)
			detected_by_set = _normalize_automatic_detection_payload(detected, selected_set)
			gene_count_per_set = {}
			total_genes_per_set = {}
			detected_gene_sets = {}
			all_genes = []
			detected_features = []
			core_to_features = {}

			for key, value in detected_by_set.items():
				if not isinstance(value, dict):
					continue
				features = value.get("features", [])
				genes = value.get("genes", [])
				if features == [] and genes == []:
					continue
				gene_count_per_set[key] = len(genes)
				detected_gene_sets[key] = sorted(set(genes))
				all_genes.extend(genes)
				detected_features.extend(features)
				for core_gene, features in value.get("matched_features_by_gene", {}).items():
					core_to_features.setdefault(core_gene, set()).update(features)

			all_genes = sorted(set(all_genes))
			detected_features = sorted(set(detected_features))
			unique_count = len(all_genes)
			total_genes_per_set = {name: unique_count for name in gene_count_per_set}
			detected_gene_entries = _build_isg_gene_entries(all_genes, core_to_features)

			update_state_store(
				isg_detection={
					"isg_genes": ", ".join(all_genes),
					"isg_features": ", ".join(detected_features),
					"detected_gene_sets": {
						set_name: sorted(set(genes))
						for set_name, genes in sorted(detected_gene_sets.items())
					},
					"matched_features_by_gene": {
						gene: sorted(features)
						for gene, features in sorted(core_to_features.items())
					},
				}
			)

			if unique_count > 0:
				with _isg_detection_tracker.phase("visualization"):
					result = create_isg_detection_results(
						"\u2713",
						"#198754",
						gene_count_per_set,
						total_genes_per_set,
						detected_features,
						detected_gene_entries,
						unique_count,
						detected_gene_sets=detected_gene_sets,
					)
				return (
					"Automatic ISG detection completed",
					result,
				)

			with _isg_detection_tracker.phase("visualization"):
				result = create_isg_detection_results(
					"\u2717",
					"#dc3545",
					gene_count_per_set,
					total_genes_per_set,
					detected_features,
					[],
					0,
					detected_gene_sets=detected_gene_sets,
				)
			return (
				"Automatic ISG detection completed, but no ISGs were detected.",
				result,
			)

		if not custom_gene_list:
			return (
				"Custom ISG detection requires a gene list.",
				html.P("Provide a custom ISG list and run detection again.", style={"color": "#6c757d", "fontSize": "14px"}),
			)

		detected = find_custom_isg_genes(custom_gene_list)
		all_genes = sorted(set(detected["genes"]))
		detected_features = sorted(set(detected["features"]))
		unique_count = len(all_genes)
		core_to_features = {
			gene: set(features)
			for gene, features in detected.get("matched_features_by_gene", {}).items()
		}
		detected_gene_sets = {"Custom List": all_genes}
		detected_gene_entries = _build_isg_gene_entries(all_genes, core_to_features)
		gene_count_per_set = {"Custom List": unique_count}
		total_genes_per_set = {"Custom List": unique_count}
		not_found = detected["not_found"] if detected["not_found"] else None

		update_state_store(
			isg_detection={
				"isg_genes": ", ".join(all_genes),
				"isg_features": ", ".join(detected_features),
				"detected_gene_sets": {"Custom List": all_genes},
				"matched_features_by_gene": {
					gene: sorted(features)
					for gene, features in sorted(core_to_features.items())
				},
			}
		)

		if unique_count > 0:
			with _isg_detection_tracker.phase("visualization"):
				result = create_isg_detection_results(
					"\u2713",
					"#198754",
					gene_count_per_set,
					total_genes_per_set,
					detected_features,
					detected_gene_entries,
					unique_count,
					not_found,
					detected_gene_sets,
				)
			return (
				"Custom ISG detection completed",
				result,
			)

		with _isg_detection_tracker.phase("visualization"):
			result = create_isg_detection_results(
				"\u2717",
				"#dc3545",
				gene_count_per_set,
				total_genes_per_set,
				detected_features,
				[],
				0,
				not_found,
				detected_gene_sets,
			)
		return (
			"Custom ISG detection completed, but no ISGs were detected.",
			result,
		)

	@app.callback(
		Output("detected-isg-genes-dropdown", "options"),
		Output("detected-isg-genes-dropdown", "value"),
		Output("current-detected-isg-gene-list", "children"),
		Output("detected-isg-unique-count", "children"),
		Output("append-isg-genes-alert", "children"),
		Output("append-isg-genes-alert", "is_open"),
		Input("remove-detected-isg-genes-button", "n_clicks"),
		Input("append-isg-genes-button", "n_clicks"),
		State("detected-isg-genes-dropdown", "value"),
		State("append-isg-genes-input", "value"),
		State("isg-set-select-dropdown", "value"),
		prevent_initial_call=True,
	)
	def curate_isg_gene_list(remove_clicks, append_clicks, selected_genes, appended_text, selected_set):
		"""Apply interactive remove/add edits to detected ISG list.

		Use:
		- Removes selected genes or appends validated genes that match dataset features.
		- Keeps genes/features/mappings synchronized in state.

		Interacts with:
		- get_state_store/update_state_store, get_dataset/get_working_dataset,
		  load_isg_set, normalize_gene_name, _build_isg_gene_entries.

		Inputs:
		- remove/add click counts, selected genes, appended text, selected set.

		Outputs:
		- updated dropdown options/value, list display, count text, alert text/open state.
		"""
		history = get_state_store()
		isg_state = history.get("isg_detection", {})

		genes = set([g.strip() for g in isg_state.get("isg_genes", "").split(",") if g.strip()])
		features = set([f.strip() for f in isg_state.get("isg_features", "").split(",") if f.strip()])
		matched_features_by_gene = {
			gene: set(gene_features)
			for gene, gene_features in isg_state.get("matched_features_by_gene", {}).items()
		}
		detected_gene_sets = {
			set_name: set(set_genes)
			for set_name, set_genes in isg_state.get("detected_gene_sets", {}).items()
		}

		if not detected_gene_sets and genes:
			detected_gene_sets = {"Detected": set(genes)}

		alert_text = ""
		alert_open = False
		triggered = ctx.triggered_id

		if triggered == "remove-detected-isg-genes-button":
			for gene in (selected_genes or []):
				genes.discard(gene)
				for feature in matched_features_by_gene.get(gene, set()):
					features.discard(feature)
				matched_features_by_gene.pop(gene, None)
				for set_name in list(detected_gene_sets.keys()):
					detected_gene_sets[set_name].discard(gene)

		elif triggered == "append-isg-genes-button":
			adata = get_dataset() or get_working_dataset()
			if adata is None:
				alert_text = "No dataset available. Please upload a dataset before appending genes."
				alert_open = True
			else:
				set_names = list_isg_sets() if selected_set in (None, "__auto__") else [selected_set]
				alias_to_core = {}
				set_to_core_genes = {}
				for set_name in set_names:
					set_data = load_isg_set(set_name)
					alias_to_core.update(set_data["alias_to_core"])
					set_to_core_genes[set_name] = set(set_data["core_genes"])

				not_found = []
				for query in _split_input_genes(appended_text):
					query_norm = normalize_gene_name(query)
					core_target = alias_to_core.get(query_norm, query_norm)

					matched_for_query = []
					for feature in adata.var_names:
						feature_norm = normalize_gene_name(feature)
						feature_core = alias_to_core.get(feature_norm, feature_norm)
						if feature_norm == query_norm or feature_core == core_target:
							matched_for_query.append(feature)

					if matched_for_query:
						genes.add(core_target)
						matched_features_by_gene.setdefault(core_target, set()).update(matched_for_query)
						features.update(matched_for_query)
						matched_sets = [
							set_name
							for set_name, core_genes in set_to_core_genes.items()
							if core_target in core_genes
						]
						if matched_sets:
							for set_name in matched_sets:
								detected_gene_sets.setdefault(set_name, set()).add(core_target)
						else:
							detected_gene_sets.setdefault("Manual Append", set()).add(core_target)
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
			set_name: sorted(set_genes)
			for set_name, set_genes in sorted(detected_gene_sets.items())
			if set_genes
		}

		update_state_store(
			isg_detection={
				"isg_genes": ", ".join(sorted(genes)),
				"isg_features": ", ".join(sorted(features)),
				"detected_gene_sets": detected_gene_sets,
				"matched_features_by_gene": {
					gene: sorted(gene_features)
					for gene, gene_features in sorted(matched_features_by_gene.items())
				},
			}
		)

		entries = _build_isg_gene_entries(genes, matched_features_by_gene)
		list_text = _format_detected_grouped_list(detected_gene_sets, "No ISGs currently selected.")

		return (
			entries,
			[],
			html.P(list_text, style={"fontSize": "13px", "color": "#6c757d", "whiteSpace": "pre-wrap"}),
			str(len(genes)),
			alert_text,
			alert_open,
		)

	@app.callback(
		Output("isg-summary-loading-signal", "children"),
		Output("isg-summary-results-container", "children"),
		Input("run-isg-summary-button", "n_clicks"),
		prevent_initial_call=True,
	)
	@_isg_summary_tracker.track
	def run_isg_summary_stats(n_clicks):
		"""Compute and display ISG score summary statistics.

		Use:
		- Scores cells using currently curated ISG features and returns summary table.

		Interacts with:
		- state store ISG feature values, sc.tl.score_genes,
		  set_working_dataset, sync_state_with_dataset, isg_summary_results.

		Inputs:
		- n_clicks (int).

		Outputs:
		- loading text and summary results component/message.
		"""
		if n_clicks is None or n_clicks == 0:
			return no_update, no_update

		history = get_state_store()
		isg_features = history.get("isg_detection", {}).get("isg_features", "")
		if not isg_features:
			return "done", "No ISG features detected. Please run ISG detection first."

		if isinstance(isg_features, str):
			isg_genes = [f.strip() for f in isg_features.split(",") if f.strip()]
		else:
			isg_genes = [f for f in isg_features if f]

		if not isg_genes:
			return "done", "No valid ISG features found. Please run ISG detection first."

		adata = get_working_dataset()
		if adata is None:
			adata = get_dataset()
		if adata is None:
			return "done", "No dataset available for ISG summary statistics."

		with _isg_summary_tracker.phase("analysis"):
			sc.tl.score_genes(
				adata,
				gene_list=isg_genes,
				score_name="ISG_score",
			)
			adata.obs["ISG_score"] = adata.obs["ISG_score"]

		set_working_dataset(adata)
		sync_state_with_dataset(adata)

		with _isg_summary_tracker.phase("processing"):
			result = isg_summary_results(adata)
		with _isg_summary_tracker.phase("visualization"):
			cache_results(**{"isg-summary-results-container": result})
		return "done", result
