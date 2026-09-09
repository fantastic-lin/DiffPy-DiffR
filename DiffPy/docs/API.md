# API reference

The definitive parameter and return documentation lives in each function's
NumPy-style docstring. In an interactive session, run `help(function_name)`.

## Matrix container

### `LabeledMatrix(values, row_names, col_names)`

Stores a dense NumPy or sparse SciPy matrix together with unique string labels.
It provides `toarray()`, `to_dataframe()`, `subset_rows()`, and
`subset_columns()` while preserving sparse storage where possible.

### `as_labeled_matrix(matrix, row_names=None, col_names=None)`

Normalizes pandas, NumPy, SciPy, and existing labeled inputs to one validated
representation.

## Network integration

### `integrate_expression_network(expression, network, ...)`

Requires nonnegative, library-size-normalized, unlogged expression; always
applies `log2(x + 1.1)`; matches genes; checks identifier overlap; extracts the
largest connected network component; and returns
`IntegrationResult(expression, adjacency)`. Sparse expression becomes dense
because `log2(1.1)` is nonzero. A maximum below 100 emits a warning that the
input may already be log-transformed, but the transformation is still applied.

## Potency estimates

### `compute_signaling_entropy(integrated, ...)`

Returns `SignalingEntropyResult` with:

- `signaling_entropy`
- `stationary_distribution`
- `local_entropy`
- `normalized_local_entropy`
- `maximum_entropy`

### `compute_ccat(expression, network, ...)`

Returns one connectome-transcriptome Pearson correlation per cell.

### `infer_potency_states(potency, ...)`

Returns `PotencyStateResult` with ordered state assignments, optional phenotype
summaries, transformed values, and the selected Gaussian mixture model.

## Trajectory

### `infer_diffusion_root(potency, expression, ...)`

Returns `TrajectoryResult` with:

- `diffusion_coordinates`
- `transition_matrix`
- `root_index` (zero-based)
- `root_cell`
- `potency_colors`
- `metadata`

### `compute_diffusion_pseudotime(trajectory, ...)`

Calculates destiny-style DPT distances from the inferred root, automatically
identifies hierarchical branches and three top-level tips, and returns
Gaussian-smoothed paths through diffusion-coordinate space. The result is a
`DiffusionPseudotimeResult` containing:

- `pseudotime`
- `distance_matrix`
- `branches`
- `tips`
- `tip_indices` and `tip_cells`
- `paths`
- `root_index` and `root_cell`

The default `paths_to=(1, 2, 3)` corresponds to the R vignette's
`paths_to=1:3`.

## Regulatory activity

### `estimate_regulatory_activity(expression, tissue, ...)`

Returns a transcription-factor-by-cell DataFrame, or a transcription-factor
Series for a single named expression profile. When the maximum expression
value exceeds 100, the function emits a warning and applies the sparse-safe
transformation `log2(x + 1)` before normalization. The required `tissue`
argument selects one of the nine packaged regulons; custom regulon input is not
supported.

The returned DataFrame has transcription factors in rows and cells in columns.
Calculate the average TF activity for every cell by averaging across rows:

```python
activity = estimate_regulatory_activity(
    expression,
    tissue="stomach",
    normalization="zscore",
)
average_tf_activity = activity.mean(axis=0)
```

## Dataset access

- `available_datasets()`
- `load_dataset(name, *, data_dir=None)`
- `load_ppi(version, *, data_dir=None)`
- `load_regulon(tissue)`

The Chu, liver, and stomach cohorts are distributed in
`DiffPy_test_data.zip`. Extract it and set `DIFFPY_DATA_DIR` to the extraction
directory, or pass that directory through `data_dir`. Both PPI networks and
all nine tissue regulons remain bundled with DiffPy.
