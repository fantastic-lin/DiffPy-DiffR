# DiffPy datasets

All matrices are compressed SciPy CSR resources with explicit row and column
identifiers. Annotation vectors load as pandas Series; categorical annotations
retain their declared category order. PPI networks and regulons are packaged
with DiffPy; the three example expression cohorts are a separate download.

## External test data

Extract ``DiffPy_test_data.zip`` beneath one directory, then set
``DIFFPY_DATA_DIR`` to that directory. You can instead pass the directory as
``data_dir`` to ``load_dataset``.

```python
from diffpy import load_dataset

chu = load_dataset("chu", data_dir="/path/to/diffpy-data")
```

Annotation vectors correspond to expression columns by position. Their Series
indices are not guaranteed to contain expression cell identifiers: Chu and
Stomach annotations use positional indices, while the Liver stage-code index
contains stage labels. Preserve the original order rather than joining an
annotation to expression columns by Series index.

| Dataset | Objects | Description |
|---|---|---|
| `chu` | `expression`, `phenotype`, `signaling_entropy` | Separate test data: embryonic stem and endothelial progenitor cells |
| `liver` | `expression`, `stage_codes`, `plot_colors` | Separate test data: developing liver time course |
| `stomach` | `expression`, `differentiation_state` | Separate test data: two stomach epithelial differentiation states |
| `ppi_2012` | `adjacency` | Bundled 2012 PPI network |
| `ppi_2016` | `adjacency` | Bundled 2016 PPI network |
| `TF-regulon_network_breast` | `regulon` | Breast epithelial regulon |
| `TF-regulon_network_colon` | `regulon` | Colon epithelial regulon |
| `TF-regulon_network_esophagus` | `regulon` | Esophageal epithelial regulon |
| `TF-regulon_network_stomach` | `regulon` | Stomach epithelial regulon |
| `TF-regulon_network_kidney` | `regulon` | Kidney epithelial regulon |
| `TF-regulon_network_liver` | `regulon` | Liver epithelial regulon |
| `TF-regulon_network_lung` | `regulon` | Lung epithelial regulon |
| `TF-regulon_network_pancreas` | `regulon` | Pancreatic epithelial regulon |
| `TF-regulon_network_skin` | `regulon` | Skin epithelial regulon |

## Expression preprocessing

The test expression cohorts have different preprocessing histories. Do not
apply another log transformation to the liver or stomach matrices.

| Dataset | Library-size normalized | Log-transformed | Transformation |
|---|---|---|---|
| Chu | Yes | No | None |
| Liver | Yes | Yes | `log2(x + 1.1)` |
| Stomach | Yes, scale factor 10,000 | Yes | `log2(x + 1)` |

`integrate_expression_network()` always applies `log2(x + 1.1)` and therefore
requires unlogged, library-size-normalized expression. The test-data Chu matrix
meets this requirement and can be supplied directly. The test-data Liver and
Stomach matrices are already logged and should not be supplied directly to
this function; use corresponding unlogged, library-size-normalized matrices
instead. Test-data Liver remains suitable for the demonstrated CCAT and
trajectory workflow, and test-data Stomach remains suitable for the
regulatory-activity workflow.
