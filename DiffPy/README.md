# DiffPy

DiffPy is a Python toolkit for estimating differentiation potency,
quantifying signaling entropy, constructing potency-guided diffusion
trajectories, and inferring tissue-specific transcription-factor regulatory activity in epithelial cells from single-cell RNA-sequencing data.

The package is designed for labeled gene-by-cell matrices and supports both
dense NumPy arrays and SciPy sparse matrices. All four PPI networks and nine tissue
regulons are included. The larger example expression cohorts are distributed
separately in ``DiffPy_test_data.zip``.

## Installation

Install the package from its source directory or a built distribution:

```bash
python -m pip install .
```

Install optional tools for regenerating documentation figures:

```bash
python -m pip install ".[docs]"
```

Python 3.10 or newer is required.

### Install the separate test data

Extract the test-data archive into a data directory. It already contains the
``chu``, ``liver``, and ``stomach`` top-level folders:

```bash
mkdir -p diffpy-data
unzip DiffPy_test_data.zip -d diffpy-data
export DIFFPY_DATA_DIR="$PWD/diffpy-data"
```

The environment variable may contain multiple directories separated by the
platform path separator. Alternatively, pass the parent directory explicitly:

```python
from diffpy import load_dataset

chu = load_dataset("chu", data_dir="/path/to/diffpy-data")
```

## Analysis functions

| Function | Purpose |
|---|---|
| `DoIntegPPI` | Apply `log2(x + 1.1)`, match expression genes to a PPI network, and retain the largest connected component |
| `CompSRana` | Calculate normalized global and gene-local signaling entropy |
| `CompCCAT` | Estimate potency rapidly from expression/connectivity correlation |
| `InferPotencyStates` | Fit BIC-selected Gaussian mixtures and order discrete potency states |
| `InferDMAPandRoot` | Construct adaptive diffusion coordinates and infer a high-potency root cell |
| `compute_diffusion_pseudotime` | Calculate root-based DPT, branches, tips, and smoothed lineage paths |
| `SciraEstRegAct` | Estimate tissue-specific transcription-factor activity from regulon targets |

Every public function has a detailed NumPy-style docstring:

```python
from diffpy import CompCCAT

help(CompCCAT)
```

## Matrix conventions

Expression matrices always have genes in rows and cells in columns. A pandas
DataFrame supplies identifiers through its index and columns. NumPy and SciPy
matrices must be wrapped in `LabeledMatrix` or accompanied by explicit names.

```python
from scipy import sparse
from diffpy import LabeledMatrix

expression = LabeledMatrix(
    sparse.csr_matrix([[1.2, 0.0], [0.4, 2.1]]),
    row_names=("gene_a", "gene_b"),
    col_names=("cell_1", "cell_2"),
)
```

`DoIntegPPI()` requires nonnegative, library-size-normalized,
unlogged expression and always applies `log2(x + 1.1)`. The transformed matrix
is dense because every sparse structural zero becomes `log2(1.1)`, which is
nonzero. If the input maximum is below 100, the function emits a warning that
the matrix may already be log-transformed; this threshold is only a heuristic.

The test-data Chu expression matrix meets this input requirement and can be
supplied directly. The test-data Liver and Stomach matrices are already
log-transformed and should not be passed directly to
`DoIntegPPI()`; use corresponding unlogged,
library-size-normalized matrices instead. Liver remains suitable for the CCAT
and trajectory examples, and Stomach remains suitable for the
regulatory-activity example.

CCAT and regulatory-activity estimation use a separate sparse-safe
`log2(x + 1)` fallback when their input maximum exceeds 100. That transformation
preserves sparse zeros and emits a `RuntimeWarning` when applied.

## Potency workflow

```python
from diffpy import (
    CompCCAT,
    CompSRana,
    InferPotencyStates,
    DoIntegPPI,
)
from diffpy.datasets import load_dataset, load_ppi

chu = load_dataset("chu")
network = load_ppi("ppi_PC_2016")

integrated = DoIntegPPI(chu["expression"], network)
entropy = CompSRana(integrated, n_jobs=4)
ccat = CompCCAT(chu["expression"], network)

states = InferPotencyStates(
    entropy.signaling_entropy,
    score_type="signaling_entropy",
    phenotype=chu["phenotype"],
)

print(entropy.signaling_entropy.head())
print(states.distribution)
```

## Diffusion and root inference

```python
from diffpy import (
    CompCCAT,
    compute_diffusion_pseudotime,
    InferDMAPandRoot,
)
from diffpy.datasets import load_dataset, load_ppi

liver = load_dataset("liver")
potency = CompCCAT(liver["expression"], load_ppi("ppi_PC_2016"))

trajectory = InferDMAPandRoot(
    potency,
    liver["expression"],
    k_neighbors=30,
    top_fraction=0.05,
)
dpt = compute_diffusion_pseudotime(
    trajectory,
    paths_to=(1, 2, 3),
)

print(trajectory.root_cell)
print(trajectory.diffusion_coordinates.iloc[:, :2].head())
print(dpt.pseudotime.head())
```

`dpt.paths` contains smoothed diffusion-coordinate paths for the requested
one-based branch numbers, reproducing the role of `paths_to=1:3` in the R
analysis.

Diffusion eigenvectors have arbitrary signs. An axis can be reversed without
changing pairwise geometry, neighborhoods, or biological interpretation.

## Regulatory activity

```python
from diffpy import SciraEstRegAct
from diffpy.datasets import load_dataset

stomach = load_dataset("stomach")
activity = SciraEstRegAct(
    stomach["expression"],
    tissue="stomach",
    norm="z",
    n_jobs=4,
)

average_tf_activity = activity.mean(axis=0)
print(average_tf_activity.head())
```

Activity is the t-statistic of a regulon-predictor slope, not the expression
level of the transcription factor itself. The required `tissue` argument
selects one of the nine packaged regulons; custom regulon input is not
supported. Because transcription factors are rows and cells are columns,
`activity.mean(axis=0)` calculates the average TF activity for each cell.

## Datasets

The four bundled PPI networks are `ppi_PC_2012`, `ppi_PC_2016`,
`ppi_PC_2024`, and `ppi_string_2020`. Pass a full name to `load_ppi`, for example
`load_ppi("ppi_PC_2016")`. Year aliases `"2012"`, `"2016"`, `"2024"`, and
`"2020"` are supported. `load_dataset` accepts all four canonical names.

The examples below assume that ``DIFFPY_DATA_DIR`` points to the extracted
``DiffPy_test_data.zip`` directory, as shown under Installation.

```python
from diffpy.datasets import available_datasets, load_dataset, load_ppi, load_regulon

print(available_datasets())
liver = load_dataset("liver")
network = load_ppi("ppi_PC_2016")
stomach_regulon = load_regulon("stomach")
```

Loaded datasets are ordinary dictionaries:

```python
print(liver["expression"].shape)
print(liver["description"])
```

Annotation vectors correspond to expression columns by position. Their pandas
Series indices are not guaranteed to contain cell identifiers, so preserve
their original order rather than joining annotations to expression columns by
Series index.

See the standalone [HTML tutorial](docs/index.html), [API.md](docs/API.md), and
[DATASETS.md](docs/DATASETS.md) for complete usage and reference documentation.

## License and scientific provenance

DiffPy is licensed under GPL-3.0-only. The implemented methods originate
from the published SCENT and SCIRA research software and associated papers. See
`NOTICE` and `LICENSE`; cite the original method publications when using the
package in scientific work.