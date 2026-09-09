# DiffPy and DiffR

Python and R tools for single-cell differentiation potency (SCENT and CCAT),
potency-guided trajectories, and tissue-specific transcription-factor activity
(SCIRA). Choose the package that matches your analysis language.

| Package | Version | Requirements | Guide |
| --- | --- | --- | --- |
| DiffPy | 2.0.0 | Python >= 3.10 | [Python README](DiffPy/README.md) |
| DiffR | 1.0.3 | R >= 3.6; dependencies may require newer R | [R README](DiffR/README.md) |

[Online tutorials](https://fantastic-lin.github.io/DiffPy-DiffR/) include worked
analyses and figures. [Release downloads](https://github.com/fantastic-lin/DiffPy-DiffR/releases)
provide source and example data. The first documentation deployment and release
must be published before these destinations become available.

## Local installation

From this combined project directory:

```sh
python -m pip install ./DiffPy
R CMD INSTALL DiffR
```

Install the dependencies listed in each package before installing offline.
To create an R source archive while retaining prebuilt tutorials, run
`R CMD build --no-build-vignettes DiffR`. Full vignette regeneration requires
external example data and the suggested R packages.

## Python installation and first example

Run in a terminal with Python and Git installed, preferably in a virtual environment:

```sh
python -m pip install "git+https://github.com/fantastic-lin/DiffPy-DiffR.git#subdirectory=DiffPy"
```

Download `DiffPy_test_data.zip` from Releases and extract it into `diffpy-data`.
That directory should contain `chu/`, `liver/`, and `stomach/` directly.

```python
from diffpy import compute_ccat
from diffpy.datasets import load_dataset, load_ppi

chu = load_dataset("chu", data_dir="/absolute/path/to/diffpy-data")
scores = compute_ccat(chu["expression"], load_ppi("2012"))
print(scores.head())
```

This prints the first cell potency scores. For the complete
[example script](DiffPy/examples/basic_workflow.py), set `DIFFPY_DATA_DIR` to
the same extracted directory before running the script. Networks and regulons
are bundled; expression example datasets are downloaded separately.

## R installation and first example

Run in R. `marray` is installed through Bioconductor; `destiny` is needed only
for the trajectory example.

```r
install.packages(c("remotes", "BiocManager"))
BiocManager::install("marray", ask = FALSE, update = FALSE)
remotes::install_github("fantastic-lin/DiffPy-DiffR", subdir = "DiffR",
                        dependencies = NA, build_vignettes = FALSE,
                        upgrade = "never")
library(DiffR)
test_data_dir <- Sys.getenv("DIFFR_TEST_DATA_DIR")
if (!nzchar(test_data_dir)) stop("Set DIFFR_TEST_DATA_DIR to the extracted DiffR_test_data directory")
load(file.path(test_data_dir, "dataChu.rda"))
data(net13Jun12, package = "DiffR")
scores <- CompCCAT(exp.m = log2(scChuSparse.m + 1), ppiA.m = net13Jun12.m)
head(scores)
```

For trajectories, install `destiny` using
`BiocManager::install("destiny", ask = FALSE, update = FALSE)`.
Download and extract `DiffR_test_data.tar.gz`, then set `DIFFR_TEST_DATA_DIR`
to the extracted `DiffR_test_data` directory before running the example.
PPI and regulon networks are bundled; expression example datasets are separate.

## Input conventions

Expression matrices have genes in rows and cells in columns. Match gene IDs
to the selected network; the R tutorial uses unique human Entrez Gene IDs.
Read each package's preprocessing instructions carefully: Python
`integrate_expression_network` applies its own log transform, whereas the R
tutorial transforms expression before calling `DoIntegPPI`.

## Documentation and attribution

- [Python API](DiffPy/docs/API.md) and [dataset guide](DiffPy/docs/DATASETS.md)
- [R tutorial source](DiffR/vignettes/DiffR.Rmd), covering SCENT, CCAT and SCIRA
- [Issues](https://github.com/fantastic-lin/DiffPy-DiffR/issues) for support

DiffPy is GPL-3.0-only; see its [LICENSE](DiffPy/LICENSE) and [NOTICE](DiffPy/NOTICE).
DiffR declares GPL-3 in its [DESCRIPTION](DiffR/DESCRIPTION). Original author
and maintainer metadata are preserved. Cite the original methods:

- Teschendorff AE, Enver T. Nature Communications (2017), 8:15599.
  DOI: 10.1038/ncomms15599.
- Teschendorff AE, Wang N. npj Genomic Medicine (2020), 5:43.
  DOI: 10.1038/s41525-020-00151-y.

Package sources were extracted from `DiffPy.zip` and `DiffR_1.0.3.tar.gz`.
See [validation details](VALIDATION.md) for checks performed during reconstruction.
