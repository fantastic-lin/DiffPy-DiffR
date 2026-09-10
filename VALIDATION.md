# Reconstruction validation

Sources: `DiffPy.zip` and `DiffR_1.0.3.tar.gz`; checksums are recorded in
`SOURCE_ARCHIVES.sha256`. At initial reconstruction, all 108 archive files were verified byte-for-byte
against the extracted package directories. Package implementations were unchanged at reconstruction; subsequent dataset updates are listed below.
The shared README, release notes, documentation landing page, and Pages workflow
were adapted from `DiffPy-DiffR copy`; its Git history was not copied.

## Results

- Python: 5 tests and 4 subtests passed using Python 3.12 with the isolated
  validation environment and the reconstructed `DiffPy/src` on PYTHONPATH.
- R: `R CMD INSTALL` from the extracted DiffR directory succeeded on R 4.6.0.
- `CompSRana`: known-answer test passed; Chu subset with 8,393 genes and three
  cells passed, including local entropy bounds, stationary probability sums,
  and agreement between one-core and two-core execution.
- Chu SR values: 0.9156052722, 0.9130773028, 0.9108913302.
- Local README links and the staged documentation landing-page links resolve.

Validation scripts, logs, the isolated R library, and a site preview are in
`../DiffPy-DiffR-validation/`. The Python environment there uses existing system
scientific dependencies plus locally installed igraph and pytest.

## Limits

Existing rendered tutorials were retained; full tutorials and all expression
example datasets were not rerun. Python was tested from source, not from an
installed wheel. R session information emitted a harmless timedatectl warning
because the environment does not expose a system bus. No GitHub push or Pages
deployment was performed.

## Upload preparation (2026-09-10)

The upload includes subsequent local updates: the pancreas regulon now has
20 transcription factors and 971 nonzero entries, updated R tutorial/data
help, two additional R PPI datasets (`ppi_PC_2024` and `ppi_string_2020`),
and the Python-to-R conversion tool. These files are no longer byte-for-byte
copies of the original archives.

Python syntax and the nine bundled regulon matrices were checked successfully
against manifest dimensions, nonzero counts, and label counts. `git diff --check`
passed. The largest project file is approximately 11.6 MiB; no obvious credentials
were found by a basic filename and token-pattern scan.

Runtime tests could not complete on this Mac: Python 3.13 lacks NumPy/SciPy,
and Python 3.10 lacks numba. R is missing igraph, mclust, and qlcMatrix.
The earlier reconstruction results above are historical, not a rerun of this
updated directory. The referenced reconstruction validation directory is not
present alongside this checkout.

For the published examples, separately publish `DiffPy_test_data.zip` and
`DiffR_test_data.tar.gz` as GitHub release assets. Regenerate source release
archives from the updated checkout before publishing them. The tutorial site
also requires GitHub Pages to use GitHub Actions and a successful Pages workflow.
