# Reconstruction validation

Sources: `DiffPy.zip` and `DiffR_1.0.3.tar.gz`; checksums are recorded in
`SOURCE_ARCHIVES.sha256`. All 108 archive files were verified byte-for-byte
against the extracted package directories. Package implementations were unchanged.
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
