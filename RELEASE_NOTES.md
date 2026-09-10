# Reconstructed combined project

Contains DiffPy 2.0.0 from `DiffPy.zip` and DiffR 1.0.3 from
`DiffR_1.0.3.tar.gz`. The initial reconstruction preserved package contents from these archives;
this checkout also includes subsequent pancreas regulon updates, additional
R PPI datasets, and the conversion tool (see VALIDATION.md).
Shared documentation and the Pages workflow were adapted from the existing
combined-project copy, with updated DiffR tutorial filenames and external-data instructions.

## Release assets

- `DiffPy.zip`: Python source package.
- `DiffR_1.0.3.tar.gz`: installable R source package.
- `DiffPy_test_data.zip`: separate Python expression example datasets.
- `DiffR_test_data.tar.gz`: separate R expression example datasets.

Networks and regulons are bundled in the packages. Expression example datasets
must be downloaded separately. Publish example-data assets separately and rebuild source archives from the
current checkout before creating a release.
