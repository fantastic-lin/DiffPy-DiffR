# Convert Python PPI networks to R datasets

Run from the combined project root:

```sh
python tools/convert_ppi_to_r.py
```

Requires NumPy, SciPy, `Rscript`, and the R `Matrix` package. Optionally pass
`ppi_PC_2024` or `ppi_string_2020` to convert only one network.

The script writes `DiffR/data/ppi_PC_2024.rda` and
`DiffR/data/ppi_string_2020.rda`, containing `ppi_PC_2024.m` and
`ppi_string_2020.m`, respectively. These are dense numeric matrices, matching
the representation in `ppi_PC_2012.rda`. Row/column gene identifiers, order,
and all adjacency values are preserved. The largest dense matrix alone uses
about 3.2 GB; allow extra memory for validation and R's copies.

The script validates binary values, symmetry, zero diagonals, and unique
matching labels, saves with xz compression, then reloads and checks all values
and labels before replacing the destination file. Temporary files are created
under `/tmp` and automatically removed.

After installing the updated DiffR package:

```r
data(ppi_PC_2024, package = "DiffR")
dim(ppi_PC_2024.m)
data(ppi_string_2020, package = "DiffR")
dim(ppi_string_2020.m)
```
