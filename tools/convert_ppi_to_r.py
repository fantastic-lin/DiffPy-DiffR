#!/usr/bin/env python3
"""Convert bundled SciPy PPI networks to labeled dense R matrices.

Requires Python numpy/scipy, Rscript, and the R Matrix package.
Run from any directory: python tools/convert_ppi_to_r.py
Intermediate files use the system temporary directory, never the repository.
"""
import argparse
import gzip
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import io, sparse

R_CONVERT = r'''
args <- commandArgs(TRUE)
library(Matrix)
a <- readMM(args[1])
rows <- readLines(args[2]); cols <- readLines(args[3])
stopifnot(nrow(a) == length(rows), ncol(a) == length(cols))
x <- as.matrix(a)
dimnames(x) <- list(rows, cols)
e <- new.env(parent = emptyenv())
assign(args[5], x, envir = e)
save(list = args[5], envir = e, file = args[4], compress = "xz")
rm(x, e); gc()
e <- new.env()
stopifnot(identical(load(args[4], envir=e), args[5]))
y <- e[[args[5]]]
stopifnot(is.matrix(y), identical(rownames(y), rows), identical(colnames(y), cols))
# Compare every matrix value to the original sparse matrix after reloading.
z <- as(y, "sparseMatrix")
dimnames(z) <- NULL
stopifnot(nnzero(z - a) == 0)
cat(args[5], nrow(y), "genes:", "saved and verified\n")
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'networks', nargs='*', default=['ppi_PC_2024', 'ppi_string_2020'],
        help='Select ppi_PC_2012, ppi_PC_2016, ppi_PC_2024, or ppi_string_2020; defaults to the two newer networks',
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    for name in args.networks:
        if name not in ('ppi_PC_2012', 'ppi_PC_2016', 'ppi_PC_2024', 'ppi_string_2020'):
            parser.error(f'Unsupported network: {name}')
        source = root / 'DiffPy/src/diffpy/data' / name
        a = sparse.load_npz(source / 'adjacency.npz').tocsr()
        a.sum_duplicates(); a.eliminate_zeros(); a.sort_indices()
        labels = []
        for axis in ('rows', 'cols'):
            with gzip.open(source / f'adjacency.{axis}.txt.gz', 'rt') as stream:
                labels.append(stream.read().splitlines())
        if (a.shape != (len(labels[0]), len(labels[1])) or labels[0] != labels[1]
                or len(set(labels[0])) != len(labels[0]) or not all(labels[0])
                or (a != a.T).nnz or np.any(a.diagonal())
                or not np.all(a.data == 1)):
            raise ValueError(f'{name}: expected square, symmetric binary network with unique matching labels and zero diagonal')
        # Explicit /tmp keeps intermediates outside the project even if TMPDIR is overridden.
        with tempfile.TemporaryDirectory(prefix='diff-ppi-', dir='/tmp') as temp:
            temp = Path(temp)
            io.mmwrite(temp / 'network.mtx', a, symmetry='general')
            for axis, values in zip(('rows', 'cols'), labels):
                (temp / f'{axis}.txt').write_text('\n'.join(values) + '\n')
            script = temp / 'convert.R'
            script.write_text(R_CONVERT)
            output = temp / f'{name}.rda'
            subprocess.run(['Rscript', str(script), str(temp/'network.mtx'),
                            str(temp/'rows.txt'), str(temp/'cols.txt'),
                            str(output), name+'.m'], check=True)
            # Publish only a fully validated result; no intermediate file enters the project.
            (root / 'DiffR/data' / output.name).write_bytes(output.read_bytes())


if __name__ == '__main__':
    main()
