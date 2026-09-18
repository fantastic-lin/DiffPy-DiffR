"""Fast connectome-transcriptome correlation for cellular potency analysis."""

from __future__ import annotations

from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy import sparse

from ._matrix import as_labeled_matrix, log2_plus_one, matrix_max, row_sums


def _column_correlations(matrix: Any, target: np.ndarray) -> np.ndarray:
    """Correlate every matrix column with a shared target vector.

    Parameters
    ----------
    matrix
        Dense or sparse observation-by-variable matrix.
    target
        One-dimensional vector with one value per matrix row.

    Returns
    -------
    numpy.ndarray
        Pearson correlation for every matrix column. Columns with zero
        variance, or a constant target, receive ``NaN``.

    Notes
    -----
    The calculation uses centered cross-products, avoiding materialization of
    a centered sparse matrix.
    """

    n_observations = matrix.shape[0]
    y = np.asarray(target, dtype=float)
    sum_y = y.sum()
    sum_y_squared = np.dot(y, y)
    if sparse.issparse(matrix):
        x = matrix.tocsc()
        sum_x = np.asarray(x.sum(axis=0)).ravel()
        sum_x_squared = np.asarray(x.multiply(x).sum(axis=0)).ravel()
        cross_product = np.asarray(x.T @ y).ravel()
    else:
        x = np.asarray(matrix, dtype=float)
        sum_x = x.sum(axis=0)
        sum_x_squared = np.square(x).sum(axis=0)
        cross_product = x.T @ y
    numerator = cross_product - sum_x * sum_y / n_observations
    centered_x = sum_x_squared - np.square(sum_x) / n_observations
    centered_y = sum_y_squared - sum_y * sum_y / n_observations
    denominator = np.sqrt(np.maximum(centered_x, 0.0) * max(centered_y, 0.0))
    result = np.full(matrix.shape[1], np.nan, dtype=float)
    valid = denominator > 0
    result[valid] = numerator[valid] / denominator[valid]
    return np.clip(result, -1.0, 1.0)


def CompCCAT(
    expression: Any,
    network: Any,
    *,
    expression_row_names: Iterable[Any] | None = None,
    expression_col_names: Iterable[Any] | None = None,
    network_names: Iterable[Any] | None = None,
    min_overlap: int = 5000,
) -> pd.Series:
    """Estimate cellular potency from expression and network connectivity.

    Parameters
    ----------
    expression
        Numeric gene-by-cell expression matrix. Accepted inputs are
        :class:`~diffpy.LabeledMatrix`, pandas DataFrame, NumPy array, and
        SciPy sparse matrix. Values with a maximum above 100 are transformed
        with ``log2(x + 1)``; sparse zeros remain zero.
    network
        Square gene-by-gene adjacency matrix for a PPI network. Network and
        expression genes must use the same identifiers.
    expression_row_names
        Gene names for an unlabeled expression array.
    expression_col_names
        Cell names for an unlabeled expression array.
    network_names
        Gene names for both axes of an unlabeled network array.
    min_overlap
        Minimum number of shared genes. Values below 5,000 can be useful for
        controlled tests but are not recommended for biological analysis.

    Returns
    -------
    pandas.Series
        One correlation-based potency score per cell, indexed by cell name.
        Values lie between -1 and 1; larger values indicate stronger
        expression of highly connected genes.

    Raises
    ------
    ValueError
        If labels are missing, the network is not square, or the overlap is
        smaller than ``min_overlap``.

    Warns
    -----
    RuntimeWarning
        If the maximum expression value exceeds 100 and automatic
        ``log2(x + 1)`` transformation is applied.

    Examples
    --------
    >>> from diffpy import CompCCAT
    >>> from diffpy.datasets import load_dataset, load_ppi
    >>> liver = load_dataset("liver", data_dir="/path/to/diffpy-data")
    >>> scores = CompCCAT(liver["expression"], load_ppi("2012"))
    >>> scores.shape[0] == liver["expression"].shape[1]
    True

    Notes
    -----
    Each cell is correlated with gene degree computed across the full network.
    The sparse implementation is algebraically equivalent to dense Pearson
    correlation and does not center the complete matrix in memory.
    """

    exp = as_labeled_matrix(
        expression,
        expression_row_names,
        expression_col_names,
        require_names=True,
    )
    net = as_labeled_matrix(
        network, network_names, network_names, require_names=True
    )
    if net.shape[0] != net.shape[1] or net.row_names != net.col_names:
        raise ValueError(
            "network must be square with identical row and column gene names"
        )
    if min_overlap < 1:
        raise ValueError("min_overlap must be a positive integer")
    values = log2_plus_one(exp.values) if matrix_max(exp.values) > 100 else exp.values

    expression_genes = set(exp.row_names)
    common = tuple(name for name in net.row_names if name in expression_genes)
    if len(common) < min_overlap:
        raise ValueError(
            "The expression/network gene overlap is only "
            f"{len(common)}; at least {min_overlap} genes are required."
        )
    expression_lookup = {name: i for i, name in enumerate(exp.row_names)}
    network_lookup = {name: i for i, name in enumerate(net.row_names)}
    expression_index = np.asarray([expression_lookup[gene] for gene in common])
    network_index = np.asarray([network_lookup[gene] for gene in common])
    selected_expression = values[expression_index, :]
    degree = row_sums(net.values[network_index, :])
    result = _column_correlations(selected_expression, degree)
    return pd.Series(result, index=exp.col_names, name="ccat")
