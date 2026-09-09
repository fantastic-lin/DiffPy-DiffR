"""Integrate single-cell expression with a PPI network."""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Any, Iterable

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components

from ._matrix import LabeledMatrix, as_labeled_matrix, matrix_max


@dataclass(frozen=True)
class IntegrationResult:
    """Expression and adjacency matrices aligned to one connected gene set.

    Parameters
    ----------
    expression
        Gene-by-cell expression matrix restricted to genes in the largest
        connected component of the PPI network.
    adjacency
        Square gene-by-gene adjacency matrix for the same genes and order as
        ``expression``.

    Notes
    -----
    The expression matrix is dense because ``log2(x + 1.1)`` changes every
    structural zero to a nonzero value. The adjacency matrix retains sparse or
    dense storage according to the input network.
    """

    expression: LabeledMatrix
    adjacency: LabeledMatrix

    def __post_init__(self) -> None:
        """Validate that expression and adjacency use the same ordered genes.

        Raises
        ------
        ValueError
            If the adjacency matrix is not square or its genes do not exactly
            match the expression rows.
        """

        if self.adjacency.shape[0] != self.adjacency.shape[1]:
            raise ValueError("adjacency must be square")
        if self.adjacency.row_names != self.adjacency.col_names:
            raise ValueError("adjacency row and column names must be identical")
        if self.expression.row_names != self.adjacency.row_names:
            raise ValueError("expression and adjacency must use the same ordered genes")


def integrate_expression_network(
    expression: Any,
    network: Any,
    *,
    expression_row_names: Iterable[Any] | None = None,
    expression_col_names: Iterable[Any] | None = None,
    network_names: Iterable[Any] | None = None,
    min_overlap: int = 5000,
) -> IntegrationResult:
    """Align expression with a network and retain its largest component.

    Parameters
    ----------
    expression
        Numeric gene-by-cell expression data supplied as a
        :class:`~diffpy.LabeledMatrix`, pandas DataFrame, NumPy array, or
        SciPy sparse matrix. Values must be nonnegative and library-size
        normalized but not log-transformed. This function always applies
        ``log2(x + 1.1)`` before matching expression to the network.
    network
        Square gene-by-gene PPI network adjacency matrix. Row and column names
        must contain the same unique gene identifiers.
    expression_row_names
        Gene identifiers for an unlabeled expression array. Ignored when
        ``expression`` already carries labels.
    expression_col_names
        Cell identifiers for an unlabeled expression array. Ignored when
        ``expression`` already carries labels.
    network_names
        Gene identifiers for both axes of an unlabeled network array.
    min_overlap
        Minimum number of genes shared between expression and network. The
        default of 5,000 protects against mismatched identifier systems.

    Returns
    -------
    IntegrationResult
        Aligned expression and adjacency matrices for the largest connected
        component of the overlapping network.

    Raises
    ------
    ValueError
        If expression contains negative values, labels are missing, the
        network is not square, too few genes overlap, or no connected
        component can be identified.

    Warns
    -----
    RuntimeWarning
        If the maximum expression value is below 100, because the matrix may
        already be log-transformed. The threshold is a heuristic and the
        requested ``log2(x + 1.1)`` transformation is still applied.

    Examples
    --------
    >>> from diffpy import integrate_expression_network
    >>> from diffpy.datasets import load_dataset, load_ppi
    >>> chu = load_dataset("chu", data_dir="/path/to/diffpy-data")
    >>> result = integrate_expression_network(chu["expression"], load_ppi("2012"))
    >>> result.expression.shape[1] == chu["expression"].shape[1]
    True
    >>> result.expression.is_sparse
    False

    Notes
    -----
    Matrix rows are genes and columns are cells throughout the package. Input
    expression is expected to be library-size normalized and unlogged. The
    function always applies ``log2(x + 1.1)``; this necessarily converts sparse
    expression to dense storage because ``log2(1.1)`` is nonzero. The output
    gene order follows the input network order so the expression and adjacency
    matrices can be used together without further matching.
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

    exp_values = exp.values
    stored_values = (
        exp_values.data if sparse.issparse(exp_values) else np.asarray(exp_values)
    )
    if np.any(stored_values < 0):
        raise ValueError(
            "expression values must be nonnegative, library-size normalized, "
            "and not log-transformed"
        )
    if matrix_max(exp_values) < 100:
        warnings.warn(
            "Expression maximum is below 100. The matrix may already be "
            "log-transformed. integrate_expression_network() expects "
            "library-size normalized, unlogged expression because it applies "
            "log2(x + 1.1).",
            RuntimeWarning,
            stacklevel=2,
        )
    dense_values = (
        exp_values.toarray()
        if sparse.issparse(exp_values)
        else np.asarray(exp_values)
    )
    transformed = np.log2(np.asarray(dense_values, dtype=float) + 1.1)
    exp = LabeledMatrix(transformed, exp.row_names, exp.col_names)

    expression_genes = set(exp.row_names)
    common = tuple(name for name in net.row_names if name in expression_genes)
    if len(common) < min_overlap:
        raise ValueError(
            "The expression/network gene overlap is only "
            f"{len(common)}; at least {min_overlap} genes are required. "
            "Check that both inputs use the same identifier system."
        )

    common_network = net.subset_rows(common).subset_columns(common)
    adjacency = sparse.csr_matrix(common_network.values)
    component_count, membership = connected_components(
        adjacency, directed=False, return_labels=True
    )
    if component_count < 1:
        raise ValueError("the overlapping network contains no connected component")
    counts = np.bincount(membership)
    largest = int(np.flatnonzero(counts == counts.max())[0])
    selected_index = np.flatnonzero(membership == largest)
    selected_names = tuple(common[i] for i in selected_index)

    adjacency_values = common_network.values[selected_index, :][:, selected_index]
    if not net.is_sparse:
        adjacency_values = np.asarray(adjacency_values)
    return IntegrationResult(
        expression=exp.subset_rows(selected_names),
        adjacency=LabeledMatrix(adjacency_values, selected_names, selected_names),
    )
