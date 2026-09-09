"""Fast signaling-entropy calculations for dense and sparse expression data."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from numba import njit
from scipy import sparse
from scipy.sparse.linalg import eigsh

from .integration import IntegrationResult


@dataclass(frozen=True)
class SignalingEntropyResult:
    """Cellular signaling entropy and its gene-level components.

    Parameters
    ----------
    signaling_entropy
        Normalized global entropy for every cell. Values are generally between
        zero and one, with larger values indicating a more promiscuous
        signaling configuration.
    stationary_distribution
        Gene-by-cell stationary probabilities of the expression-weighted
        random walk. Each valid column sums to one.
    local_entropy
        Unnormalized gene-by-cell entropy of outgoing transition probabilities.
    normalized_local_entropy
        Optional local entropy normalized by the logarithm of each gene's
        positive neighbor count. Present only when requested.
    maximum_entropy
        Logarithm of the leading network eigenvalue used to normalize global
        entropy.
    """

    signaling_entropy: pd.Series
    stationary_distribution: pd.DataFrame
    local_entropy: pd.DataFrame
    normalized_local_entropy: pd.DataFrame | None
    maximum_entropy: float

    def as_dict(self) -> dict[str, Any]:
        """Return all result components in a plain dictionary.

        Returns
        -------
        dict
            Mapping from descriptive Python field names to the corresponding
            pandas objects or normalization constant.
        """

        return {
            "signaling_entropy": self.signaling_entropy,
            "stationary_distribution": self.stationary_distribution,
            "local_entropy": self.local_entropy,
            "normalized_local_entropy": self.normalized_local_entropy,
            "maximum_entropy": self.maximum_entropy,
        }


@njit(cache=True, nogil=True)
def _one_cell(
    expression: np.ndarray,
    indptr: np.ndarray,
    indices: np.ndarray,
    edge_values: np.ndarray,
    maximum_entropy: float,
    include_normalized_local: bool,
) -> tuple[float, np.ndarray, np.ndarray, np.ndarray]:
    """Compute all entropy quantities for one expression profile.

    Parameters
    ----------
    expression
        Dense gene-expression vector for one cell.
    indptr, indices, edge_values
        Compressed sparse row representation of the PPI network.
    maximum_entropy
        Network-level normalization constant.
    include_normalized_local
        Whether to calculate degree-normalized local entropy.

    Returns
    -------
    tuple
        Global entropy, stationary distribution, local entropy, and normalized
        local entropy. The final vector contains ``NaN`` when it was not
        requested.

    Notes
    -----
    Numba compiles this loop once per compatible numeric signature. Explicit
    CSR traversal avoids creating a dense cell-specific transition matrix.
    """

    n_genes = expression.size
    neighbor_sum = np.zeros(n_genes)
    stationary = np.zeros(n_genes)
    local_entropy = np.zeros(n_genes)
    normalized = np.empty(n_genes)
    normalized[:] = np.nan

    for gene in range(n_genes):
        total = 0.0
        for position in range(indptr[gene], indptr[gene + 1]):
            total += edge_values[position] * expression[indices[position]]
        neighbor_sum[gene] = total
        stationary[gene] = expression[gene] * total

    normalization = stationary.sum()
    if not np.isfinite(normalization) or normalization <= 0.0:
        stationary[:] = np.nan
        local_entropy[:] = np.nan
        return np.nan, stationary, local_entropy, normalized
    stationary /= normalization

    for gene in range(n_genes):
        denominator = neighbor_sum[gene]
        if not np.isfinite(denominator) or denominator <= 0.0:
            local_entropy[gene] = 0.0
            if include_normalized_local:
                normalized[gene] = 0.0
            continue
        entropy = 0.0
        positive_neighbors = 0
        for position in range(indptr[gene], indptr[gene + 1]):
            probability = (
                edge_values[position]
                * expression[indices[position]]
                / denominator
            )
            if probability > 0.0 and np.isfinite(probability):
                entropy -= probability * np.log(probability)
                positive_neighbors += 1
        local_entropy[gene] = entropy
        if include_normalized_local:
            normalized[gene] = (
                entropy / np.log(positive_neighbors)
                if positive_neighbors > 1
                else 0.0
            )

    global_entropy = np.dot(stationary, local_entropy)
    if np.isfinite(maximum_entropy) and maximum_entropy != 0.0:
        global_entropy /= maximum_entropy
    return global_entropy, stationary, local_entropy, normalized


def _validate_integration(
    integrated: IntegrationResult,
) -> tuple[Any, sparse.csr_matrix]:
    """Validate an integration result and construct CSR adjacency.

    Parameters
    ----------
    integrated
        Output of :func:`diffpy.integrate_expression_network`.

    Returns
    -------
    tuple
        The labeled expression matrix and sorted CSR adjacency matrix.

    Raises
    ------
    TypeError
        If ``integrated`` is not an :class:`IntegrationResult`.
    ValueError
        If fewer than two genes are present.
    """

    if not isinstance(integrated, IntegrationResult):
        raise TypeError(
            "integrated must be an IntegrationResult returned by "
            "integrate_expression_network"
        )
    adjacency = sparse.csr_matrix(integrated.adjacency.values, dtype=float)
    adjacency.sort_indices()
    if adjacency.shape[0] < 2:
        raise ValueError("the integrated network must contain at least two genes")
    return integrated.expression, adjacency


def compute_signaling_entropy(
    integrated: IntegrationResult,
    *,
    include_normalized_local: bool = False,
    n_jobs: int = 1,
) -> SignalingEntropyResult:
    """Compute signaling entropy for every cell in an integrated network.

    Parameters
    ----------
    integrated
        Aligned expression and adjacency matrices returned by
        :func:`~diffpy.integrate_expression_network`. Dense and sparse
        expression and adjacency matrices are all supported.
    include_normalized_local
        Calculate local entropy normalized by the logarithm of the number of
        positive transition probabilities. This adds one gene-by-cell output
        matrix and therefore increases memory use.
    n_jobs
        Number of parallel worker threads. Use ``-1`` to request all available
        processors. Numba releases the interpreter lock inside the compiled
        per-cell calculation.

    Returns
    -------
    SignalingEntropyResult
        Global cellular entropy, stationary probabilities, local entropy, the
        optional normalized local entropy, and the network normalization
        constant.

    Raises
    ------
    TypeError
        If ``integrated`` is not an :class:`IntegrationResult`.
    ValueError
        If ``n_jobs`` is zero, the network is too small, or its leading
        eigenvalue is not positive.

    Examples
    --------
    >>> from diffpy import integrate_expression_network, compute_signaling_entropy
    >>> from diffpy.datasets import load_dataset, load_ppi
    >>> chu = load_dataset("chu", data_dir="/path/to/diffpy-data")
    >>> integrated = integrate_expression_network(chu["expression"], load_ppi("2012"))
    >>> result = compute_signaling_entropy(integrated, n_jobs=4)
    >>> result.signaling_entropy.shape[0] == chu["expression"].shape[1]
    True

    Notes
    -----
    For each cell, expression weights network transitions. The stationary
    distribution weights local Shannon entropies to produce a global entropy
    rate, which is divided by the logarithm of the leading adjacency
    eigenvalue. Expression columns are converted individually when the input is
    sparse; the complete expression matrix is never densified here.
    """

    if n_jobs == 0:
        raise ValueError("n_jobs cannot be zero")
    expression, adjacency = _validate_integration(integrated)
    leading_eigenvalue = float(
        eigsh(adjacency, k=1, which="LM", return_eigenvectors=False)[0]
    )
    if leading_eigenvalue <= 0 or not np.isfinite(leading_eigenvalue):
        raise ValueError("the leading adjacency eigenvalue must be positive")
    maximum_entropy = float(np.log(leading_eigenvalue))

    expression_csc = (
        expression.values.tocsc() if sparse.issparse(expression.values) else None
    )

    def calculate(cell: int) -> tuple[float, np.ndarray, np.ndarray, np.ndarray]:
        """Extract one cell and dispatch the compiled entropy calculation."""

        if expression_csc is not None:
            vector = np.asarray(expression_csc[:, cell].toarray()).ravel()
        else:
            vector = np.asarray(expression.values[:, cell], dtype=float).ravel()
        return _one_cell(
            vector,
            adjacency.indptr,
            adjacency.indices,
            adjacency.data,
            maximum_entropy,
            include_normalized_local,
        )

    outputs = Parallel(n_jobs=n_jobs, prefer="threads")(
        delayed(calculate)(cell) for cell in range(expression.shape[1])
    )
    global_entropy = np.asarray([output[0] for output in outputs])
    stationary = np.column_stack([output[1] for output in outputs])
    local_entropy = np.column_stack([output[2] for output in outputs])
    normalized = (
        np.column_stack([output[3] for output in outputs])
        if include_normalized_local
        else None
    )
    return SignalingEntropyResult(
        signaling_entropy=pd.Series(
            global_entropy, index=expression.col_names, name="signaling_entropy"
        ),
        stationary_distribution=pd.DataFrame(
            stationary, index=expression.row_names, columns=expression.col_names
        ),
        local_entropy=pd.DataFrame(
            local_entropy, index=expression.row_names, columns=expression.col_names
        ),
        normalized_local_entropy=(
            pd.DataFrame(
                normalized, index=expression.row_names, columns=expression.col_names
            )
            if normalized is not None
            else None
        ),
        maximum_entropy=maximum_entropy,
    )
