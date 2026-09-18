"""Estimate tissue-specific transcription-factor regulatory activity."""

from __future__ import annotations

import re
from typing import Any, Iterable

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from numba import njit
from scipy import sparse

from ._matrix import LabeledMatrix, as_labeled_matrix, log2_plus_one, matrix_max


SUPPORTED_TISSUES = (
    "stomach",
    "skin",
    "esophagus",
    "liver",
    "lung",
    "pancreas",
    "colon",
    "breast",
    "kidney",
)


def _canonical_tissue(tissue: str) -> str:
    """Normalize and validate a tissue name.

    Parameters
    ----------
    tissue
        Tissue label. Matching is case-insensitive and ignores punctuation.

    Returns
    -------
    str
        Canonical lowercase tissue label used by the dataset loader.

    Raises
    ------
    TypeError
        If ``tissue`` is not a string.
    ValueError
        If no packaged regulon is available for the requested tissue.
    """

    if not isinstance(tissue, str):
        raise TypeError("tissue must be one string")
    key = re.sub(r"[^a-z0-9]", "", tissue.lower())
    lookup = {re.sub(r"[^a-z0-9]", "", name): name for name in SUPPORTED_TISSUES}
    if key not in lookup:
        raise ValueError(
            f"No packaged SCIRA network is available for {tissue!r}. "
            f"Supported tissues: {', '.join(SUPPORTED_TISSUES)}."
        )
    return lookup[key]


def _sparse_row_stats(matrix: sparse.spmatrix) -> tuple[np.ndarray, np.ndarray]:
    """Calculate finite-aware mean and sample deviation for sparse rows.

    Parameters
    ----------
    matrix
        Gene-by-cell SciPy sparse matrix. Implicit zeros are included as real
        observations, while explicitly stored nonfinite values are excluded.

    Returns
    -------
    tuple of numpy.ndarray
        Row means and sample standard deviations using ``ddof=1``.
    """

    csr = matrix.tocsr().astype(float)
    n_cols = csr.shape[1]
    counts = np.full(csr.shape[0], n_cols, dtype=int)
    if np.any(~np.isfinite(csr.data)):
        coo = csr.tocoo()
        bad = ~np.isfinite(coo.data)
        np.add.at(counts, coo.row[bad], -1)
        csr = csr.copy()
        csr.data[~np.isfinite(csr.data)] = 0.0
        csr.eliminate_zeros()
    sums = np.asarray(csr.sum(axis=1)).ravel()
    sums_sq = np.asarray(csr.multiply(csr).sum(axis=1)).ravel()
    means = np.divide(sums, counts, out=np.full_like(sums, np.nan), where=counts > 0)
    variance = np.full_like(sums, np.nan)
    valid = counts > 1
    variance[valid] = (
        sums_sq[valid] - np.square(sums[valid]) / counts[valid]
    ) / (counts[valid] - 1)
    variance[(variance < 0) & (np.abs(variance) < np.sqrt(np.finfo(float).eps))] = 0
    return means, np.sqrt(variance)


def _dense_row_stats(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Calculate finite-aware mean and sample deviation for dense rows.

    Parameters
    ----------
    matrix
        Dense gene-by-cell numeric matrix.

    Returns
    -------
    tuple of numpy.ndarray
        Row means and sample standard deviations using ``ddof=1``. Nonfinite
        observations are excluded.
    """

    finite = np.where(np.isfinite(matrix), matrix, np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.nanmean(finite, axis=1), np.nanstd(finite, axis=1, ddof=1)


def _activity(expression: np.ndarray, regulon: np.ndarray) -> np.ndarray:
    """Compute slope t-statistics for one cell against every regulon.

    Parameters
    ----------
    expression
        Normalized expression vector with one value per overlapping gene.
    regulon
        Gene-by-transcription-factor predictor matrix.

    Returns
    -------
    numpy.ndarray
        One regression slope t-statistic per transcription factor. A value is
        ``NaN`` when fewer than three observations or no predictor variation is
        available.
    """

    result = np.full(regulon.shape[1], np.nan)
    for j in range(regulon.shape[1]):
        predictor = regulon[:, j]
        valid = np.isfinite(expression) & np.isfinite(predictor)
        n = int(valid.sum())
        if n < 3:
            continue
        x = predictor[valid]
        y = expression[valid]
        dx = x - x.mean()
        dy = y - y.mean()
        sxx = float(np.dot(dx, dx))
        if sxx <= 0:
            continue
        slope = float(np.dot(dx, dy) / sxx)
        residual = dy - slope * dx
        mse = float(np.dot(residual, residual) / (n - 2))
        standard_error = np.sqrt(mse / sxx)
        if standard_error > 0 and np.isfinite(standard_error):
            result[j] = slope / standard_error
    return result


_activity_numba = njit(cache=True, nogil=True)(_activity)


def SciraEstRegAct(
    expression: Any,
    tissue: str,
    *,
    expression_row_names: Iterable[Any] | None = None,
    expression_col_names: Iterable[Any] | None = None,
    norm: str = "z",
    n_jobs: int = 1,
) -> pd.DataFrame | pd.Series:
    """Estimate transcription-factor activity from regulon target expression.

    Parameters
    ----------
    expression
        A named pandas Series for one cell, or a numeric gene-by-cell matrix as
        :class:`~diffpy.LabeledMatrix`, pandas DataFrame, NumPy array, or
        SciPy sparse matrix. Values should be normalized by library size. If
        the maximum exceeds 100, values are treated as untransformed and
        converted with ``log2(x + 1)``. Gene identifiers should be human
        Entrez IDs when a packaged regulon is used.
    tissue
        Tissue whose packaged regulon is loaded. Supported values are
        ``"stomach"``, ``"skin"``, ``"esophagus"``, ``"liver"``, ``"lung"``,
        ``"pancreas"``, ``"colon"``, ``"breast"``, and ``"kidney"``.
    expression_row_names
        Gene identifiers for unlabeled expression arrays.
    expression_col_names
        Cell identifiers for unlabeled expression arrays.
    norm
        Gene-wise preprocessing applied across cells: ``"z"`` subtracts
        the row mean and divides by sample standard deviation; ``"c"``
        subtracts only the row mean.
    n_jobs
        Number of worker threads used for matrix input.

    Returns
    -------
    pandas.DataFrame or pandas.Series
        For matrix input, a transcription-factor-by-cell activity matrix. For
        Series input, one activity value per transcription factor.

    Raises
    ------
    ValueError
        If the tissue or normalization is unknown, labels are duplicated, or
        fewer than three genes overlap the packaged tissue regulon.

    Warns
    -----
    RuntimeWarning
        If the maximum expression value exceeds 100 and automatic
        ``log2(x + 1)`` transformation is applied.

    Examples
    --------
    >>> from diffpy import SciraEstRegAct
    >>> from diffpy.datasets import load_dataset
    >>> stomach = load_dataset("stomach", data_dir="/path/to/diffpy-data")
    >>> activity = SciraEstRegAct(
    ...     stomach["expression"], tissue="stomach", norm="z"
    ... )
    >>> activity.shape[1] == stomach["expression"].shape[1]
    True
    >>> average_tf_activity = activity.mean(axis=0)
    >>> average_tf_activity.shape[0] == stomach["expression"].shape[1]
    True

    Notes
    -----
    Activity is the t-statistic of the regulon-predictor slope in a linear
    regression with an intercept. It is not the expression of the transcription
    factor itself. Only the nine packaged tissue regulons are supported. Sparse
    expression is normalized one cell at a time, avoiding a fully centered
    dense matrix.
    """

    normalization_key = norm.lower().replace("-", "").replace("_", "")
    normalization_lookup = {
        "zscore": "z",
        "z": "z",
        "center": "c",
        "centre": "c",
        "c": "c",
    }
    if normalization_key not in normalization_lookup:
        raise ValueError("norm must be 'z' or 'c'")
    normalization_mode = normalization_lookup[normalization_key]
    if n_jobs == 0:
        raise ValueError("n_jobs cannot be zero")
    from .datasets import load_regulon

    reg = load_regulon(_canonical_tissue(tissue))

    if isinstance(expression, pd.Series):
        if expression.index.has_duplicates:
            raise ValueError("Vector Entrez Gene IDs must be unique")
        values = expression.copy().astype(float)
        if matrix_max(values.to_numpy()) > 100:
            values.iloc[:] = log2_plus_one(values.to_numpy())
        values.index = values.index.map(str)
        common = tuple(g for g in map(str, expression.index) if g in set(reg.row_names))
        if len(common) < 3:
            raise ValueError("Fewer than three genes overlap expression and regulon")
        y = values.loc[list(common)].to_numpy(dtype=float)
        x = reg.subset_rows(common).toarray()
        return pd.Series(_activity(y, x), index=reg.col_names, name=expression.name)

    exp = as_labeled_matrix(
        expression,
        expression_row_names,
        expression_col_names,
        require_names=True,
    )
    if matrix_max(exp.values) > 100:
        exp = LabeledMatrix(
            log2_plus_one(exp.values), exp.row_names, exp.col_names
        )
    reg_names = set(reg.row_names)
    common = tuple(g for g in exp.row_names if g in reg_names)
    if len(common) < 3:
        raise ValueError("Fewer than three genes overlap expression and regulon")
    exp = exp.subset_rows(common)
    reg_values = reg.subset_rows(common).toarray()
    if sparse.issparse(exp.values):
        means, stds = _sparse_row_stats(exp.values)
        exp_csc = exp.values.tocsc()

        def normalized_column(j: int) -> np.ndarray:
            """Normalize one sparse expression column without dense centering."""

            vector = exp_csc[:, j].toarray().ravel() - means
            if normalization_mode == "z":
                variable = np.isfinite(stds) & (stds > 0)
                constant = np.isfinite(stds) & (stds == 0)
                vector[variable] /= stds[variable]
                vector[constant & np.isfinite(vector)] = 0.0
            return vector

    else:
        dense = np.asarray(exp.values, dtype=float)
        means, stds = _dense_row_stats(dense)

        def normalized_column(j: int) -> np.ndarray:
            """Normalize one column from an in-memory dense expression matrix."""

            vector = dense[:, j].copy() - means
            if normalization_mode == "z":
                variable = np.isfinite(stds) & (stds > 0)
                constant = np.isfinite(stds) & (stds == 0)
                vector[variable] /= stds[variable]
                vector[constant & np.isfinite(vector)] = 0.0
            return vector

    activity_kernel = _activity if n_jobs == 1 else _activity_numba

    def activity_for_column(column: int) -> np.ndarray:
        """Normalize one cell and calculate all transcription-factor scores."""

        return activity_kernel(normalized_column(column), reg_values)

    columns = Parallel(n_jobs=n_jobs, prefer="threads")(
        delayed(activity_for_column)(column) for column in range(exp.shape[1])
    )
    return pd.DataFrame(
        np.column_stack(columns), index=reg.col_names, columns=exp.col_names
    )
