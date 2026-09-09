"""Labeled dense/sparse matrix helpers used throughout the package."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Sequence
import warnings

import numpy as np
import pandas as pd
from scipy import sparse


MatrixLike = np.ndarray | sparse.spmatrix


def _names(values: Iterable[Any] | None, size: int, axis: str) -> tuple[str, ...]:
    """Normalize optional axis labels into a validated string tuple.

    Parameters
    ----------
    values
        Candidate labels, or ``None`` to generate zero-based string labels.
    size
        Required number of labels.
    axis
        Human-readable axis name used in validation errors.

    Returns
    -------
    tuple of str
        Unique labels with exactly ``size`` elements.

    Raises
    ------
    ValueError
        If the number of labels is incorrect or labels are duplicated.
    """

    if values is None:
        return tuple(str(i) for i in range(size))
    result = tuple(str(x) for x in values)
    if len(result) != size:
        raise ValueError(f"{axis}_names has length {len(result)}; expected {size}")
    if len(set(result)) != len(result):
        raise ValueError(f"{axis}_names must be unique")
    return result


@dataclass(frozen=True)
class LabeledMatrix:
    """A numeric dense or sparse matrix with stable axis identifiers.

    Parameters
    ----------
    values
        Two-dimensional NumPy array or SciPy sparse matrix. Dense values are
        stored as floating-point arrays; sparse values are normalized to CSR.
    row_names
        Unique labels corresponding to matrix rows.
    col_names
        Unique labels corresponding to matrix columns.

    Examples
    --------
    >>> import numpy as np
    >>> from diffpy import LabeledMatrix
    >>> matrix = LabeledMatrix(
    ...     np.eye(2), ("gene_a", "gene_b"), ("cell_1", "cell_2")
    ... )
    >>> matrix.shape
    (2, 2)
    """

    values: MatrixLike
    row_names: tuple[str, ...]
    col_names: tuple[str, ...]

    def __post_init__(self) -> None:
        """Normalize numeric storage and validate dimensions and labels.

        Raises
        ------
        ValueError
            If values are not two-dimensional or axis labels have invalid
            length or duplicates.
        """

        values = self.values
        if sparse.issparse(values):
            values = values.tocsr().astype(float)
        else:
            values = np.asarray(values, dtype=float)
        if values.ndim != 2:
            raise ValueError("values must be a two-dimensional matrix")
        rows = _names(self.row_names, values.shape[0], "row")
        cols = _names(self.col_names, values.shape[1], "column")
        object.__setattr__(self, "values", values)
        object.__setattr__(self, "row_names", rows)
        object.__setattr__(self, "col_names", cols)

    @property
    def shape(self) -> tuple[int, int]:
        """Return matrix dimensions as ``(rows, columns)``."""

        return self.values.shape

    @property
    def is_sparse(self) -> bool:
        """Return whether values use a SciPy sparse representation."""

        return sparse.issparse(self.values)

    def toarray(self) -> np.ndarray:
        """Return a dense NumPy representation.

        Returns
        -------
        numpy.ndarray
            Dense floating-point matrix. Sparse inputs are materialized and
            may require substantially more memory.
        """

        return self.values.toarray() if self.is_sparse else np.asarray(self.values)

    def to_dataframe(self) -> pd.DataFrame:
        """Convert values and labels to a pandas DataFrame.

        Returns
        -------
        pandas.DataFrame
            Dense DataFrame for dense storage, or a DataFrame with pandas
            sparse columns for sparse storage.
        """

        if self.is_sparse:
            return pd.DataFrame.sparse.from_spmatrix(
                self.values, index=self.row_names, columns=self.col_names
            )
        return pd.DataFrame(self.values, index=self.row_names, columns=self.col_names)

    def subset_rows(self, names: Sequence[str]) -> "LabeledMatrix":
        """Select and reorder rows by label.

        Parameters
        ----------
        names
            Requested row labels in output order.

        Returns
        -------
        LabeledMatrix
            New matrix containing the selected rows.

        Raises
        ------
        KeyError
            If any requested row label is absent.
        """

        index = _index_for(self.row_names, names, "row")
        return LabeledMatrix(self.values[index, :], tuple(names), self.col_names)

    def subset_columns(self, names: Sequence[str]) -> "LabeledMatrix":
        """Select and reorder columns by label.

        Parameters
        ----------
        names
            Requested column labels in output order.

        Returns
        -------
        LabeledMatrix
            New matrix containing the selected columns.

        Raises
        ------
        KeyError
            If any requested column label is absent.
        """

        index = _index_for(self.col_names, names, "column")
        return LabeledMatrix(self.values[:, index], self.row_names, tuple(names))


def _index_for(existing: Sequence[str], requested: Sequence[str], axis: str) -> np.ndarray:
    """Translate requested labels into integer positions.

    Parameters
    ----------
    existing
        Complete ordered axis labels.
    requested
        Labels to locate, in desired output order.
    axis
        Axis description used in errors.

    Returns
    -------
    numpy.ndarray
        Integer positions corresponding to ``requested``.

    Raises
    ------
    KeyError
        If one or more requested labels are not present.
    """

    lookup = {name: i for i, name in enumerate(existing)}
    missing = [str(name) for name in requested if str(name) not in lookup]
    if missing:
        preview = ", ".join(missing[:5])
        raise KeyError(f"Unknown {axis} name(s): {preview}")
    return np.asarray([lookup[str(name)] for name in requested], dtype=int)


def as_labeled_matrix(
    matrix: Any,
    row_names: Iterable[Any] | None = None,
    col_names: Iterable[Any] | None = None,
    *,
    require_names: bool = False,
) -> LabeledMatrix:
    """Convert common matrix containers into a :class:`LabeledMatrix`.

    Parameters
    ----------
    matrix
        Existing :class:`LabeledMatrix`, pandas DataFrame, NumPy array, or
        SciPy sparse matrix.
    row_names
        Optional row labels for unlabeled arrays.
    col_names
        Optional column labels for unlabeled arrays.
    require_names
        Require explicit labels for unlabeled input instead of generating
        zero-based strings.

    Returns
    -------
    LabeledMatrix
        Validated floating-point matrix with unique string labels.

    Raises
    ------
    ValueError
        If labels conflict with an already labeled input, required labels are
        absent, dimensions are invalid, or labels are duplicated.

    Examples
    --------
    >>> import pandas as pd
    >>> from diffpy import as_labeled_matrix
    >>> frame = pd.DataFrame([[1, 2]], index=["gene"], columns=["a", "b"])
    >>> as_labeled_matrix(frame).row_names
    ('gene',)
    """

    if isinstance(matrix, LabeledMatrix):
        if row_names is not None or col_names is not None:
            raise ValueError("Do not pass names with an existing LabeledMatrix")
        return matrix
    if isinstance(matrix, pd.DataFrame):
        if row_names is not None or col_names is not None:
            raise ValueError("A DataFrame already supplies row and column names")
        if any(isinstance(dtype, pd.SparseDtype) for dtype in matrix.dtypes):
            values = matrix.sparse.to_coo().tocsr().astype(float)
        else:
            values = matrix.to_numpy(dtype=float, copy=False)
        return LabeledMatrix(values, tuple(map(str, matrix.index)), tuple(map(str, matrix.columns)))
    if sparse.issparse(matrix):
        values = matrix.tocsr().astype(float)
    else:
        values = np.asarray(matrix, dtype=float)
    if values.ndim != 2:
        raise ValueError("matrix must be two-dimensional")
    if require_names and (row_names is None or col_names is None):
        raise ValueError("row_names and col_names are required for unlabeled matrix input")
    return LabeledMatrix(
        values,
        _names(row_names, values.shape[0], "row"),
        _names(col_names, values.shape[1], "column"),
    )


def matrix_max(matrix: MatrixLike) -> float:
    """Return the maximum stored or dense matrix value.

    Parameters
    ----------
    matrix
        Dense or sparse numeric matrix.

    Returns
    -------
    float
        Maximum value, treating an empty sparse matrix as all zeros.
    """

    if sparse.issparse(matrix):
        return float(matrix.max()) if matrix.nnz else 0.0
    return float(np.nanmax(matrix))


def log2_plus_one(matrix: MatrixLike) -> MatrixLike:
    """Compute ``log2(x + 1)`` while preserving sparse zeros.

    Parameters
    ----------
    matrix
        Dense or sparse nonnegative numeric matrix.

    Returns
    -------
    numpy.ndarray or scipy.sparse.spmatrix
        Transformed values in the same storage family as the input.

    Warns
    -----
    RuntimeWarning
        Every time the transformation is applied, so automatic preprocessing
        is visible to the caller.
    """

    warnings.warn(
        "Applying automatic log2(x + 1) transformation to expression values.",
        RuntimeWarning,
        stacklevel=2,
    )
    if sparse.issparse(matrix):
        result = matrix.copy().astype(float)
        result.data = np.log2(result.data + 1.0)
        result.eliminate_zeros()
        return result
    return np.log2(np.asarray(matrix, dtype=float) + 1.0)


def row_sums(matrix: MatrixLike) -> np.ndarray:
    """Return row sums as a one-dimensional dense array.

    Parameters
    ----------
    matrix
        Dense or sparse numeric matrix.

    Returns
    -------
    numpy.ndarray
        One sum per matrix row.
    """

    return np.asarray(matrix.sum(axis=1)).ravel()


def column_sums(matrix: MatrixLike) -> np.ndarray:
    """Return column sums as a one-dimensional dense array.

    Parameters
    ----------
    matrix
        Dense or sparse numeric matrix.

    Returns
    -------
    numpy.ndarray
        One sum per matrix column.
    """

    return np.asarray(matrix.sum(axis=0)).ravel()
