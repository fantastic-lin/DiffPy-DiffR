"""Infer diffusion geometry and a potency-guided trajectory root cell."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial.distance import pdist, squareform
from scipy.sparse.linalg import eigsh
from sklearn.neighbors import NearestNeighbors

from ._matrix import LabeledMatrix, as_labeled_matrix


@dataclass(frozen=True)
class TrajectoryResult:
    """Diffusion-map outputs and inferred root information.

    Parameters
    ----------
    metadata
        Eigenvalues, density normalization, selected genes, and candidate
        root-state indices used during inference.
    diffusion_coordinates
        Cell-by-component diffusion embedding.
    transition_matrix
        Symmetric cell-by-cell diffusion operator stored as a labeled sparse
        matrix.
    root_index
        Zero-based integer position of the inferred root cell.
    root_cell
        Identifier of the inferred root cell.
    potency_colors
        Plotting colors ordered from low to high potency.
    """

    metadata: dict[str, Any]
    diffusion_coordinates: pd.DataFrame
    transition_matrix: LabeledMatrix
    root_index: int
    root_cell: str
    potency_colors: pd.Series

    def as_dict(self) -> dict[str, Any]:
        """Return all trajectory outputs in a plain dictionary.

        Returns
        -------
        dict
            Metadata, diffusion coordinates, transition matrix, root index,
            root-cell identifier, and potency colors.
        """

        return {
            "metadata": self.metadata,
            "diffusion_coordinates": self.diffusion_coordinates,
            "transition_matrix": self.transition_matrix,
            "root_index": self.root_index,
            "root_cell": self.root_cell,
            "potency_colors": self.potency_colors,
        }


@dataclass(frozen=True)
class DiffusionPseudotimeResult:
    """Diffusion pseudotime, branches, tips, and smoothed lineage paths.

    Parameters
    ----------
    pseudotime
        DPT distance from the selected root to every cell.
    distance_matrix
        Symmetric cell-by-cell DPT distance matrix.
    branches
        One-based hierarchical branch assignments. Missing entries indicate
        that a cell was not subdivided at that hierarchy level.
    tips
        Boolean table marking the three tips used at each branch level.
    tip_indices
        Zero-based indices of the three top-level tips, beginning with the
        selected root.
    tip_cells
        Cell identifiers corresponding to ``tip_indices``.
    paths
        Smoothed paths in diffusion-coordinate space, keyed by one-based
        top-level branch number.
    root_index
        Zero-based root-cell position.
    root_cell
        Root-cell identifier.
    """

    pseudotime: pd.Series
    distance_matrix: pd.DataFrame
    branches: pd.DataFrame
    tips: pd.DataFrame
    tip_indices: tuple[int, int, int]
    tip_cells: tuple[str, str, str]
    paths: dict[int, pd.DataFrame]
    root_index: int
    root_cell: str

    def as_dict(self) -> dict[str, Any]:
        """Return every DPT output in a plain dictionary."""

        return {
            "pseudotime": self.pseudotime,
            "distance_matrix": self.distance_matrix,
            "branches": self.branches,
            "tips": self.tips,
            "tip_indices": self.tip_indices,
            "tip_cells": self.tip_cells,
            "paths": self.paths,
            "root_index": self.root_index,
            "root_cell": self.root_cell,
        }


def _row_mean_sd(matrix: Any) -> tuple[np.ndarray, np.ndarray]:
    """Calculate per-gene means and sample deviations.

    Parameters
    ----------
    matrix
        Dense or sparse gene-by-cell matrix.

    Returns
    -------
    tuple of numpy.ndarray
        Mean and sample standard deviation for every gene.

    Raises
    ------
    ValueError
        If fewer than two cells are present.
    """

    n = matrix.shape[1]
    if n < 2:
        raise ValueError("expression must contain at least two cells")
    if sparse.issparse(matrix):
        mean = np.asarray(matrix.mean(axis=1)).ravel()
        sum_sq = np.asarray(matrix.multiply(matrix).sum(axis=1)).ravel()
        variance = np.maximum(0.0, (sum_sq - n * mean * mean) / (n - 1))
        return mean, np.sqrt(variance)
    values = np.asarray(matrix, dtype=float)
    return np.mean(values, axis=1), np.std(values, axis=1, ddof=1)


def _diffusion_operator(
    cell_by_gene: Any, k: int
) -> tuple[sparse.csr_matrix, np.ndarray, np.ndarray]:
    """Construct a density-normalized adaptive Gaussian diffusion operator.

    Parameters
    ----------
    cell_by_gene
        Dense or sparse cell-by-gene matrix containing selected variable genes.
    k
        Number of nearest neighbors used for local bandwidth estimation.

    Returns
    -------
    tuple
        Symmetric sparse diffusion operator, raw affinity density, and
        density-corrected degree.

    Raises
    ------
    ValueError
        If ``k`` is outside the interval from one to ``n_cells - 1``.

    Notes
    -----
    Each cell's bandwidth is its farthest retained neighbor. The kernel is
    symmetrized, anisotropically density-corrected, and normalized into a
    symmetric operator suitable for stable eigendecomposition.
    """

    n_cells = cell_by_gene.shape[0]
    if not 1 <= k < n_cells:
        raise ValueError(f"k must be between 1 and n_cells - 1 ({n_cells - 1})")
    neighbors = NearestNeighbors(n_neighbors=k + 1).fit(cell_by_gene)
    distances, indices = neighbors.kneighbors(cell_by_gene)
    distances = distances[:, 1:]
    indices = indices[:, 1:]
    sigma = np.maximum(distances[:, -1], np.finfo(float).eps)
    rows = np.repeat(np.arange(n_cells), k)
    cols = indices.ravel()
    d2 = np.square(distances.ravel())
    si = np.repeat(sigma, k)
    sj = sigma[cols]
    scale2 = np.square(si) + np.square(sj)
    weights = np.sqrt(2.0 * si * sj / scale2) * np.exp(-d2 / scale2)
    affinity = sparse.coo_matrix((weights, (rows, cols)), shape=(n_cells, n_cells)).tocsr()
    affinity = affinity.maximum(affinity.T)
    affinity.setdiag(1.0)
    affinity.eliminate_zeros()

    density = np.asarray(affinity.sum(axis=1)).ravel()
    inv_density = sparse.diags(1.0 / np.maximum(density, np.finfo(float).eps))
    kernel = inv_density @ affinity @ inv_density
    degree = np.asarray(kernel.sum(axis=1)).ravel()
    inv_sqrt_degree = sparse.diags(1.0 / np.sqrt(np.maximum(degree, np.finfo(float).eps)))
    symmetric = (inv_sqrt_degree @ kernel @ inv_sqrt_degree).tocsr()
    return symmetric, density, degree


def _root_colors(values: np.ndarray) -> np.ndarray:
    """Map potency quantiles to a light-blue-to-black plotting palette.

    Parameters
    ----------
    values
        One finite potency value per cell.

    Returns
    -------
    numpy.ndarray
        String color code for every input value.
    """

    quantiles = np.quantile(values, np.arange(0.05, 1.0, 0.1))
    low = np.array([135, 206, 235], dtype=float)
    high = np.array([0, 0, 255], dtype=float)
    palette = [
        "#%02x%02x%02x" % tuple(np.rint(low + t * (high - low)).astype(int))
        for t in np.linspace(0.0, 1.0, 10)
    ]
    colors = np.full(values.size, "black", dtype=object)
    for i in range(len(quantiles) - 1, -1, -1):
        colors[values <= quantiles[i]] = palette[i]
    return colors


def _gaussian_smooth(values: np.ndarray, window: float) -> np.ndarray:
    """Reproduce destiny's normalized Gaussian moving-window smoother."""

    length = len(values)
    window_length = int(max(abs(window * length if 0 < window < 1 else window), 1))
    if window_length >= length:
        raise ValueError("path smoothing window must be smaller than the path")
    half = abs(window_length / 2.0)
    offsets = np.arange(window_length) - int(half)
    weights = np.exp(-0.5 * np.square(2.5 * offsets / half))
    weights /= weights.sum()
    left = window_length // 2
    right = window_length - left
    smoothed = np.empty(length, dtype=float)
    for position in range(length):
        start = max(0, position - left)
        stop = min(length, position + right)
        weight_start = start - (position - left)
        selected_weights = weights[weight_start : weight_start + stop - start]
        selected_weights = selected_weights / selected_weights.sum()
        smoothed[position] = np.dot(values[start:stop], selected_weights)
    return smoothed


def _tip_statistics(
    distances: np.ndarray, cells: np.ndarray, supplied_tips: tuple[int, ...]
) -> tuple[tuple[int, int, int], np.ndarray, np.ndarray, float]:
    """Choose three DPT tips and calculate destiny's branching statistic."""

    first = supplied_tips[0]
    first_distances = distances[first, cells]
    second = (
        supplied_tips[1]
        if len(supplied_tips) >= 2
        else int(cells[np.argmax(first_distances)])
    )
    second_distances = distances[second, cells]
    third = (
        supplied_tips[2]
        if len(supplied_tips) == 3
        else int(cells[np.argmax(first_distances + second_distances)])
    )
    combined = first_distances + second_distances
    minimum = float(np.min(combined))
    statistic = float(np.max(combined) / minimum) if minimum > 0 else np.inf
    return (first, second, third), first_distances, second_distances, statistic


def _branch_cut(
    tip_distances: np.ndarray, orders: np.ndarray, branch: int, window: float
) -> np.ndarray:
    """Port destiny's Kendall-based cut for one of three candidate branches."""

    n_cells = len(orders)
    other = [index for index in range(3) if index != branch]
    branch_order = orders[:, branch]
    first = tip_distances[branch_order, other[0]]
    second = tip_distances[branch_order, other[1]]
    scores = np.empty(n_cells - 1, dtype=float)
    for split in range(1, n_cells):
        left_first = first[:split]
        left_second = second[:split]
        right_first = first[split:]
        right_second = second[split:]
        left = np.where(left_first >= first[split], 1.0, -1.0)
        left_other = np.where(left_second >= second[split], 1.0, -1.0)
        right = np.where(right_first >= first[split - 1], 1.0, -1.0)
        right_other = np.where(right_second >= second[split - 1], 1.0, -1.0)
        scores[split - 1] = (
            np.dot(left, left_other) / split
            - np.dot(right, right_other) / (n_cells - split)
        )
    cut = int(np.argmax(_gaussian_smooth(scores, window))) + 1
    return branch_order[:cut]


def _organize_branches(branches: list[np.ndarray]) -> list[np.ndarray]:
    """Remove shared cells while retaining each branch's first ordered cell."""

    intersections: set[int] = set()
    for first in range(3):
        for second in range(first + 1, 3):
            intersections.update(
                np.intersect1d(branches[first], branches[second]).tolist()
            )
    organized: list[np.ndarray] = []
    for branch in branches:
        retained = [int(cell) for cell in branch if int(cell) not in intersections]
        if int(branch[0]) not in retained:
            retained.append(int(branch[0]))
        organized.append(np.asarray(retained, dtype=int))
    return organized


def _automatic_branches(
    distances: np.ndarray,
    cells: np.ndarray,
    supplied_tips: tuple[int, ...],
    window: float,
    *,
    minimum_cells: int = 10,
    minimum_statistic: float = 1.1,
) -> tuple[np.ndarray, np.ndarray, tuple[int, int, int]]:
    """Recursively port destiny's automatic DPT branch assignment."""

    tip_indices, _, _, statistic = _tip_statistics(
        distances, cells, supplied_tips
    )
    if len(cells) < minimum_cells or statistic < minimum_statistic:
        raise ValueError(
            "DPT branch inference requires at least 10 cells and a branching "
            "statistic of at least 1.1"
        )
    tip_distances = distances[np.ix_(cells, np.asarray(tip_indices))]
    orders = np.argsort(tip_distances, axis=0, kind="stable")
    local_branches = _organize_branches(
        [
            _branch_cut(tip_distances, orders, branch, window)
            for branch in range(3)
        ]
    )
    cell_to_local = {int(cell): position for position, cell in enumerate(cells)}
    branch_values = np.full((len(cells), 1), np.nan)
    tips = np.zeros((len(cells), 1), dtype=bool)
    for branch_number, local_positions in enumerate(local_branches, start=1):
        branch_values[local_positions, 0] = branch_number
    for tip in tip_indices:
        tips[cell_to_local[tip], 0] = True

    subdivisions: list[tuple[np.ndarray, np.ndarray, np.ndarray] | None] = []
    for local_positions, tip in zip(local_branches, tip_indices, strict=True):
        subset_cells = cells[local_positions]
        if len(subset_cells) < minimum_cells or tip not in subset_cells:
            subdivisions.append(None)
            continue
        sub_tip_indices, _, _, sub_statistic = _tip_statistics(
            distances, subset_cells, (tip,)
        )
        if sub_statistic < minimum_statistic:
            subdivisions.append(None)
            continue
        try:
            sub_branches, sub_tips, _ = _automatic_branches(
                distances,
                subset_cells,
                sub_tip_indices,
                window,
                minimum_cells=minimum_cells,
                minimum_statistic=minimum_statistic,
            )
        except ValueError:
            subdivisions.append(None)
        else:
            subdivisions.append((local_positions, sub_branches, sub_tips))

    populated = [subdivision for subdivision in subdivisions if subdivision is not None]
    if populated:
        extra_levels = max(subdivision[1].shape[1] for subdivision in populated)
        branch_values = np.column_stack(
            [branch_values, np.full((len(cells), extra_levels), np.nan)]
        )
        tips = np.column_stack(
            [tips, np.zeros((len(cells), extra_levels), dtype=bool)]
        )
        branch_offset = int(np.nanmax(branch_values[:, 0]))
        for local_positions, sub_branches, sub_tips in populated:
            width = sub_branches.shape[1]
            adjusted = sub_branches.copy()
            finite = np.isfinite(adjusted)
            adjusted[finite] += branch_offset
            branch_values[np.ix_(local_positions, np.arange(1, width + 1))] = adjusted
            tips[np.ix_(local_positions, np.arange(1, width + 1))] = sub_tips
            branch_offset = int(np.nanmax(branch_values))
    return branch_values, tips, tip_indices


def compute_diffusion_pseudotime(
    trajectory: TrajectoryResult,
    *,
    root: int | str | None = None,
    paths_to: Iterable[int] = (1, 2, 3),
    window_width: float = 0.1,
) -> DiffusionPseudotimeResult:
    """Compute destiny-style diffusion pseudotime and lineage paths.

    Parameters
    ----------
    trajectory
        Result returned by :func:`infer_diffusion_root`.
    root
        Root cell as a zero-based integer index or cell identifier. The
        trajectory's inferred root is used by default.
    paths_to
        One-based top-level branch numbers for which smoothed paths are
        returned. The default reproduces ``paths_to=1:3`` in the R vignette.
    window_width
        Gaussian smoothing-window width. Values between zero and one are
        interpreted as a fraction of the number of cells, matching destiny.

    Returns
    -------
    DiffusionPseudotimeResult
        Root-based pseudotime, all pairwise DPT distances, hierarchical branch
        and tip tables, and smoothed diffusion-coordinate paths.

    Raises
    ------
    TypeError
        If ``trajectory`` is not a :class:`TrajectoryResult`.
    ValueError
        If root, path, eigenvalue, or smoothing arguments are invalid.

    Examples
    --------
    >>> from diffpy import (
    ...     compute_ccat,
    ...     compute_diffusion_pseudotime,
    ...     infer_diffusion_root,
    ... )
    >>> from diffpy.datasets import load_dataset, load_ppi
    >>> liver = load_dataset("liver", data_dir="/path/to/diffpy-data")
    >>> potency = compute_ccat(liver["expression"], load_ppi("2012"))
    >>> trajectory = infer_diffusion_root(potency, liver["expression"])
    >>> result = compute_diffusion_pseudotime(trajectory)
    >>> result.pseudotime.loc[result.root_cell]
    0.0

    Notes
    -----
    DPT distances use the destiny weighting
    ``lambda / (1 - lambda)`` for every nontrivial diffusion eigenvector.
    Automatic tips, branch cuts, and Gaussian-smoothed paths follow destiny's
    DPT and plotting algorithms. Small numerical differences from R can remain
    because the underlying diffusion operators and eigensolvers differ.
    """

    if not isinstance(trajectory, TrajectoryResult):
        raise TypeError("trajectory must be returned by infer_diffusion_root")
    cells = trajectory.diffusion_coordinates.index.astype(str).tolist()
    if root is None:
        root_index = trajectory.root_index
    elif isinstance(root, str):
        if root not in cells:
            raise ValueError(f"unknown root cell {root!r}")
        root_index = cells.index(root)
    else:
        root_index = int(root)
        if root_index < 0 or root_index >= len(cells):
            raise ValueError("root index is outside the cell range")
    if window_width <= 0:
        raise ValueError("window_width must be positive")

    coordinates = trajectory.diffusion_coordinates.to_numpy(dtype=float)
    eigenvalues = np.asarray(trajectory.metadata["eigenvalues"], dtype=float)[
        1 : coordinates.shape[1] + 1
    ]
    usable = (
        np.isfinite(eigenvalues)
        & (np.abs(eigenvalues) > np.finfo(float).eps)
        & (np.abs(1.0 - eigenvalues) > np.sqrt(np.finfo(float).eps))
    )
    if not np.any(usable):
        raise ValueError("trajectory has no usable nontrivial diffusion eigenvalues")
    eigenvectors = coordinates[:, usable] / eigenvalues[usable]
    weighted = eigenvectors * (
        eigenvalues[usable] / (1.0 - eigenvalues[usable])
    )
    distances = squareform(pdist(weighted, metric="euclidean"))
    all_cells = np.arange(len(cells), dtype=int)
    branches, tips, tip_indices = _automatic_branches(
        distances, all_cells, (root_index,), window_width
    )

    root_branch = int(branches[root_index, 0])
    branch_tip_positions = np.flatnonzero(
        (branches[:, 0] == root_branch) & tips[:, 0]
    )
    reference = (
        int(branch_tip_positions[0])
        if branch_tip_positions.size
        else int(np.flatnonzero(branches[:, 0] == root_branch)[0])
    )
    pseudotime_values = distances[reference, :]
    requested_paths = tuple(int(branch) for branch in paths_to)
    available = set(int(value) for value in branches[:, 0] if np.isfinite(value))
    if any(branch not in available for branch in requested_paths):
        raise ValueError(
            f"paths_to must use available one-based branches {sorted(available)}"
        )
    paths: dict[int, pd.DataFrame] = {}
    for branch in requested_paths:
        selected = np.isin(branches[:, 0], (root_branch, branch))
        order = np.argsort(pseudotime_values[selected], kind="stable")
        path_coordinates = coordinates[selected, :][order, :]
        smoothed = np.column_stack(
            [
                _gaussian_smooth(path_coordinates[:, component], window_width)
                for component in range(path_coordinates.shape[1])
            ]
        )
        paths[branch] = pd.DataFrame(
            smoothed, columns=trajectory.diffusion_coordinates.columns
        )

    branch_frame = pd.DataFrame(
        branches,
        index=cells,
        columns=[f"Branch{level}" for level in range(1, branches.shape[1] + 1)],
    ).astype("Int64")
    tip_frame = pd.DataFrame(
        tips,
        index=cells,
        columns=[f"Tips{level}" for level in range(1, tips.shape[1] + 1)],
    )
    return DiffusionPseudotimeResult(
        pseudotime=pd.Series(
            pseudotime_values, index=cells, name="diffusion_pseudotime"
        ),
        distance_matrix=pd.DataFrame(distances, index=cells, columns=cells),
        branches=branch_frame,
        tips=tip_frame,
        tip_indices=tip_indices,
        tip_cells=tuple(cells[index] for index in tip_indices),
        paths=paths,
        root_index=root_index,
        root_cell=cells[root_index],
    )


def infer_diffusion_root(
    potency: Iterable[float] | pd.Series,
    expression: Any,
    *,
    expression_row_names: Iterable[Any] | None = None,
    expression_col_names: Iterable[Any] | None = None,
    mean_threshold: float = 1.0,
    sd_threshold: float = 0.25,
    k_neighbors: int = 30,
    top_fraction: float = 0.05,
    n_components: int = 10,
) -> TrajectoryResult:
    """Construct a diffusion map and infer a robust high-potency root cell.

    Parameters
    ----------
    potency
        One finite continuous potency score per expression column. A pandas
        Series is accepted, but cell order—not the Series index—is used for
        alignment.
    expression
        Normalized, log-transformed gene-by-cell expression matrix supplied as
        :class:`~diffpy.LabeledMatrix`, pandas DataFrame, NumPy array, or
        SciPy sparse matrix.
    expression_row_names
        Gene identifiers for an unlabeled expression array.
    expression_col_names
        Cell identifiers for an unlabeled expression array.
    mean_threshold
        Genes must have mean expression greater than this value.
    sd_threshold
        Genes must have sample standard deviation greater than this value.
    k_neighbors
        Nearest-neighbor count used to construct the adaptive diffusion kernel.
    top_fraction
        Fraction of highest-potency cells considered as root-state candidates.
        At least ten candidates are required.
    n_components
        Maximum number of nontrivial diffusion coordinates returned.

    Returns
    -------
    TrajectoryResult
        Diffusion coordinates, transition operator, root-state metadata,
        zero-based root index, root-cell name, and potency colors.

    Raises
    ------
    ValueError
        If dimensions do not align, potency contains nonfinite values, too few
        root candidates are selected, no genes pass filtering, or neighbor and
        component parameters are invalid.
    ImportError
        If ``python-igraph`` is unavailable for walk-trap community detection.

    Examples
    --------
    >>> from diffpy import compute_ccat, infer_diffusion_root
    >>> from diffpy.datasets import load_dataset, load_ppi
    >>> liver = load_dataset("liver", data_dir="/path/to/diffpy-data")
    >>> potency = compute_ccat(liver["expression"], load_ppi("2012"))
    >>> result = infer_diffusion_root(
    ...     potency, liver["expression"], k_neighbors=30, top_fraction=0.05
    ... )
    >>> result.root_cell in liver["expression"].col_names
    True

    Notes
    -----
    Walk-trap clustering is applied only to the highest-potency candidate
    subgraph. The largest community defines the candidate root state, and the
    selected root minimizes mean absolute distance from that state's
    component-wise median. Diffusion-axis signs are arbitrary; sign flips do
    not change geometry or biological interpretation. This implementation does
    not perform an implicit PCA step.
    """

    exp = as_labeled_matrix(
        expression,
        expression_row_names,
        expression_col_names,
        require_names=True,
    )
    pot = potency.astype(float) if isinstance(potency, pd.Series) else pd.Series(potency, dtype=float)
    if len(pot) != exp.shape[1]:
        raise ValueError("potency length must equal the number of expression columns")
    if not np.all(np.isfinite(pot.to_numpy())):
        raise ValueError("potency values must all be finite")
    if not 0 < top_fraction <= 1:
        raise ValueError("top_fraction must be in (0, 1]")
    n_top = int(np.floor(top_fraction * len(pot)))
    if n_top < 10:
        raise ValueError(
            "increase top_fraction so at least 10 cells define the candidate root state"
        )
    candidate = np.argsort(-pot.to_numpy(), kind="stable")[:n_top]

    means, stds = _row_mean_sd(exp.values)
    selected = np.flatnonzero(
        (means > mean_threshold) & (stds > sd_threshold)
    )
    if selected.size == 0:
        raise ValueError("no genes pass mean_threshold and sd_threshold")
    genes_by_cells = exp.values[selected, :]
    cells_by_genes = genes_by_cells.T.tocsr() if sparse.issparse(genes_by_cells) else np.asarray(genes_by_cells).T
    transition, density, degree = _diffusion_operator(
        cells_by_genes, int(k_neighbors)
    )

    max_components = min(int(n_components) + 1, exp.shape[1] - 1)
    if max_components < 2:
        raise ValueError("At least three cells are needed for diffusion coordinates")
    eigenvalues, eigenvectors = eigsh(transition, k=max_components, which="LA")
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    coordinates = eigenvectors[:, 1:] * eigenvalues[1:]
    dc_names = tuple(f"DC{i}" for i in range(1, coordinates.shape[1] + 1))

    try:
        import igraph as ig
    except ImportError as error:
        raise ImportError("infer_diffusion_root requires python-igraph>=0.10") from error
    sub = transition[candidate, :][:, candidate].tocoo()
    edge_mask = sub.row <= sub.col
    edges = list(zip(sub.row[edge_mask].tolist(), sub.col[edge_mask].tolist()))
    graph = ig.Graph(n=n_top, edges=edges, directed=False)
    weights = sub.data[edge_mask].tolist()
    walk = graph.community_walktrap(weights=weights, steps=max(1, int(np.floor(0.25 * n_top))))
    membership = np.asarray(walk.as_clustering().membership)
    counts = np.bincount(membership)
    largest = int(np.flatnonzero(counts == counts.max())[0])
    root_state = candidate[membership == largest]

    positions = coordinates[root_state, :]
    medians = np.nanmedian(positions, axis=0)
    # Column-wise subtraction: each diffusion component uses its own median.
    distances = np.nanmean(np.abs(positions - medians[np.newaxis, :]), axis=1)
    root = int(root_state[np.nanargmin(distances)])
    cell_names = exp.col_names
    return TrajectoryResult(
        metadata={
            "eigenvalues": eigenvalues,
            "density": density,
            "degree": degree,
            "selected_gene_indices": selected,
            "root_state_indices": root_state,
        },
        diffusion_coordinates=pd.DataFrame(
            coordinates, index=cell_names, columns=dc_names
        ),
        transition_matrix=LabeledMatrix(transition, cell_names, cell_names),
        root_index=root,
        root_cell=cell_names[root],
        potency_colors=pd.Series(
            _root_colors(pot.to_numpy()), index=cell_names, name="potency_color"
        ),
    )
