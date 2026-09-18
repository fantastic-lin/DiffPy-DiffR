"""Load the curated expression, PPI network, and regulon datasets."""

from __future__ import annotations

import csv
import gzip
import json
import os
import re
from functools import lru_cache
from importlib.resources import files
from pathlib import Path
from typing import Any

import pandas as pd
from scipy import sparse

from ._matrix import LabeledMatrix


DataObject = LabeledMatrix | pd.Series | str
DataDirectory = str | os.PathLike[str] | None

@lru_cache(maxsize=1)
def _manifest() -> dict[str, Any]:
    """Read and cache packaged dataset metadata.

    Returns
    -------
    dict
        Parsed manifest indexed by public dataset name.
    """

    resource = files("diffpy").joinpath("data", "manifest.json")
    with resource.open("r", encoding="utf-8") as stream:
        return json.load(stream)


def available_datasets() -> tuple[str, ...]:
    """List all dataset identifiers supported by the package.

    Returns
    -------
    tuple of str
        Separately distributed expression cohorts plus bundled PPI networks
        and regulons, in manifest order.

    Examples
    --------
    >>> from diffpy.datasets import available_datasets
    >>> "liver" in available_datasets()
    True
    """

    return tuple(_manifest())


def _gzip_lines(resource: Any) -> tuple[str, ...]:
    """Read newline-delimited labels from a packaged gzip resource.

    Parameters
    ----------
    resource
        ``importlib.resources`` object pointing to a gzip text file.

    Returns
    -------
    tuple of str
        Lines without trailing newline characters.
    """

    with resource.open("rb") as raw, gzip.open(
        raw, "rt", encoding="utf-8"
    ) as stream:
        return tuple(line.rstrip("\n") for line in stream)


def _external_roots(data_dir: DataDirectory) -> tuple[Path, ...]:
    """Return explicit and environment-configured external data roots."""

    if data_dir is not None:
        return (Path(data_dir).expanduser(),)
    configured = os.environ.get("DIFFPY_DATA_DIR", "")
    return tuple(
        Path(item).expanduser()
        for item in configured.split(os.pathsep)
        if item
    )


def _dataset_base(dataset: str, data_dir: DataDirectory = None) -> Any:
    """Locate a dataset externally first, then inside the installed package."""

    searched: list[str] = []
    for root in _external_roots(data_dir):
        candidates = ((root,) if root.name == dataset else ()) + (root / dataset,)
        for candidate in candidates:
            searched.append(str(candidate))
            if candidate.is_dir():
                return candidate

    bundled = files("diffpy").joinpath("data", dataset)
    if bundled.is_dir():
        return bundled

    searched_text = (
        ", ".join(searched)
        if searched
        else "no external directory configured"
    )
    raise FileNotFoundError(
        f"Dataset {dataset!r} is not bundled with DiffPy and was not found "
        f"externally ({searched_text}). Extract its ZIP so that the directory "
        f"{dataset!r} is inside DIFFPY_DATA_DIR, or pass that parent directory "
        "as data_dir."
    )


def _load_matrix(base: Any, spec: dict[str, Any]) -> LabeledMatrix:
    """Load a compressed sparse matrix and its labels.

    Parameters
    ----------
    base
        Resolved package or external dataset directory.
    spec
        Matrix manifest entry containing value and label filenames.

    Returns
    -------
    LabeledMatrix
        CSR-backed matrix with preserved row and column identifiers.
    """

    with base.joinpath(spec["file"]).open("rb") as stream:
        values = sparse.load_npz(stream).tocsr()
    rows = _gzip_lines(base.joinpath(spec["rows"]))
    columns = _gzip_lines(base.joinpath(spec["columns"]))
    return LabeledMatrix(values, rows, columns)


def _load_vector(base: Any, name: str, spec: dict[str, Any]) -> pd.Series:
    """Load a compressed annotation or score vector.

    Parameters
    ----------
    base
        Resolved package or external dataset directory.
    name
        Public object name assigned to the resulting Series.
    spec
        Vector manifest entry including data type and optional factor levels.

    Returns
    -------
    pandas.Series
        Numeric, string, logical, or categorical vector with original names
        preserved as its index when available.
    """

    with base.joinpath(spec["file"]).open("rb") as raw, gzip.open(
        raw, "rt", encoding="utf-8", newline=""
    ) as stream:
        records = list(csv.DictReader(stream))
    raw_values = [record["value"] for record in records]
    dtype = spec["dtype"]
    if dtype in {"double", "integer"}:
        values: Any = pd.to_numeric(raw_values, errors="coerce")
    elif dtype == "logical":
        values = [
            value == "TRUE" if value != "NA" else pd.NA for value in raw_values
        ]
    elif dtype == "factor":
        values = pd.Categorical(
            raw_values, categories=spec.get("levels"), ordered=False
        )
    else:
        values = raw_values
    names = [record["name"] for record in records]
    index = names if any(names) else pd.RangeIndex(len(records))
    return pd.Series(values, index=index, name=name)


def load_dataset(
    name: str, *, data_dir: DataDirectory = None
) -> dict[str, DataObject]:
    """Load a complete curated dataset by its Python identifier.

    Parameters
    ----------
    name
        One value returned by :func:`available_datasets`, such as ``"chu"``,
        ``"liver"``, ``"stomach"``, ``"ppi_PC_2012"``, or
        ``"TF-regulon_network_stomach"``.
    data_dir
        Optional external data root. The root may contain a directory named
        after the dataset, or may itself be that directory. When omitted,
        directories listed in the ``DIFFPY_DATA_DIR`` environment variable
        are searched before installed package data.

    Returns
    -------
    dict
        Dictionary containing ``name``, ``description``, and the dataset's
        labeled matrices and annotation vectors. Annotation vectors correspond
        to expression columns by position; their Series indices are not
        guaranteed to contain expression cell identifiers.

    Raises
    ------
    KeyError
        If ``name`` is not present in the packaged manifest.
    FileNotFoundError
        If a separately distributed dataset is not found in an external data
        directory or in a legacy installation that bundled it.

    Examples
    --------
    >>> from diffpy.datasets import load_dataset
    >>> chu = load_dataset("chu", data_dir="/path/to/diffpy-data")
    >>> chu["expression"].is_sparse
    True
    >>> len(chu["phenotype"]) == chu["expression"].shape[1]
    True

    Notes
    -----
    Preserve annotation order when associating vector values with expression
    columns. Do not assume that an annotation Series index contains cell
    identifiers or use index-based joins without first assigning the
    expression column names.
    """

    manifest = _manifest()
    if name not in manifest:
        raise KeyError(
            f"Unknown dataset {name!r}; choose from {', '.join(manifest)}"
        )
    dataset_spec = manifest[name]
    base = _dataset_base(name, data_dir)
    objects: dict[str, DataObject] = {
        "name": name,
        "description": dataset_spec.get("description", ""),
    }
    for object_name, spec in dataset_spec["objects"].items():
        objects[object_name] = (
            _load_matrix(base, spec)
            if spec["kind"] == "matrix"
            else _load_vector(base, object_name, spec)
        )
    return objects


_REGULON_DATASETS = {
    "stomach": "TF-regulon_network_stomach",
    "skin": "TF-regulon_network_skin",
    "esophagus": "TF-regulon_network_esophagus",
    "liver": "TF-regulon_network_liver",
    "lung": "TF-regulon_network_lung",
    "pancreas": "TF-regulon_network_pancreas",
    "colon": "TF-regulon_network_colon",
    "breast": "TF-regulon_network_breast",
    "kidney": "TF-regulon_network_kidney",
}


def load_regulon(tissue: str) -> LabeledMatrix:
    """Load a tissue-specific transcription-factor regulon.

    Parameters
    ----------
    tissue
        Tissue name. Matching is case-insensitive and ignores punctuation.

    Returns
    -------
    LabeledMatrix
        Sparse gene-by-transcription-factor regulon matrix.

    Raises
    ------
    KeyError
        If the requested tissue is not one of the nine packaged regulons.

    Examples
    --------
    >>> from diffpy.datasets import load_regulon
    >>> load_regulon("stomach").shape
    (18165, 32)
    """

    key = re.sub(r"[^a-z0-9]", "", tissue.lower())
    lookup = {
        re.sub(r"[^a-z0-9]", "", tissue_name): dataset
        for tissue_name, dataset in _REGULON_DATASETS.items()
    }
    if key not in lookup:
        raise KeyError(
            f"Unknown tissue {tissue!r}; choose from "
            f"{', '.join(_REGULON_DATASETS)}"
        )
    return load_dataset(lookup[key])["regulon"]  # type: ignore[return-value]


def load_ppi(
    version: str = "2016", *, data_dir: DataDirectory = None
) -> LabeledMatrix:
    """Load a bundled or externally overridden PPI adjacency matrix.

    Parameters
    ----------
    version
        Canonical network name: ``"ppi_PC_2012"``, ``"ppi_PC_2016"``,
        ``"ppi_PC_2024"``, or ``"ppi_string_2020"``. Year aliases ``"2012"``,
        ``"2016"``, ``"2024"``, and ``"2020"`` are also accepted.
    data_dir
        Optional directory containing a canonical network folder. When
        omitted, ``DIFFPY_DATA_DIR`` is searched before the bundled network.

    Returns
    -------
    LabeledMatrix
        Sparse square gene-by-gene adjacency matrix.

    Raises
    ------
    KeyError
        If the identifier does not identify a supported network.
    FileNotFoundError
        If the selected network is absent from both external and bundled data.

    Examples
    --------
    >>> from diffpy.datasets import load_ppi
    >>> load_ppi("2012").shape
    (8434, 8434)
    """

    canonical_names = ("ppi_PC_2012", "ppi_PC_2016", "ppi_PC_2024", "ppi_string_2020")
    lookup = {name.lower(): name for name in canonical_names}
    lookup.update(dict(zip(("2012", "2016", "2024", "2020"), canonical_names)))
    key = str(version).strip().lower()
    dataset = lookup.get(key)
    if dataset is None:
        raise KeyError(f"Unknown PPI network {version!r}; choose from {', '.join(canonical_names)}")
    return load_dataset(dataset, data_dir=data_dir)["adjacency"]  # type: ignore[return-value]
