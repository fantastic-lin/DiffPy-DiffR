"""Tests for separately distributed DiffPy datasets."""

from __future__ import annotations

import gzip

from scipy import sparse

from diffpy.datasets import load_dataset, load_ppi


def test_load_ppi_from_external_root(tmp_path) -> None:
    """An extracted versioned network folder is usable through data_dir."""

    dataset = tmp_path / "ppi_2012"
    dataset.mkdir()
    sparse.save_npz(dataset / "adjacency.npz", sparse.eye(8434, format="csr"))
    labels = "".join(f"gene_{index}\n" for index in range(8434))
    for filename in ("adjacency.rows.txt.gz", "adjacency.cols.txt.gz"):
        with gzip.open(dataset / filename, "wt", encoding="utf-8") as stream:
            stream.write(labels)

    network = load_ppi("2012", data_dir=tmp_path)

    assert network.shape == (8434, 8434)
    assert network.row_names[0] == "gene_0"


def test_missing_external_dataset_has_setup_instructions(tmp_path) -> None:
    """A missing separate download produces an actionable error."""

    try:
        load_dataset("chu", data_dir=tmp_path)
    except FileNotFoundError as error:
        message = str(error)
        assert "chu" in message
        assert "DIFFPY_DATA_DIR" in message
    else:
        raise AssertionError("missing external data should raise FileNotFoundError")
