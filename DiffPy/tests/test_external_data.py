"""Tests for separately distributed DiffPy datasets."""

from __future__ import annotations

import gzip

import pytest

from scipy import sparse

from diffpy.datasets import load_dataset, load_ppi


@pytest.mark.parametrize("direct", [False, True])
def test_load_ppi_from_external_root(tmp_path, direct) -> None:
    """An extracted versioned network folder is usable through data_dir."""

    dataset = tmp_path / "ppi_PC_2016"
    dataset.mkdir()
    sparse.save_npz(dataset / "adjacency.npz", sparse.eye(11751, format="csr"))
    labels = "".join(f"gene_{index}\n" for index in range(11751))
    for filename in ("adjacency.rows.txt.gz", "adjacency.cols.txt.gz"):
        with gzip.open(dataset / filename, "wt", encoding="utf-8") as stream:
            stream.write(labels)

    network = load_ppi("2016", data_dir=dataset if direct else tmp_path)

    assert network.shape == (11751, 11751)
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


@pytest.mark.parametrize("name", ["ppi_PC_2012", "ppi_PC_2016", "ppi_PC_2024", "ppi_string_2020"])
def test_canonical_external_network_precedence(tmp_path, monkeypatch, name):
    dataset = tmp_path / name
    dataset.mkdir()
    sparse.save_npz(dataset / "adjacency.npz", sparse.csr_matrix([[0., 1.], [1., 0.]]))
    for axis in ("rows", "cols"):
        with gzip.open(dataset / f"adjacency.{axis}.txt.gz", "wt") as stream:
            stream.write("external_a\nexternal_b\n")
    monkeypatch.setenv("DIFFPY_DATA_DIR", str(tmp_path))
    network = load_ppi(name)
    assert network.row_names == ("external_a", "external_b")
    assert network.shape == (2, 2)
    assert load_dataset(name, data_dir=dataset)["adjacency"].row_names == network.row_names
