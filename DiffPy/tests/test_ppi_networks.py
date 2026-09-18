"""Bundled network registration, compatibility, and preprocessing regressions."""
import json
from importlib.resources import files
import warnings

import numpy as np
import pytest
from scipy import sparse

from diffpy import DoIntegPPI, LabeledMatrix
from diffpy.datasets import available_datasets, load_dataset, load_ppi

NETWORKS = ("ppi_PC_2012", "ppi_PC_2016", "ppi_PC_2024", "ppi_string_2020")


@pytest.mark.parametrize("name,alias", list(zip(NETWORKS, ("2012", "2016", "2024", "2020"))))
def test_bundled_network(name, alias, monkeypatch):
    monkeypatch.delenv("DIFFPY_DATA_DIR", raising=False)
    network = load_ppi(name)
    by_alias = load_ppi(alias)
    assert name in available_datasets()
    assert network.row_names == network.col_names == by_alias.row_names
    assert len(set(network.row_names)) == network.shape[0]
    assert (network.values != by_alias.values).nnz == 0
    assert (network.values != network.values.T).nnz == 0
    assert np.all(network.values.data == 1)
    assert not np.any(network.values.diagonal())
    manifest = json.loads(files("diffpy").joinpath("data/manifest.json").read_text())
    spec = manifest[name]["objects"]["adjacency"]
    assert list(network.shape) == spec["shape"]
    assert network.values.nnz == spec["nnz"]
    by_dataset = load_dataset(name)["adjacency"]
    assert by_dataset.row_names == network.row_names
    assert (by_dataset.values != network.values).nnz == 0


@pytest.mark.parametrize("year", ["2012", "2016"])
def test_removed_legacy_names_are_rejected(year):
    removed_name = "ppi_" + year
    with pytest.raises(KeyError, match="Unknown dataset"):
        load_dataset(removed_name)
    with pytest.raises(KeyError, match="Unknown PPI network"):
        load_ppi(removed_name)


def test_unknown_network():
    with pytest.raises(KeyError, match="Unknown PPI network"):
        load_ppi("unknown")


@pytest.mark.parametrize("maximum", [10, 100, 101])
@pytest.mark.parametrize("sparse_input", [False, True])
def test_integration_preprocessing_is_unchanged(maximum, sparse_input):
    genes = ("a", "b", "c")
    values = np.array([[0., maximum], [1., 2.], [3., 4.]])
    expression = LabeledMatrix(sparse.csr_matrix(values) if sparse_input else values, genes, ("x", "y"))
    network = LabeledMatrix(sparse.csr_matrix(np.ones((3, 3)) - np.eye(3)), genes, genes)
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        result = DoIntegPPI(expression, network, min_overlap=3)
    assert not sparse.issparse(result.expression.values)
    np.testing.assert_allclose(result.expression.values, np.log2(values + 1.1))
    assert any("below 100" in str(w.message) for w in captured) == (maximum < 100)
