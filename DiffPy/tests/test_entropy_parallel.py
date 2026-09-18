"""Regression tests for signaling-entropy parallel execution."""

from __future__ import annotations

import numpy as np
from scipy import sparse
import unittest

from diffpy import (
    IntegrationResult,
    LabeledMatrix,
    CompSRana,
)


GENES = ("g1", "g2", "g3", "g4")
CELLS = ("cell_3", "cell_1", "cell_4", "cell_2", "cell_5")
EXPRESSION = np.array(
    [
        [1.0, 2.0, 3.0, 4.0, 5.0],
        [2.0, 4.0, 1.0, 3.0, 6.0],
        [3.0, 2.0, 5.0, 1.0, 4.0],
        [4.0, 1.0, 2.0, 6.0, 3.0],
    ]
)
ADJACENCY = np.array(
    [
        [0.0, 1.0, 1.0, 0.0],
        [1.0, 0.0, 1.0, 1.0],
        [1.0, 1.0, 0.0, 1.0],
        [0.0, 1.0, 1.0, 0.0],
    ]
)


def integrated(*, sparse_expression: bool, sparse_adjacency: bool) -> IntegrationResult:
    expression = sparse.csc_matrix(EXPRESSION) if sparse_expression else EXPRESSION
    adjacency = sparse.csr_matrix(ADJACENCY) if sparse_adjacency else ADJACENCY
    return IntegrationResult(
        expression=LabeledMatrix(expression, GENES, CELLS),
        adjacency=LabeledMatrix(adjacency, GENES, GENES),
    )


def assert_results_equal(left, right) -> None:
    assert left.signaling_entropy.index.tolist() == list(CELLS)
    assert right.signaling_entropy.index.tolist() == list(CELLS)
    assert left.stationary_distribution.columns.tolist() == list(CELLS)
    assert right.stationary_distribution.columns.tolist() == list(CELLS)
    np.testing.assert_allclose(
        left.maximum_entropy, right.maximum_entropy, rtol=1e-14, atol=1e-14
    )
    np.testing.assert_allclose(
        left.signaling_entropy, right.signaling_entropy, rtol=1e-14, atol=1e-14
    )
    np.testing.assert_allclose(
        left.stationary_distribution,
        right.stationary_distribution,
        rtol=1e-14,
        atol=1e-14,
    )
    np.testing.assert_allclose(
        left.local_entropy, right.local_entropy, rtol=1e-14, atol=1e-14
    )
    np.testing.assert_allclose(
        left.normalized_local_entropy,
        right.normalized_local_entropy,
        rtol=1e-14,
        atol=1e-14,
    )


class ParallelEntropyTests(unittest.TestCase):
    def test_parallel_matches_single_worker_and_preserves_order(self) -> None:
        for sparse_expression in (False, True):
            for sparse_adjacency in (False, True):
                with self.subTest(
                    sparse_expression=sparse_expression,
                    sparse_adjacency=sparse_adjacency,
                ):
                    data = integrated(
                        sparse_expression=sparse_expression,
                        sparse_adjacency=sparse_adjacency,
                    )
                    single = CompSRana(
                        data, include_normalized_local=True, n_jobs=1
                    )
                    parallel = CompSRana(
                        data, include_normalized_local=True, n_jobs=2
                    )
                    assert_results_equal(single, parallel)

    def test_parallel_without_normalized_local_matches_single_worker(self) -> None:
        data = integrated(sparse_expression=False, sparse_adjacency=True)
        single = CompSRana(data, n_jobs=1)
        parallel = CompSRana(data, n_jobs=2)
        self.assertIsNone(single.normalized_local_entropy)
        self.assertIsNone(parallel.normalized_local_entropy)
        np.testing.assert_allclose(
            single.signaling_entropy,
            parallel.signaling_entropy,
            rtol=1e-14,
            atol=1e-14,
        )
        np.testing.assert_allclose(
            single.stationary_distribution,
            parallel.stationary_distribution,
            rtol=1e-14,
            atol=1e-14,
        )
        np.testing.assert_allclose(
            single.local_entropy, parallel.local_entropy, rtol=1e-14, atol=1e-14
        )

    def test_zero_workers_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "n_jobs cannot be zero"):
            CompSRana(
                integrated(sparse_expression=False, sparse_adjacency=False),
                n_jobs=0,
            )


if __name__ == "__main__":
    unittest.main()
