"""Single-cell potency, signaling entropy, trajectory, and regulon analysis.

The top-level namespace exposes the stable, Python-first public API. Dataset
loaders are included for convenience; advanced matrix utilities remain
available from their defining modules.
"""

from ._matrix import LabeledMatrix, as_labeled_matrix
from .ccat import compute_ccat
from .datasets import (
    available_datasets,
    load_dataset,
    load_ppi,
    load_regulon,
)
from .entropy import SignalingEntropyResult, compute_signaling_entropy
from .integration import IntegrationResult, integrate_expression_network
from .potency import PotencyStateResult, infer_potency_states
from .scira import estimate_regulatory_activity
from .trajectory import (
    DiffusionPseudotimeResult,
    TrajectoryResult,
    compute_diffusion_pseudotime,
    infer_diffusion_root,
)

__version__ = "2.0.0"

__all__ = [
    "LabeledMatrix",
    "as_labeled_matrix",
    "available_datasets",
    "load_dataset",
    "load_ppi",
    "load_regulon",
    "IntegrationResult",
    "integrate_expression_network",
    "compute_ccat",
    "SignalingEntropyResult",
    "compute_signaling_entropy",
    "PotencyStateResult",
    "infer_potency_states",
    "TrajectoryResult",
    "infer_diffusion_root",
    "DiffusionPseudotimeResult",
    "compute_diffusion_pseudotime",
    "estimate_regulatory_activity",
]
