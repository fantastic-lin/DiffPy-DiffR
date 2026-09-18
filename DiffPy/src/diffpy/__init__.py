"""Single-cell potency, signaling entropy, trajectory, and regulon analysis.

The top-level namespace exposes the stable, Python-first public API. Dataset
loaders are included for convenience; advanced matrix utilities remain
available from their defining modules.
"""

from ._matrix import LabeledMatrix, as_labeled_matrix
from .CompCCAT import CompCCAT
from .datasets import (
    available_datasets,
    load_dataset,
    load_ppi,
    load_regulon,
)
from .CompSRana import SignalingEntropyResult, CompSRana
from .DoIntegPPI import IntegrationResult, DoIntegPPI
from .InferPotencyStates import PotencyStateResult, InferPotencyStates
from .SciraEstRegAct import SciraEstRegAct
from .InferDMAPandRoot import (
    DiffusionPseudotimeResult,
    TrajectoryResult,
    compute_diffusion_pseudotime,
    InferDMAPandRoot,
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
    "DoIntegPPI",
    "CompCCAT",
    "SignalingEntropyResult",
    "CompSRana",
    "PotencyStateResult",
    "InferPotencyStates",
    "TrajectoryResult",
    "InferDMAPandRoot",
    "DiffusionPseudotimeResult",
    "compute_diffusion_pseudotime",
    "SciraEstRegAct",
]
