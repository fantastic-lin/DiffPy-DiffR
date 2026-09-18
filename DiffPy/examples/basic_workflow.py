"""Run a compact potency and regulatory-activity workflow."""

from diffpy import (
    CompCCAT,
    compute_diffusion_pseudotime,
    CompSRana,
    SciraEstRegAct,
    InferDMAPandRoot,
    InferPotencyStates,
    DoIntegPPI,
)
from diffpy.datasets import load_dataset, load_ppi


chu = load_dataset("chu")
network = load_ppi("ppi_PC_2016")

# Network integration expects library-size-normalized, unlogged, nonnegative
# expression and always applies log2(x + 1.1). Its expression output is dense.
integrated = DoIntegPPI(chu["expression"], network)
entropy = CompSRana(integrated, n_jobs=1)

ccat = CompCCAT(chu["expression"], network)
states = InferPotencyStates(
    ccat,
    score_type="ccat",
    phenotype=chu["phenotype"],
)

print("CCAT range:", float(ccat.min()), float(ccat.max()))
print(
    "Signaling-entropy range:",
    float(entropy.signaling_entropy.min()),
    float(entropy.signaling_entropy.max()),
)
print("Potency-state distribution:")
print(states.distribution)

liver = load_dataset("liver")
liver_ccat = CompCCAT(liver["expression"], network)
trajectory = InferDMAPandRoot(
    liver_ccat,
    liver["expression"],
    k_neighbors=30,
    top_fraction=0.05,
)
dpt = compute_diffusion_pseudotime(
    trajectory,
    paths_to=(1, 2, 3),
)

print("Inferred diffusion root:", trajectory.root_cell)
print("Diffusion-coordinate shape:", trajectory.diffusion_coordinates.shape)
print(
    "Diffusion-pseudotime range:",
    float(dpt.pseudotime.min()),
    float(dpt.pseudotime.max()),
)
print("Diffusion branch tips:", dpt.tip_cells)

stomach = load_dataset("stomach")
activity = SciraEstRegAct(
    stomach["expression"],
    tissue="stomach",
    norm="z",
    n_jobs=1,
)
average_tf_activity = activity.mean(axis=0)
print("Regulatory-activity shape:", activity.shape)
print("Average TF activity by cell:")
print(average_tf_activity.head())