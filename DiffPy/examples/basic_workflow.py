"""Run a compact potency and regulatory-activity workflow."""

from diffpy import (
    compute_ccat,
    compute_diffusion_pseudotime,
    compute_signaling_entropy,
    estimate_regulatory_activity,
    infer_diffusion_root,
    infer_potency_states,
    integrate_expression_network,
)
from diffpy.datasets import load_dataset, load_ppi


chu = load_dataset("chu")
network = load_ppi("2012")

# Network integration expects library-size-normalized, unlogged, nonnegative
# expression and always applies log2(x + 1.1). Its expression output is dense.
integrated = integrate_expression_network(chu["expression"], network)
entropy = compute_signaling_entropy(integrated, n_jobs=1)

ccat = compute_ccat(chu["expression"], network)
states = infer_potency_states(
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
liver_ccat = compute_ccat(liver["expression"], network)
trajectory = infer_diffusion_root(
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
activity = estimate_regulatory_activity(
    stomach["expression"],
    tissue="stomach",
    normalization="zscore",
    n_jobs=1,
)
average_tf_activity = activity.mean(axis=0)
print("Regulatory-activity shape:", activity.shape)
print("Average TF activity by cell:")
print(average_tf_activity.head())
