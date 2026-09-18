"""Recalculate tutorial results and render all documentation figures.

Run from the project root with the package source on ``PYTHONPATH``::

    PYTHONPATH=src python tools/generate_tutorial_plots.py --n-jobs 4

Every biological result is computed through the public ``diffpy`` API.
Matplotlib is used only to render the returned values.
"""

from __future__ import annotations

import argparse
import base64
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from diffpy import (
    CompCCAT,
    compute_diffusion_pseudotime,
    CompSRana,
    SciraEstRegAct,
    InferDMAPandRoot,
    InferPotencyStates,
    DoIntegPPI,
    load_dataset,
    load_ppi,
)


COLORS = {
    "brand": "#157b6e",
    "dark": "#0d5149",
    "accent": "#e07a38",
    "blue": "#4c78a8",
    "red": "#d1495b",
}


def _finish_figure(figure: plt.Figure, path: Path) -> None:
    """Apply common layout, save a high-resolution PNG, and close it."""

    figure.tight_layout()
    figure.savefig(path, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def _grouped_boxplot(
    values: pd.Series,
    groups: pd.Series,
    order: list[str],
    *,
    title: str,
    xlabel: str,
    ylabel: str,
    path: Path,
    colors: list[str] | None = None,
) -> None:
    """Draw a labeled boxplot with jittered cell-level observations."""

    arrays = [values.to_numpy()[groups.astype(str).to_numpy() == group] for group in order]
    figure, axis = plt.subplots(figsize=(7.2, 4.5))
    plot = axis.boxplot(arrays, tick_labels=order, patch_artist=True, showfliers=False)
    palette = colors or [COLORS["brand"]] * len(order)
    for patch, color in zip(plot["boxes"], palette, strict=True):
        patch.set_facecolor(color)
        patch.set_alpha(0.72)
    rng = np.random.default_rng(260721)
    for position, observations in enumerate(arrays, start=1):
        jitter = rng.normal(position, 0.045, size=len(observations))
        axis.scatter(jitter, observations, s=9, alpha=0.28, color=COLORS["dark"], linewidths=0)
    axis.set(title=title, xlabel=xlabel, ylabel=ylabel)
    axis.grid(axis="y", alpha=0.22)
    _finish_figure(figure, path)


def _summary_by_group(values: pd.Series, groups: pd.Series) -> dict[str, dict[str, float | int]]:
    """Return JSON-safe counts, medians, and means for each group."""

    table = pd.DataFrame({"value": values.to_numpy(), "group": groups.astype(str).to_numpy()})
    result: dict[str, dict[str, float | int]] = {}
    for group, frame in table.groupby("group", sort=False):
        result[str(group)] = {
            "n": int(len(frame)),
            "median": float(frame["value"].median()),
            "mean": float(frame["value"].mean()),
        }
    return result


def embed_images(html_path: Path) -> None:
    """Embed generated PNG files into the tutorial as Base64 data URIs.

    Each tutorial image retains its relative path in ``data-source``. This
    allows subsequent runs to replace an existing data URI with a newly
    generated image while keeping the HTML completely standalone.
    """

    html = html_path.read_text(encoding="utf-8")
    pattern = re.compile(
        r'(<img\b[^>]*\bdata-source="([^"]+)"[^>]*\bsrc=")[^"]*(")'
    )
    embedded = 0

    def replacement(match: re.Match[str]) -> str:
        nonlocal embedded
        image_path = html_path.parent / match.group(2)
        if not image_path.is_file():
            raise FileNotFoundError(f"Tutorial image not found: {image_path}")
        payload = base64.b64encode(image_path.read_bytes()).decode("ascii")
        embedded += 1
        return f'{match.group(1)}data:image/png;base64,{payload}{match.group(3)}'

    html = pattern.sub(replacement, html)
    if embedded != 8:
        raise RuntimeError(f"Expected to embed 8 tutorial images, embedded {embedded}")
    html_path.write_text(html, encoding="utf-8")


def update_captions(html_path: Path, summary: dict) -> None:
    """Keep numerical captions synchronized with the generated results."""
    sr = summary["chu_signaling_entropy"]
    ccat = summary["chu_ccat"]
    liver = summary["liver_ccat"]
    trajectory = summary["liver_trajectory"]
    dpt = summary["liver_dpt"]
    stomach = summary["stomach_regulatory_activity"]
    states = summary["chu_potency_states"]
    state_text = "; ".join(
        f"{group}: " + ", ".join(f"state {state} = {count}" for state, count in counts.items())
        for group, counts in states.items()
    )
    captions = {
        "chu_signaling_entropy": f'<strong>Calculated result.</strong> Median signaling entropy is {sr["hESC"]["median"]:.4f} in {sr["hESC"]["n"]} hESCs and {sr["EC"]["median"]:.4f} in {sr["EC"]["n"]} endothelial progenitor cells.',
        "chu_ccat": f'<strong>Calculated result.</strong> Median CCAT is {ccat["hESC"]["median"]:.4f} in hESCs and {ccat["EC"]["median"]:.4f} in endothelial progenitor cells.',
        "chu_potency_states": f'<strong>Calculated result.</strong> Potency-state assignments (state 1 is highest potency): {state_text}.',
        "liver_ccat_by_stage": f'<strong>Calculated result.</strong> Median CCAT is {liver["E10"]["median"]:.4f} at E10 and {liver["E17"]["median"]:.4f} at E17.',
        "liver_diffusion_by_stage": f'<strong>Stage view.</strong> The calculation inferred the {trajectory["root_stage"]} cell <code>{trajectory["root_cell"]}</code> as the root.',
        "liver_diffusion_pseudotime": f'<strong>DPT view.</strong> Root: <code>{dpt["root_cell"]}</code>. Pseudotime ranges from {dpt["minimum"]:.4f} to {dpt["maximum"]:.4f} (median {dpt["median"]:.4f}). The paths connect the root branch to the two differentiated branches.',
        "stomach_regulatory_activity": f'<strong>Calculated result.</strong> Across the 32 stomach-regulon transcription factors, median cell-level average activity is {stomach["undiffEpi"]["median"]:.4f} in {stomach["undiffEpi"]["n"]} undifferentiated epithelial cells and {stomach["diffEpi"]["median"]:.4f} in {stomach["diffEpi"]["n"]} differentiated epithelial cells.',
    }
    html = html_path.read_text(encoding="utf-8")
    for name, caption in captions.items():
        pattern = rf'(data-source="images/{re.escape(name)}\.png"[^>]*>\s*<figcaption>).*?(</figcaption>)'
        html, count = re.subn(pattern, lambda match: match[1] + caption + match[2], html, flags=re.DOTALL)
        if count != 1:
            raise RuntimeError(f"Expected one caption for {name}, found {count}")
    html_path.write_text(html, encoding="utf-8")


def generate(output_dir: Path, summary_path: Path, n_jobs: int) -> None:
    """Calculate all tutorial analyses and write figures plus a summary."""

    output_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.titleweight": "bold",
        }
    )
    summary: dict[str, object] = {
        "generator": "tools/generate_tutorial_plots.py",
        "random_seed": 260721,
        "ppi_network": "ppi_PC_2016",
    }

    network = load_ppi("ppi_PC_2016")
    chu = load_dataset("chu")
    integrated = DoIntegPPI(chu["expression"], network)
    entropy_result = CompSRana(integrated, n_jobs=n_jobs)
    entropy = entropy_result.signaling_entropy
    chu_groups = chu["phenotype"]
    _grouped_boxplot(
        entropy,
        chu_groups,
        ["hESC", "EC"],
        title="Signaling-entropy potency estimates",
        xlabel="Cell type",
        ylabel="Normalized signaling entropy",
        path=output_dir / "chu_signaling_entropy.png",
    )
    summary["chu_signaling_entropy"] = _summary_by_group(entropy, chu_groups)

    chu_ccat = CompCCAT(chu["expression"], network)
    _grouped_boxplot(
        chu_ccat,
        chu_groups,
        ["hESC", "EC"],
        title="CCAT potency estimates",
        xlabel="Cell type",
        ylabel="CCAT",
        path=output_dir / "chu_ccat.png",
    )
    summary["chu_ccat"] = _summary_by_group(chu_ccat, chu_groups)

    states = InferPotencyStates(
        chu_ccat,
        score_type="ccat",
        phenotype=chu_groups,
        random_state=0,
    )
    state_counts = states.distribution
    if state_counts is None:
        raise RuntimeError("phenotype-level state counts were not returned")
    figure, axis = plt.subplots(figsize=(7.2, 4.5))
    state_counts.plot(
        kind="bar",
        stacked=True,
        color=[COLORS["brand"], COLORS["accent"], COLORS["blue"], COLORS["red"]][
            : state_counts.shape[1]
        ],
        ax=axis,
    )
    axis.set(
        title="Inferred potency states by cell type",
        xlabel="Cell type",
        ylabel="Number of cells",
    )
    axis.legend(title="Potency state", frameon=False)
    axis.tick_params(axis="x", rotation=0)
    axis.grid(axis="y", alpha=0.22)
    _finish_figure(figure, output_dir / "chu_potency_states.png")
    summary["chu_potency_states"] = {
        str(group): {str(state): int(count) for state, count in row.items()}
        for group, row in state_counts.iterrows()
    }

    liver = load_dataset("liver")
    liver_ccat = CompCCAT(liver["expression"], network)
    liver_stages = pd.Series(liver["stage_codes"].index, index=liver_ccat.index)
    stage_order = ["E10", "E11", "E12", "E13", "E14", "E15", "E17"]
    stage_palette = [str(color) for color in liver["plot_colors"].to_numpy()]
    _grouped_boxplot(
        liver_ccat,
        liver_stages,
        stage_order,
        title="Liver CCAT potency across embryonic stages",
        xlabel="Embryonic stage",
        ylabel="CCAT",
        path=output_dir / "liver_ccat_by_stage.png",
        colors=stage_palette,
    )
    summary["liver_ccat"] = _summary_by_group(liver_ccat, liver_stages)

    trajectory = InferDMAPandRoot(
        liver_ccat,
        liver["expression"],
        k_neighbors=30,
        top_fraction=0.05,
    )
    coordinates = trajectory.diffusion_coordinates
    figure, axis = plt.subplots(figsize=(7.2, 5.2))
    for stage, color in zip(stage_order, stage_palette, strict=True):
        selected = liver_stages.to_numpy() == stage
        axis.scatter(
            coordinates.loc[selected, "DC1"],
            coordinates.loc[selected, "DC2"],
            s=22,
            alpha=0.78,
            color=color,
            label=stage,
            linewidths=0,
        )
    root = coordinates.loc[trajectory.root_cell]
    axis.scatter(root["DC1"], root["DC2"], marker="*", s=260, color="red", edgecolor="black", label="Root")
    axis.set(title="Liver diffusion map colored by stage", xlabel="DC1", ylabel="DC2")
    axis.legend(frameon=False, ncol=2)
    _finish_figure(figure, output_dir / "liver_diffusion_by_stage.png")

    figure, axis = plt.subplots(figsize=(7.2, 5.2))
    points = axis.scatter(
        coordinates["DC1"],
        coordinates["DC2"],
        c=liver_ccat.to_numpy(),
        cmap="viridis",
        s=22,
        alpha=0.82,
        linewidths=0,
    )
    axis.scatter(root["DC1"], root["DC2"], marker="*", s=260, color="red", edgecolor="black", label="Root")
    axis.set(title="Liver diffusion map colored by CCAT", xlabel="DC1", ylabel="DC2")
    axis.legend(frameon=False)
    figure.colorbar(points, ax=axis, label="CCAT potency")
    _finish_figure(figure, output_dir / "liver_diffusion_by_ccat.png")
    dpt = compute_diffusion_pseudotime(
        trajectory, paths_to=(1, 2, 3), window_width=0.1
    )
    figure, axis = plt.subplots(figsize=(7.2, 5.2))
    points = axis.scatter(
        coordinates["DC1"],
        coordinates["DC2"],
        c=dpt.pseudotime.to_numpy(),
        cmap="viridis",
        s=22,
        alpha=0.82,
        linewidths=0,
    )
    for (branch, path), color in zip(
        dpt.paths.items(), ("darkgreen", "green", "yellow"), strict=True
    ):
        axis.plot(
            path["DC1"],
            path["DC2"],
            color=color,
            linewidth=3,
            label=f"Path to branch {branch}",
        )
    tip_coordinates = coordinates.iloc[list(dpt.tip_indices)]
    axis.scatter(
        tip_coordinates["DC1"],
        tip_coordinates["DC2"],
        marker="D",
        s=58,
        color="red",
        edgecolor="black",
        label="DPT tips",
    )
    axis.set(
        title="Liver diffusion pseudotime and lineage paths",
        xlabel="DC1",
        ylabel="DC2",
    )
    axis.legend(frameon=False)
    figure.colorbar(points, ax=axis, label="Diffusion pseudotime")
    _finish_figure(figure, output_dir / "liver_diffusion_pseudotime.png")
    summary["liver_trajectory"] = {
        "root_cell": trajectory.root_cell,
        "root_index": int(trajectory.root_index),
        "root_stage": str(liver_stages.iloc[trajectory.root_index]),
        "selected_genes": int(len(trajectory.metadata["selected_gene_indices"])),
        "root_state_cells": int(len(trajectory.metadata["root_state_indices"])),
    }
    summary["liver_dpt"] = {
        "root_cell": dpt.root_cell,
        "tip_cells": list(dpt.tip_cells),
        "minimum": float(dpt.pseudotime.min()),
        "maximum": float(dpt.pseudotime.max()),
        "median": float(dpt.pseudotime.median()),
        "top_level_branch_counts": {
            str(branch): int(count)
            for branch, count in dpt.branches["Branch1"]
            .value_counts(dropna=False)
            .sort_index()
            .items()
        },
    }

    stomach = load_dataset("stomach")
    activity = SciraEstRegAct(
        stomach["expression"], tissue="stomach", norm="z", n_jobs=n_jobs
    )
    average_activity = activity.mean(axis=0, skipna=True)
    stomach_groups = stomach["differentiation_state"]
    _grouped_boxplot(
        average_activity,
        stomach_groups,
        ["undiffEpi", "diffEpi"],
        title="Stomach epithelial TF-regulatory activity",
        xlabel="Epithelial differentiation state",
        ylabel="Average TF activity",
        path=output_dir / "stomach_regulatory_activity.png",
    )
    summary["stomach_regulatory_activity"] = _summary_by_group(
        average_activity, stomach_groups
    )

    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    update_captions(output_dir.parent / "index.html", summary)
    embed_images(output_dir.parent / "index.html")
    from render_tutorial_outputs import render_outputs

    render_outputs(output_dir.parent / "index.html")


def main() -> None:
    """Parse command-line options and generate the tutorial artifacts."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("docs/images"))
    parser.add_argument("--summary", type=Path, default=Path("docs/generated_results.json"))
    parser.add_argument("--n-jobs", type=int, default=1)
    arguments = parser.parse_args()
    generate(arguments.output_dir, arguments.summary, arguments.n_jobs)


if __name__ == "__main__":
    main()
