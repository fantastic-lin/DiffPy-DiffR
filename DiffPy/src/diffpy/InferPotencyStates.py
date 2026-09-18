"""Infer discrete cellular potency states with Gaussian mixture models."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture


@dataclass(frozen=True)
class PotencyStateResult:
    """Discrete potency assignments and phenotype-level summaries.

    Parameters
    ----------
    states
        One-indexed state assignment for every cell. State 1 has the greatest
        fitted mean potency.
    distribution
        Optional phenotype-by-state cell counts.
    probabilities
        Optional within-phenotype state probabilities.
    heterogeneity
        Optional normalized Shannon entropy of potency-state probabilities for
        each phenotype.
    model
        BIC-selected scikit-learn Gaussian mixture model.
    transformed_values
        Values used to fit the mixture after log-odds transformation.
    """

    states: pd.Series
    distribution: pd.DataFrame | None
    probabilities: pd.DataFrame | None
    heterogeneity: pd.Series | None
    model: GaussianMixture
    transformed_values: pd.Series

    def as_dict(self) -> dict[str, Any]:
        """Return the biological result tables in a plain dictionary.

        Returns
        -------
        dict
            State assignments, distribution, probabilities, and heterogeneity.
            The fitted model and transformed values remain available as object
            attributes.
        """

        return {
            "states": self.states,
            "distribution": self.distribution,
            "probabilities": self.probabilities,
            "heterogeneity": self.heterogeneity,
        }


def InferPotencyStates(
    potency: Iterable[float] | pd.Series,
    *,
    score_type: str = "signaling_entropy",
    phenotype: Iterable[Any] | pd.Series | None = None,
    variable_variance: bool = True,
    max_states: int = 5,
    random_state: int | None = 0,
) -> PotencyStateResult:
    """Infer ordered potency states from continuous cellular scores.

    Parameters
    ----------
    potency
        One finite potency estimate per cell. A pandas Series preserves cell
        identifiers in every returned table.
    score_type
        Transformation appropriate for the supplied scores. Use
        ``"signaling_entropy"`` for values strictly between 0 and 1, or
        ``"ccat"`` for correlation scores strictly between -1 and 1.
    phenotype
        Optional categorical annotation with one value per cell. When given,
        the result includes phenotype-by-state counts, probabilities, and
        normalized Shannon heterogeneity.
    variable_variance
        Compare mixtures with component-specific and shared variances. When
        false, only shared-variance models are considered.
    max_states
        Maximum mixture-component count evaluated by Bayesian information
        criterion (BIC).
    random_state
        Seed controlling mixture initialization. The default makes repeated
        calls reproducible.

    Returns
    -------
    PotencyStateResult
        Ordered state assignments, optional phenotype summaries, the selected
        model, and transformed fitting values.

    Raises
    ------
    ValueError
        If scores are nonfinite, outside the domain required by ``score_type``,
        the phenotype length differs, or ``max_states`` is invalid.

    Examples
    --------
    >>> from diffpy import InferPotencyStates
    >>> from diffpy.datasets import load_dataset
    >>> chu = load_dataset("chu", data_dir="/path/to/diffpy-data")
    >>> result = InferPotencyStates(
    ...     chu["signaling_entropy"],
    ...     score_type="signaling_entropy",
    ...     phenotype=chu["phenotype"],
    ... )
    >>> result.states.min()
    1

    Notes
    -----
    Signaling entropy is transformed as ``log2(s / (1 - s))``. CCAT is
    transformed as ``log2((1 + c) / (1 - c))``. Component labels are reordered
    after fitting so lower state numbers always represent greater potency.
    """

    series = (
        potency.astype(float)
        if isinstance(potency, pd.Series)
        else pd.Series(potency, dtype=float)
    )
    kind = score_type.lower().replace("-", "_").replace(" ", "_")
    values = series.to_numpy(dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError("potency values must all be finite")
    if kind in {"signaling_entropy", "entropy", "sr"}:
        if np.any((values <= 0) | (values >= 1)):
            raise ValueError(
                "signaling-entropy values must lie strictly between 0 and 1"
            )
        transformed_values = np.log2(values / (1.0 - values))
        canonical_kind = "signaling_entropy"
    elif kind == "ccat":
        if np.any((values <= -1) | (values >= 1)):
            raise ValueError("CCAT values must lie strictly between -1 and 1")
        transformed_values = np.log2((1.0 + values) / (1.0 - values))
        canonical_kind = "ccat"
    else:
        raise ValueError("score_type must be 'signaling_entropy' or 'ccat'")

    data = transformed_values.reshape(-1, 1)
    covariance_types = ("full", "tied") if variable_variance else ("tied",)
    upper = min(int(max_states), len(values))
    if upper < 1:
        raise ValueError("max_states must be at least one")
    candidates: list[GaussianMixture] = []
    for state_count in range(1, upper + 1):
        for covariance_type in covariance_types:
            candidates.append(
                GaussianMixture(
                    n_components=state_count,
                    covariance_type=covariance_type,
                    random_state=random_state,
                    n_init=10,
                ).fit(data)
            )
    model = min(candidates, key=lambda fitted: fitted.bic(data))
    labels = model.predict(data)
    high_to_low = np.argsort(-model.means_.ravel())
    rank = np.empty_like(high_to_low)
    rank[high_to_low] = np.arange(1, len(high_to_low) + 1)
    states = pd.Series(rank[labels], index=series.index, name="potency_state")

    distribution = probabilities = None
    heterogeneity = None
    if phenotype is not None:
        phenotype_series = (
            phenotype
            if isinstance(phenotype, pd.Series)
            else pd.Series(phenotype, index=series.index)
        )
        if len(phenotype_series) != len(series):
            raise ValueError("phenotype and potency must have equal length")
        phenotype_series = pd.Series(
            phenotype_series.to_numpy(), index=series.index, name="phenotype"
        )
        distribution = pd.crosstab(
            phenotype_series, states, dropna=False
        ).reindex(columns=range(1, model.n_components + 1), fill_value=0)
        probabilities = distribution.div(distribution.sum(axis=1), axis=0)
        if model.n_components == 1:
            heterogeneity = pd.Series(
                0.0, index=probabilities.index, name="heterogeneity"
            )
        else:
            safe_probabilities = probabilities.where(probabilities > 0, 1.0)
            heterogeneity = -(
                probabilities * np.log(safe_probabilities)
            ).sum(axis=1) / np.log(model.n_components)
            heterogeneity.name = "heterogeneity"
    transformed = pd.Series(
        transformed_values,
        index=series.index,
        name=f"transformed_{canonical_kind}",
    )
    return PotencyStateResult(
        states=states,
        distribution=distribution,
        probabilities=probabilities,
        heterogeneity=heterogeneity,
        model=model,
        transformed_values=transformed,
    )
