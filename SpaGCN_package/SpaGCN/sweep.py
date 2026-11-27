"""Run a SpaGCN p sweep with one PCA fit per tissue slice."""

import math
import random
from typing import NamedTuple

import numpy as np
import torch

from .SpaGCN import SpaGCN, _check_slice, _pca_embedding, _training_device
from .util import _search_l, _search_resolution, calculate_p

__all__ = ["SweepResult", "sweep_p"]


class SweepResult(NamedTuple):
    """Predictions for one p value, in the input spot order."""

    l: float
    resolution: float
    labels: np.ndarray
    probabilities: np.ndarray


def sweep_p(adata, adj, p_values, n_clusters, *, device="cpu", seed=100,
            num_pcs=50, search_epochs=20, max_epochs=200, lr=0.05, tol=5e-3,
            res_start=0.7, res_step=0.1, res_max_run=10, l_start=0.01,
            l_end=1000, l_tol=0.01, l_max_run=100):
    """Cluster one preprocessed slice for each p, reusing its PCA and distances.

    Supply the dense distances from calculate_adj_matrix in adata.obs_names
    order. Returns a dict keyed by p. As in search_res, n_clusters is a target.
    """
    device = _training_device(device)
    _check_slice(adata, adj)
    values = tuple(float(p) for p in p_values)
    if any(not math.isfinite(p) or p < 0 for p in values):
        raise ValueError("p_values must contain finite, nonnegative values")
    if not values:
        return {}
    if len(set(values)) != len(values):
        raise ValueError("p_values must be unique")

    np.random.seed(seed)
    embed = _pca_embedding(adata, num_pcs)
    # The upstream seeded PCA may consume random numbers before cluster fitting.
    after_pca = np.random.get_state()
    evaluated_lengths = {}

    def evaluate_length(length):
        if length not in evaluated_lengths:
            evaluated_lengths[length] = calculate_p(adj, length)
        return evaluated_lengths[length]

    def fit_predict(length, affinity, resolution, epochs):
        random.seed(seed)
        torch.manual_seed(seed)
        np.random.set_state(after_pca)
        model = SpaGCN()
        model.set_l(length)
        model._fit_prepared(embed, affinity, lr=lr, max_epochs=epochs,
                            res=resolution, tol=tol, device=device)
        return model.predict()

    results = {}
    for p in values:
        length = _search_l(p, evaluate_length, l_start, l_end, l_tol, l_max_run)
        if length is None:
            raise ValueError(f"No length scale found for p={p}. Check l_start, l_end and l_tol.")
        affinity = np.exp(-1 * (adj ** 2) / (2 * (length ** 2)))

        def count_clusters(resolution):
            labels, _ = fit_predict(length, affinity, resolution, search_epochs)
            return len(set(labels))

        resolution = _search_resolution(n_clusters, count_clusters, res_start,
                                        res_step, res_max_run)
        labels, probabilities = fit_predict(length, affinity, resolution, max_epochs)
        results[p] = SweepResult(length, resolution, labels, probabilities)
    return results
