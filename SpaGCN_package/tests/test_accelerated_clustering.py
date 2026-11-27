import importlib
import random
import unittest
from unittest.mock import patch

from anndata import AnnData
import numpy as np
from scipy.sparse import csr_matrix
from scipy.spatial.distance import cdist
import SpaGCN as spg
import torch


def seed():
    random.seed(100)
    np.random.seed(100)
    torch.manual_seed(100)


class AcceleratedClusteringTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        generator = np.random.RandomState(4)
        expression = generator.normal(0, 0.1, (64, 60)).astype(np.float32)
        expression[:32, :10] += 5
        expression[32:, 10:20] += 5
        cls.adata = AnnData(csr_matrix(expression))
        coordinates = np.column_stack(np.unravel_index(np.arange(64), (8, 8)))
        cls.distances = cdist(coordinates, coordinates).astype(np.float32)

    def test_sweep_reuses_pca_and_matches_individual_fits(self):
        values = [0.5, 0.9]
        sweep = importlib.import_module("SpaGCN.sweep")
        with patch.object(sweep, "_pca_embedding", wraps=sweep._pca_embedding) as pca:
            results = spg.sweep_p(self.adata, self.distances, values, 2,
                                  search_epochs=4, max_epochs=4, res_max_run=0)
        pca.assert_called_once()

        for p in values:
            length = spg.search_l(p, self.distances)
            resolution = spg.search_res(self.adata, self.distances, length, 2,
                                        start=0.7, max_epochs=4, max_run=0)
            seed()
            model = spg.SpaGCN()
            model.set_l(length)
            model.train(self.adata, self.distances, res=resolution,
                        tol=0.005, lr=0.05, max_epochs=4)
            labels, probabilities = model.predict()
            self.assertEqual(next(model.model.parameters()).device.type, "cpu")
            self.assertEqual(results[p].l, length)
            self.assertEqual(results[p].resolution, resolution)
            np.testing.assert_array_equal(results[p].labels, labels)
            np.testing.assert_allclose(results[p].probabilities, probabilities)

    @unittest.skipUnless(torch.cuda.is_available(), "CUDA is unavailable")
    def test_cuda_training_and_prediction(self):
        predictions = {}
        for device in ("cpu", "cuda"):
            seed()
            model = spg.SpaGCN()
            model.set_l(0.8)
            model.train(self.adata, self.distances, num_pcs=8,
                        init="kmeans", n_clusters=2, max_epochs=4, device=device)
            self.assertEqual(next(model.model.parameters()).device.type, device)
            predictions[device] = model.predict()
        np.testing.assert_array_equal(predictions["cpu"][0], predictions["cuda"][0])
        np.testing.assert_allclose(predictions["cpu"][1], predictions["cuda"][1],
                                   rtol=2e-5, atol=2e-5)
