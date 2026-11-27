"""CUDA backend for the single-slice SpaGCN clustering model."""

import numpy as np
import pandas as pd
import scanpy as sc
import torch
from sklearn.cluster import KMeans

from .models import simple_GC_DEC


class CudaSimpleGCDEC(simple_GC_DEC):
    def __init__(self, n_features, n_hidden, device):
        super().__init__(n_features, n_hidden)
        self.device = device

    def fit(self, X, adj, lr=0.001, max_epochs=5000, update_interval=3,
            trajectory_interval=50, weight_decay=5e-4, opt="sgd",
            init="louvain", n_neighbors=10, res=0.4, n_clusters=10,
            init_spa=True, tol=1e-3):
        # Preserve the CPU initialization and its seeded graph convolution.
        self.cpu()
        x_cpu = torch.as_tensor(X, dtype=torch.float32)
        adj_cpu = torch.as_tensor(adj, dtype=torch.float32)
        with torch.no_grad():
            features = self.gc(x_cpu, adj_cpu).numpy()
        init_features = features if init_spa else X
        if init == "kmeans":
            labels = KMeans(n_clusters, n_init=20).fit_predict(init_features)
        elif init == "louvain":
            adata = sc.AnnData(init_features)
            sc.pp.neighbors(adata, n_neighbors=n_neighbors)
            sc.tl.louvain(adata, resolution=res)
            labels = adata.obs["louvain"].astype(int).to_numpy()
        else:
            raise ValueError(f"Unknown initializer {init!r}")

        self.n_clusters = len(np.unique(labels))
        self.trajectory = [labels]
        centers = pd.DataFrame(features).groupby(labels).mean().to_numpy()
        self.mu = torch.nn.Parameter(torch.as_tensor(centers, dtype=torch.float32))
        self.train()
        self.to(self.device)
        self._fit_input, self._fit_adj = X, adj
        self._x_cuda = x_cpu.to(self.device)
        self._adj_cuda = adj_cpu.to(self.device)

        # Upstream creates the optimizer before mu, so it updates only gc.
        if opt == "sgd":
            optimizer = torch.optim.SGD(self.gc.parameters(), lr=lr, momentum=0.9,
                                        foreach=False)
        elif opt == "admin":
            optimizer = torch.optim.Adam(self.gc.parameters(), lr=lr,
                                         weight_decay=weight_decay, foreach=False,
                                         fused=False)
        else:
            raise ValueError(f"Unknown optimizer {opt!r}")

        last_labels = torch.as_tensor(labels, device=self.device)
        for epoch in range(max_epochs):
            if epoch % update_interval == 0:
                with torch.no_grad():
                    _, q = self(self._x_cuda, self._adj_cuda)
                    target = self.target_distribution(q)
            optimizer.zero_grad()
            _, q = self(self._x_cuda, self._adj_cuda)
            self.loss_function(target, q).backward()
            optimizer.step()
            current_labels = q.detach().argmax(dim=1)
            if epoch % trajectory_interval == 0:
                self.trajectory.append(current_labels.cpu().numpy())
            if tol > 0 and epoch > 0 and (epoch - 1) % update_interval == 0:
                changes = (current_labels != last_labels).count_nonzero().item()
                if np.float32(changes) / self._x_cuda.shape[0] < tol:
                    break
            last_labels = current_labels

    def predict(self, X, adj):
        x_tensor = self._x_cuda if X is self._fit_input else torch.as_tensor(
            X, dtype=torch.float32, device=self.device)
        adj_tensor = self._adj_cuda if adj is self._fit_adj else torch.as_tensor(
            adj, dtype=torch.float32, device=self.device)
        return self(x_tensor, adj_tensor)
