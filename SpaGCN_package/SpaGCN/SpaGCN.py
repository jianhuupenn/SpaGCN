#cd Desktop/GCN/pygcn/env2/
import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
from scipy.sparse import issparse
from anndata import AnnData
import torch
from sklearn.decomposition import PCA
import math
import matplotlib.colors as clr
import matplotlib.pyplot as plt
from . models import *
from .cuda import CudaSimpleGCDEC as _CudaSimpleGCDEC


def _pca_embedding(adata, num_pcs):
    matrix = adata.X.toarray() if issparse(adata.X) else adata.X
    return PCA(n_components=num_pcs).fit(matrix).transform(matrix)


def _check_slice(adata, adj):
    if adj.shape != (adata.n_obs, adata.n_obs):
        raise ValueError("adj must be a square matrix in adata.obs_names order")


def _training_device(device):
    resolved = torch.device(device)
    if resolved.type == "cpu":
        return resolved
    if resolved.type != "cuda":
        raise ValueError("device must be 'cpu', 'cuda', or 'cuda:<index>'")
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is unavailable. Use CUDA-enabled PyTorch in a GPU allocation or pass device='cpu'.")
    index = torch.cuda.current_device() if resolved.index is None else resolved.index
    if index >= torch.cuda.device_count():
        raise ValueError(f"CUDA device {index} is not visible to PyTorch")
    return torch.device("cuda", index)


class SpaGCN(object):
    def __init__(self):
        super(SpaGCN, self).__init__()
        self.l=None

    def set_l(self, l):
        self.l=l

    def train(self,adata,adj, 
            num_pcs=50, 
            lr=0.005,
            max_epochs=2000,
            weight_decay=0,
            opt="admin",
            init_spa=True,
            init="louvain", #louvain or kmeans
            n_neighbors=10, #for louvain
            n_clusters=None, #for kmeans
            res=0.4, #for louvain
            tol=1e-3, *, device="cpu"):
        """Train SpaGCN using preprocessed expression and a distance matrix."""
        device = _training_device(device)
        _check_slice(adata, adj)
        if self.l is None:
            raise ValueError('l should not be set before fitting the model!')
        embed = _pca_embedding(adata, num_pcs)
        adj_exp = np.exp(-1 * (adj ** 2) / (2 * (self.l ** 2)))
        self._fit_prepared(embed, adj_exp, lr=lr, max_epochs=max_epochs,
                           weight_decay=weight_decay, opt=opt, init_spa=init_spa,
                           init=init, n_neighbors=n_neighbors, n_clusters=n_clusters,
                           res=res, tol=tol, device=device)

    def _fit_prepared(self, embed, adj_exp, *, lr, max_epochs, res, tol, device,
                      weight_decay=0, opt="admin", init_spa=True, init="louvain",
                      n_neighbors=10, n_clusters=None):
        self.num_pcs = embed.shape[1]
        self.device = device
        self.res=res
        self.lr=lr
        self.max_epochs=max_epochs
        self.weight_decay=weight_decay
        self.opt=opt
        self.init_spa=init_spa
        self.init=init
        self.n_neighbors=n_neighbors
        self.n_clusters=n_clusters
        self.tol=tol
        if device.type == "cuda":
            self.model = _CudaSimpleGCDEC(self.num_pcs, self.num_pcs, device)
        else:
            self.model = simple_GC_DEC(self.num_pcs, self.num_pcs)
        self.model.fit(embed,adj_exp,lr=self.lr,max_epochs=self.max_epochs,weight_decay=self.weight_decay,opt=self.opt,init_spa=self.init_spa,init=self.init,n_neighbors=self.n_neighbors,n_clusters=self.n_clusters,res=self.res, tol=self.tol)
        self.embed=embed
        self.adj_exp=adj_exp

    def predict(self):
        with torch.no_grad():
            _, probabilities = self.model.predict(self.embed, self.adj_exp)
            labels = probabilities.argmax(dim=1)
        return labels.cpu().numpy(), probabilities.cpu().numpy()



class multiSpaGCN(object):
    def __init__(self):
        super(multiSpaGCN, self).__init__()
        self.l=None

    def train(self,adata_list,adj_list, l_list,
            num_pcs=50, 
            lr=0.005,
            max_epochs=2000,
            weight_decay=0,
            opt="admin",
            init_spa=True,
            init="louvain", #louvain or kmeans
            n_neighbors=10, #for louvain
            n_clusters=None, #for kmeans
            res=0.4, #for louvain
            tol=1e-3):
        self.num_pcs=num_pcs
        self.res=res
        self.lr=lr
        self.max_epochs=max_epochs
        self.weight_decay=weight_decay
        self.opt=opt
        self.init_spa=init_spa
        self.init=init
        self.n_neighbors=n_neighbors
        self.n_clusters=n_clusters
        self.res=res
        self.tol=tol
        num_spots=0
        for i in adata_list: 
            num_spots+=i.shape[0]
        adj_exp_all=np.zeros((num_spots, num_spots))
        start=0
        for i in range(len(l_list)):
            l=l_list[i]
            adj=adj_list[i]
            adj_exp=np.exp(-1*(adj**2)/(2*(l**2)))
            adj_exp_all[start:start+adj_exp.shape[0],start:start+adj_exp.shape[0]]=adj_exp
            start+=adj_exp.shape[0]
        batch_cat=[str(i) for i in range(len(l_list))]
        self.adata_all=AnnData.concatenate(*adata_list,join='inner',batch_key="dataset_batch",batch_categories=batch_cat)
        pca = PCA(n_components=self.num_pcs)
        if issparse(self.adata_all.X):
            pca.fit(self.adata_all.X.A)
            embed=pca.transform(self.adata_all.X.A)
        else:
            pca.fit(self.adata_all.X)
            embed=pca.transform(self.adata_all.X)
        #----------Train model----------
        self.model=simple_GC_DEC(embed.shape[1],embed.shape[1])
        self.model.fit(embed,adj_exp_all,lr=self.lr,max_epochs=self.max_epochs,weight_decay=self.weight_decay,opt=self.opt,init_spa=self.init_spa,init=self.init,n_neighbors=self.n_neighbors,n_clusters=self.n_clusters,res=self.res, tol=self.tol)
        self.embed=embed
        self.adj_exp=adj_exp_all

    def predict(self):
        z,q=self.model.predict(self.embed,self.adj_exp)
        y_pred = torch.argmax(q, dim=1).data.cpu().numpy()
        # Max probability plot
        prob=q.detach().numpy()
        return y_pred, prob

