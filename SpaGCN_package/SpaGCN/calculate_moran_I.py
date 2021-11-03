import os,csv,re
import pandas as pd
import numpy as np
import statsmodels.api as sm
from scipy import stats
import scanpy as sc
from sklearn.neighbors import NearestNeighbors
from . calculate_adj import *
from . util import *

def Moran_I(genes_exp,x, y, k=5, knn=True):
    XYmap=pd.DataFrame({"x": x, "y":y})
    if knn:
        XYnbrs = NearestNeighbors(n_neighbors=k, algorithm='auto',metric = 'euclidean').fit(XYmap)
        XYdistances, XYindices = XYnbrs.kneighbors(XYmap)
        W = np.zeros((genes_exp.shape[0],genes_exp.shape[0]))
        for i in range(0,genes_exp.shape[0]):
            W[i,XYindices[i,:]]=1
        for i in range(0,genes_exp.shape[0]):
            W[i,i]=0
    else:
        W=calculate_adj_matrix(x=x,y=y, histology=False)
    I = pd.Series(index=genes_exp.columns, dtype="float64")
    for k in genes_exp.columns:
        X_minus_mean = np.array(genes_exp[k] - np.mean(genes_exp[k]))
        X_minus_mean = np.reshape(X_minus_mean,(len(X_minus_mean),1))
        Nom = np.sum(np.multiply(W,np.matmul(X_minus_mean,X_minus_mean.T)))
        Den = np.sum(np.multiply(X_minus_mean,X_minus_mean))
        I[k] = (len(genes_exp[k])/np.sum(W))*(Nom/Den)
    return I




def Geary_C(genes_exp,x, y, k=5, knn=True):
    XYmap=pd.DataFrame({"x": x, "y":y})
    if knn:
        XYnbrs = NearestNeighbors(n_neighbors=k, algorithm='auto',metric = 'euclidean').fit(XYmap)
        XYdistances, XYindices = XYnbrs.kneighbors(XYmap)
        W = np.zeros((genes_exp.shape[0],genes_exp.shape[0]))
        for i in range(0,genes_exp.shape[0]):
            W[i,XYindices[i,:]]=1
        for i in range(0,genes_exp.shape[0]):
            W[i,i]=0
    else:
        W=calculate_adj_matrix(x=x,y=y, histology=False)
    C = pd.Series(index=genes_exp.columns, dtype="float64")
    for k in genes_exp.columns:
        X=np.array(genes_exp[k])
        X_minus_mean = X - np.mean(X)
        X_minus_mean = np.reshape(X_minus_mean,(len(X_minus_mean),1))
        Xij=np.array([X,]*X.shape[0]).transpose()-np.array([X,]*X.shape[0])
        Nom = np.sum(np.multiply(W,np.multiply(Xij,Xij)))
        Den = np.sum(np.multiply(X_minus_mean,X_minus_mean))
        C[k] = (len(genes_exp[k])/(2*np.sum(W)))*(Nom/Den)
    return C













