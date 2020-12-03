import scanpy as sc
import pandas as pd
import numpy as np
import scipy
import os
from anndata import AnnData,read_csv,read_text,read_mtx
from scipy.sparse import issparse

def prefilter_cells(adata,min_counts=None,max_counts=None,min_genes=200,max_genes=None):
    if min_genes is None and min_counts is None and max_genes is None and max_counts is None:
        raise ValueError('Provide one of min_counts, min_genes, max_counts or max_genes.')
    id_tmp=np.asarray([True]*adata.shape[0],dtype=bool)
    id_tmp=np.logical_and(id_tmp,sc.pp.filter_cells(adata.X,min_genes=min_genes)[0]) if min_genes is not None  else id_tmp
    id_tmp=np.logical_and(id_tmp,sc.pp.filter_cells(adata.X,max_genes=max_genes)[0]) if max_genes is not None  else id_tmp
    id_tmp=np.logical_and(id_tmp,sc.pp.filter_cells(adata.X,min_counts=min_counts)[0]) if min_counts is not None  else id_tmp
    id_tmp=np.logical_and(id_tmp,sc.pp.filter_cells(adata.X,max_counts=max_counts)[0]) if max_counts is not None  else id_tmp
    adata._inplace_subset_obs(id_tmp)
    adata.raw=sc.pp.log1p(adata,copy=True) #check the rowname 
    print("the var_names of adata.raw: adata.raw.var_names.is_unique=:",adata.raw.var_names.is_unique)
   

def prefilter_genes(adata,min_counts=None,max_counts=None,min_cells=10,max_cells=None):
    if min_cells is None and min_counts is None and max_cells is None and max_counts is None:
        raise ValueError('Provide one of min_counts, min_genes, max_counts or max_genes.')
    id_tmp=np.asarray([True]*adata.shape[1],dtype=bool)
    id_tmp=np.logical_and(id_tmp,sc.pp.filter_genes(adata.X,min_cells=min_cells)[0]) if min_cells is not None  else id_tmp
    id_tmp=np.logical_and(id_tmp,sc.pp.filter_genes(adata.X,max_cells=max_cells)[0]) if max_cells is not None  else id_tmp
    id_tmp=np.logical_and(id_tmp,sc.pp.filter_genes(adata.X,min_counts=min_counts)[0]) if min_counts is not None  else id_tmp
    id_tmp=np.logical_and(id_tmp,sc.pp.filter_genes(adata.X,max_counts=max_counts)[0]) if max_counts is not None  else id_tmp
    adata._inplace_subset_var(id_tmp)


def prefilter_specialgenes(adata,Gene1Pattern="ERCC",Gene2Pattern="MT-"):
    id_tmp1=np.asarray([not str(name).startswith(Gene1Pattern) for name in adata.var_names],dtype=bool)
    id_tmp2=np.asarray([not str(name).startswith(Gene2Pattern) for name in adata.var_names],dtype=bool)
    id_tmp=np.logical_and(id_tmp1,id_tmp2)
    adata._inplace_subset_var(id_tmp)


def test_l(adj, list_l):
    for l in list_l:
        adj_exp=np.exp(-1*adj/(2*(l**2)))
        print("l is ",str(l),"Percentage of total expression contributed by neighborhoods:",np.mean(np.sum(adj_exp,1))-1)

def find_l(p, adj, start=0.5, end=2,sep=0.01, tol=0.01):
    for i in np.arange(start, end, sep):
        adj_exp=np.exp(-1*adj/(2*(i**2)))
        q=np.mean(np.sum(adj_exp,1))-1
        print("L=", str(i), "P=", str(round(q,5)))
        if np.abs(p-q)<tol:
            return i
    print("l not found, try bigger range or smaller sep!")


def find_neighbor_clusters(target_cluster,cell_id, x, y, pred,radius, ratio=1/2):
    cluster_num = dict()
    for i in pred:
        cluster_num[i] = cluster_num.get(i, 0) + 1
    df = {'cell_id': cell_id, 'x': x, "y":y, "pred":pred}
    df = pd.DataFrame(data=df)
    df.index=df['cell_id']
    target_df=df[df["pred"]==target_cluster]
    nbr_num={}
    row_index=0
    num_nbr=[]
    for index, row in target_df.iterrows():
        x=row["x"]
        y=row["y"]
        tmp_nbr=df[((df["x"]-x)**2+(df["y"]-y)**2)<=(radius**2)]
        #tmp_nbr=df[(df["x"]<x+radius) & (df["x"]>x-radius) & (df["y"]<y+radius) & (df["y"]>y-radius)]
        num_nbr.append(tmp_nbr.shape[0])
        for p in tmp_nbr["pred"]:
            nbr_num[p]=nbr_num.get(p,0)+1
    del nbr_num[target_cluster]
    nbr_num=[(k, v)  for k, v in nbr_num.items() if v>(ratio*cluster_num[k])]
    nbr_num.sort(key=lambda x: -x[1])
    print("radius=", radius, "average number of neighbors for each spot is", np.mean(num_nbr))
    print(" Cluster",target_cluster, "has neighbors:")
    for t in nbr_num:
        print("Dmain ", t[0], ": ",t[1])
    ret=[t[0] for t in nbr_num]
    if len(ret)==0:
        print("No neighbor domain found, try bigger radius or smaller ratio")
    else:
        return ret


def rank_genes_groups(input_adata, target_cluster,nbr_list, label_col, adj_nbr=True, log=False):
    if adj_nbr:
        nbr_list=nbr_list+[target_cluster]
        adata=input_adata[input_adata.obs[label_col].isin(nbr_list)]
    else:
        adata=input_adata.copy()
    adata.var_names_make_unique()
    adata.obs["target"]=((adata.obs[label_col]==target_cluster)*1).astype('category')
    sc.tl.rank_genes_groups(adata, groupby="target",reference="rest", n_genes=adata.shape[1],method='wilcoxon')
    pvals_adj=[i[0] for i in adata.uns['rank_genes_groups']["pvals_adj"]]
    genes=[i[1] for i in adata.uns['rank_genes_groups']["names"]]
    if issparse(adata.X):
        obs_tidy=pd.DataFrame(adata.X.A)
    else:
        obs_tidy=pd.DataFrame(adata.X)
    obs_tidy.index=adata.obs["target"].tolist()
    obs_tidy.columns=adata.var.index.tolist()
    obs_tidy=obs_tidy.loc[:,genes]
    # 1. compute mean value
    mean_obs = obs_tidy.groupby(level=0).mean()
    # 2. compute fraction of cells having value >0
    obs_bool = obs_tidy.astype(bool)
    fraction_obs = obs_bool.groupby(level=0).sum() / obs_bool.groupby(level=0).count()
    # compute fold change.
    if log: #The adata already logged
        fold_change=np.exp((mean_obs.loc[1] - mean_obs.loc[0]).values)
    else:
        fold_change = (mean_obs.loc[1] / (mean_obs.loc[0]+ 1e-9)).values
    df = {'genes': genes, 'in_group_fraction': fraction_obs.loc[1].tolist(), "out_group_fraction":fraction_obs.loc[0].tolist(),"in_out_group_ratio":(fraction_obs.loc[1]/fraction_obs.loc[0]).tolist(),"in_group_mean_exp": mean_obs.loc[1].tolist(), "out_group_mean_exp": mean_obs.loc[0].tolist(),"fold_change":fold_change.tolist(), "pvals_adj":pvals_adj}
    df = pd.DataFrame(data=df)
    return df

def relative_func(expres):
    #expres: an array counts expression for a gene
    maxd = np.max(expres) - np.min(expres)
    min_exp=np.min(expres)
    rexpr = (expres - min_exp)/maxd
    return rexpr

def plot_relative_exp(input_adata, gene, x_name, y_name,color,use_raw=False, spot_size=200000):
    adata=input_adata.copy()
    if use_raw:
        X=adata.X
    else:
        X=adata.raw.X
    if issparse(X):
        X=pd.DataFrame(X.A)
    else:
        X=pd.DataFrame(X)
    X.index=adata.obs.index
    X.columns=adata.var.index
    rexpr=relative_func(X.loc[:,gene])
    adata.obs["rexpr"]=rexpr
    fig=sc.pl.scatter(adata,x=x_name,y=y_name,color="rexpr",title=gene+"_rexpr",color_map=color,show=False,size=spot_size/adata.shape[0])
    return fig

def plot_log_exp(input_adata, gene, x_name, y_name,color,use_raw=False):
    adata=input_adata.copy()
    if use_raw:
        X=adata.X
    else:
        X=adata.raw.X
    if issparse(X):
        X=pd.DataFrame(X.A)
    else:
        X=pd.DataFrame(X)
    X.index=adata.obs.index
    X.columns=adata.var.index
    adata.obs["log"]=np.log((X.loc[:,gene]+1).tolist())
    fig=sc.pl.scatter(adata,x=x_name,y=y_name,color="log",title=gene+"_log",color_map=color,show=False,size=200000/adata.shape[0])
    return fig

def detect_subclusters(cell_id, x, y, pred, target_cluster, radius=3, res=0.2):
    df = {'cell_id': cell_id, 'x': x, "y":y, "pred":pred}
    df = pd.DataFrame(data=df)
    df.index=df['cell_id']
    target_df=df[df["pred"]==target_cluster]
    nbr=np.zeros([target_df.shape[0],len(set(df["pred"]))],dtype=int)
    num_nbr=[]
    row_index=0
    for index, row in target_df.iterrows():
        x=row["x"]
        y=row["y"]
        tmp_nbr=df[(df["x"]<x+radius) & (df["x"]>x-radius) & (df["y"]<y+radius) & (df["y"]>y-radius)]
        num_nbr.append(tmp_nbr.shape[0])
        for p in tmp_nbr["pred"]:
            nbr[row_index,int(p)]+=1
        row_index+=1
    #Minus out the cell itself
    nbr[:,target_cluster]=nbr[:,target_cluster]-1
    nbr=sc.AnnData(nbr)
    sc.pp.neighbors(nbr, n_neighbors=10)
    sc.tl.louvain(nbr,resolution=res)
    sub_cluster=nbr.obs['louvain'].astype(int).to_numpy()
    target_df["sub_cluster"]=sub_cluster
    target_df["sub_cluster"]=target_df["sub_cluster"].astype('category')
    tmp=[]
    for j in df.index:
        if j in target_df.index:
            tmp.append(target_df.loc[j,"sub_cluster"])
        else:
            tmp.append("-1")
    #ret = {'cell_id': cell_id, 'sub_cluster_'+str(target_cluster): tmp}
    #ret = pd.DataFrame(data=ret)
    #ret.index=ret['cell_id']
    ret=tmp
    return ret

def find_meta_gene(input_adata,
                    pred,
                    target_domain,
                    start_gene,
                    mean_diff=0,
                    early_stop=True,
                    max_iter=5,
                    use_raw=False):
    meta_name=start_gene
    adata=input_adata.copy()
    adata.obs["meta"]=adata.X[:,adata.var.index==start_gene]
    adata.obs["pred"]=pred
    num_non_target=adata.shape[0]
    for i in range(max_iter):
        #Select cells
        tmp=adata[((adata.obs["meta"]>np.mean(adata.obs[adata.obs["pred"]==target_domain]["meta"]))|(adata.obs["pred"]==target_domain))]
        tmp.obs["target"]=((tmp.obs["pred"]==target_domain)*1).astype('category').copy()
        #DE
        sc.tl.rank_genes_groups(tmp, groupby="target",reference="rest", n_genes=1,method='wilcoxon')
        adj_g=tmp.uns['rank_genes_groups']["names"][0][0]
        add_g=tmp.uns['rank_genes_groups']["names"][0][1]
        meta_name_cur=meta_name+"+"+add_g+"-"+adj_g
        print("Add gene: ", add_g)
        print("Minus gene: ", adj_g)
        #Meta gene
        adata.obs[add_g]=adata.X[:,adata.var.index==add_g]
        adata.obs[adj_g]=adata.X[:,adata.var.index==adj_g]
        adata.obs["meta_cur"]=(adata.obs["meta"]+adata.obs[add_g]-adata.obs[adj_g])
        adata.obs["meta_cur"]=adata.obs["meta_cur"]-np.min(adata.obs["meta_cur"])
        mean_diff_cur=np.mean(adata.obs["meta_cur"][adata.obs["pred"]==target_domain])-np.mean(adata.obs["meta_cur"][adata.obs["pred"]!=target_domain])
        num_non_target_cur=np.sum(tmp.obs["target"]==0)
        if (early_stop==False) | ((num_non_target>=num_non_target_cur) & (mean_diff<=mean_diff_cur)):
            num_non_target=num_non_target_cur
            mean_diff=mean_diff_cur
            print("Absolute mean change:", mean_diff)
            print("Number of non-target spots reduced to:",num_non_target)
        else:
            print("Stopped!", "Previous Number of non-target spots",num_non_target, num_non_target_cur, mean_diff,mean_diff_cur)
            print("Previous Number of non-target spots",num_non_target, num_non_target_cur, mean_diff,mean_diff_cur)
            print("Previous Number of non-target spots",num_non_target)
            print("Current Number of non-target spots",num_non_target_cur)
            print("Absolute mean change", mean_diff)
            print("===========================================================================")
            print("Meta gene: ", meta_name)
            print("===========================================================================")
            return meta_name, adata.obs["meta"].tolist()
        meta_name=meta_name_cur
        adata.obs["meta"]=adata.obs["meta_cur"]
        print("===========================================================================")
        print("Meta gene is: ", meta_name)
        print("===========================================================================")
    
    return meta_name, adata.obs["meta"].tolist()




