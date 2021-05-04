#--------------------------------------------------------------------------------------------------------#
#-------------------------------1. Load packages--------------------------------------------------#
#--------------------------------------------------------------------------------------------------------#
import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc

import math
#pip3 install scikit-image
from skimage import io, color
#conda config --add channels bioconda
#conda install -n pytorch scanpy
"""
pip3 uninstall SpaGCN
conda activate pytorch
python3 setup.py build
python3 setup.py install
conda deactivate
"""
import SpaGCN as spg

#--------------------------------------------------------------------------------------------------------#
#-------------------------------2. Read in data--------------------------------------------------#
#--------------------------------------------------------------------------------------------------------#
"""
#Read original data and save it to h5ad
from scanpy import read_10x_h5
adata = read_10x_h5("../tutorial/data/expression_matrix.h5")
spatial=pd.read_csv("../tutorial/data/positions.txt",sep=",",header=None,na_filter=False,index_col=0) 
adata.obs["x1"]=spatial[1]
adata.obs["x2"]=spatial[2]
adata.obs["x3"]=spatial[3]
adata.obs["x4"]=spatial[4]
adata.obs["x5"]=spatial[5]
#Select captured samples
adata=adata[adata.obs["x1"]==1]
adata.var_names=[i.upper() for i in list(adata.var_names)]
adata.var["genename"]=adata.var.index.astype("str")
adata.write_h5ad("../tutorial/data/sample_data.h5ad")
"""
#Read in gene expression and spatial location
adata=sc.read("../tutorial/data/sample_data.h5ad")
#Read in hitology image
image=io.imread("../tutorial/data/histology.tif")
#--------------------------------------------------------------------------------------------------------#
#-------------------------------3. Calculate adjacent matrix--------------------------------------------------#
#--------------------------------------------------------------------------------------------------------#
b=49
a=1
x2=adata.obs["x2"].tolist()
x3=adata.obs["x3"].tolist()
x4=adata.obs["x4"].tolist()
x5=adata.obs["x5"].tolist()
adj=spg.calculate_adj_matrix(x=x2,y=x3, x_pixel=x4, y_pixel=x5, image=image, beta=b, alpha=a, histology=True)
np.savetxt('../tutorial/data/adj.csv', adj, delimiter=',')


#--------------------------------------------------------------------------------------------------------#
#-------------------------------------4. Run SpaGCN-----------------------------------------------------#
#--------------------------------------------------------------------------------------------------------#
adata.var_names_make_unique()
spg.prefilter_genes(adata,min_cells=3) # avoiding all genes are zeros
spg.prefilter_specialgenes(adata)
#Normalize and take log for UMI-------
sc.pp.normalize_per_cell(adata, min_adata=0)
sc.pp.log1p(adata)
#Set percentage of total expression contributed by neighborhoods
p=0.5
#Find l
#spg.test_l(adj,[0.5,0.8,1,1.2])
"""
l is  0.5 Percentage of total expression contributed by neighborhoods: 0.022694894646745123
l is  0.8 Percentage of total expression contributed by neighborhoods: 0.5607284420349321
l is  1 Percentage of total expression contributed by neighborhoods: 2.0466495056643024
l is  1.2 Percentage of total expression contributed by neighborhoods: 5.558170663311434
"""
#Therefore,search l around 0.8
l=spg.find_l(p=p,adj=adj,start=0.75, end=0.8,sep=0.001, tol=0.01)
res=0.6
clf=spg.SpaGCN()
clf.set_l(l)
#Init using louvain
clf.train(adata,adj,init_spa=True,init="louvain",res=res, tol=5e-3)
#Or init using kmeans
#clf.train(adata,adj,init_spa=True,init="kmeans",n_clusters=7, tol=5e-3)
y_pred, prob=clf.predict()

adata.obs["pred"]= y_pred
adata.obs["pred"]=adata.obs["pred"].astype('category')
#Save
adata.write_h5ad("../tutorial/results.h5ad")

#Plot
colors_use=['#111010', '#FFFF00', '#4a6fe3', '#bb7784', '#bec1d4', '#ff9896', '#98df8a', '#ffbb78', '#2ca02c', '#ff7f0e', '#1f77b4', '#800080', '#959595', '#ffff00', '#014d01', '#0000ff', '#ff0000', '#000000']
num_celltype=len(adata.obs["pred"].unique())
adata.uns["pred_colors"]=list(colors_use[:num_celltype])
fig=sc.pl.scatter(adata,alpha=1,x="x5",y="x4",color="pred",show=False,size=150000/adata.shape[0])
fig.set_aspect('equal', 'box')
fig.figure.savefig("../tutorial/domains.png", dpi=300)

#--------------------------------------------------------------------------------------------------------#
#-------------------------------------5. Identify SVGs -------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------#
#Use domain 1 as an example
min_in_group_fraction=0.8
min_in_out_group_ratio=1
min_fold_change=1.2
target=1
nbr_domians=spg.find_neighbor_clusters(target_cluster=target,
                                   cell_id=adata.obs.index.tolist(), 
                                   x=adata.obs["x2"].tolist(), 
                                   y=adata.obs["x3"].tolist(), 
                                   pred=adata.obs["pred"].tolist(),
                                   radius=2,
                                   ratio=1/2)
nbr_domians=nbr_domians[0:3]
de_genes_info=spg.rank_genes_groups(input_adata=adata,
                                target_cluster=target,
                                nbr_list=nbr_domians, 
                                label_col="pred", 
                                adj_nbr=True, 
                                log=True)
de_genes_info=de_genes_info[(de_genes_info["pvals_adj"]<0.05)]
#de_genes_info=de_genes_info.sort_values(by="in_group_fraction", ascending=False)
filtered_info=de_genes_info
filtered_info=filtered_info[(filtered_info["pvals_adj"]<0.05) &
(filtered_info["in_out_group_ratio"]>=min_in_out_group_ratio) &
(filtered_info["in_group_fraction"]>=min_in_group_fraction) &
(filtered_info["fold_change"]>=min_fold_change)]
filtered_info=filtered_info.sort_values(by="in_group_fraction", ascending=False)
SVG=filtered_info["genes"].tolist()

print("Domain", target,":", SVG)

#Plot SVG
import matplotlib.colors as clr
import matplotlib.pyplot as plt
color_self = clr.LinearSegmentedColormap.from_list('pink_green', ['#3AB370',"#EAE7CC","#FD1593"], N=256)
#Plot some SVGs
for g in [SVG]:
    fig=plot_relative_exp(adata, g, "x5", "x4", use_raw=True,color=color_self)
    fig.set_aspect('equal', 'box')
    fig.figure.savefig("../tutorial/"+g+"_raw_relative.png", dpi=300)







#----------------------------------------------------------------------------------------------------
#--------------------------Considering Neighborhoods-------------------------------------------------
#----------------------------------------------------------------------------------------------------
t=2
adata_cluster=detect_adataclusters(cell_id=adata.obs.index.tolist(), x=adata.obs["x2"].tolist(), y=adata.obs["x3"].tolist(), pred=adata.obs["pred"].tolist(), target_cluster=t, window_size=3, res=0.2)
adata.obs["adata"]=adata_cluster["adata_cluster_2"]
fig=sc.pl.scatter(adata,alpha=1,x="x5",y="x4",color="adata",show=False,size=200000/adata.shape[0])
#fig.axis('equal')
fig.set_aspect('equal', 'box')
fig.figure.savefig("./figures/adata2.png")



"""
#Spark relative expr
adata <- read.table("./raw_data/Rep11_MOB_count_matrix-1.tsv", check.names = F)
rn <- rownames(adata)
info <- cbind.data.frame(x = as.numeric(sapply(strsplit(rn, split = "x"), "[", 1)), 
                         y = as.numeric(sapply(strsplit(rn, split = "x"), "[", 2)))
rownames(info) <- rn

gene_plot <- c("Reln", "Cldn5", "Gng4", "Doc2g", "Kctd12", "Penk")
vst_ct <- var_stabilize(t(adata)) # R function in funcs.R
sig_vst_ct <- vst_ct[gene_plot, ]
rel_vst_ct <- apply(sig_vst_ct, 1, relative_func)


relative_func <- function(expres) {
    maxd = max(expres) - min(expres)
    rexpr = (expres - min(expres))/maxd
    return(rexpr)
}# end func

"""



#DE genes
{0: ['NEFM', 'NEFH', 'HAPLN4', 'GPX3', 'EPDR1', 'VAMP1'], 1: ['PCP4'], 2: [], 3: ['HPCAL1'], 4: ['PLP1', 'CNP', 'GFAP', 'PTGDS', 'CRYAB', 'S100B', 'MOBP', 'TF', 'SCD', 'SEPT4', 'CLDND1', 'GPM6B', 'MARCKSL1', 'APOD', 'PPP1R14A', 'SPP1', 'MAG', 'PLEKHB1', 'C4ORF48', 'QKI', 'RNASE1', 'QDPR', 'FEZ1', 'GPRC5B', 'PAQR6', 'CSRP1', 'TSC22D4', 'SELENOP', 'GSN', 'CLDN11', 'ABCA2', 'ERMN', 'BCAS1', 'LINC00844', 'SIRT2', 'EDIL3', 'H3F3B', 'NENF', 'CNTN2', 'OLIG1', 'SEPT8', 'MOG', 'ENPP2', 'CARNS1', 'RHOB', 'NCAM1', 'PRDX1', 'FIS1', 'MAL', 'SLC44A1', 'LAMP2', 'PMP22', 'LHPP', 'TMEM144', 'HTRA1', 'HSPA2', 'BIN1', 'PIP4K2A', 'CBR1', 'SUN2', 'EFHD1', 'C1ORF122', 'PLLP', 'PLA2G16', 'OPALIN', 'PTP4A2', 'NDRG1', 'FGFR2', 'CD9', 'RAPGEF5', 'ELOVL1', 'NKX6-2', 'TP53INP2', 'CERCAM', 'SLAIN1', 'PADI2', 'HIPK2', 'TTYH2', 'SLC48A1', 'AMER2', 'SHTN1', 'CNDP1', 'MYRF', 'KLK6', 'RNF130', 'HEPACAM', 'PLEKHH1', 'AATK', 'RASSF2', 'DDR1', 'PREX1', 'SEMA4D', 'FMNL2', 'UGT8', 'HAPLN2', 'SEMA3B', 'SLC12A2', 'TPPP3', 'GPR37', 'SPOCK3', 'ANKRD40', 'CFL2', 'C1ORF198', 'FGF1', 'ZEB2', 'RGCC', 'FAM107B', 'SOX8', 'CERS2', 'VWA1', 'CDKN1C', 'JAM3', 'ANLN', 'DAAM2', 'BOK', 'CDC42EP1', 'FA2H', 'CMTM5', 'PHLDB1', 'EVI2A', 'GLTP', 'GLDN', 'TMEM63A'], 5: [], 6: []}



