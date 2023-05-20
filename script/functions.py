import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
from scipy.sparse import issparse
import random, torch
import warnings
warnings.filterwarnings("ignore")
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import SpaGCN as spg
import cv2
from anndata import AnnData
from scanpy import read_10x_h5

def CheckFolder(name):
    if not os.path.exists(name):
        os.makedirs(name)

def ConvertH5File(matrix, positions):
    # Read the input 10x_h5 data and save it to h5ad
    adata = read_10x_h5(matrix)
    
    # Read the positions file as a pandas dataframe
    spatial=pd.read_csv(positions,sep=",",header=None,na_filter=False,index_col=0) 
    
    # Add the spatial information to the h5ad object
    adata.obs["x1"]=spatial[1]
    adata.obs["x2"]=spatial[2]
    adata.obs["x3"]=spatial[3]
    adata.obs["x4"]=spatial[4]
    adata.obs["x5"]=spatial[5]
    adata.obs["x_array"]=adata.obs["x2"]
    adata.obs["y_array"]=adata.obs["x3"]
    adata.obs["x_pixel"]=adata.obs["x4"]
    adata.obs["y_pixel"]=adata.obs["x5"]
    
    # Select captured samples
    adata=adata[adata.obs["x1"]==1]
    
    # Convert the gene names to uppercase
    adata.var_names=[i.upper() for i in list(adata.var_names)]
    
    # Add a "genename" column to the var dataframe
    adata.var["genename"]=adata.var.index.astype("str")
    
    # Get the filename from the input path
    split = matrix.split("/")
    name = split[len(split)-1].split(".")
    
    # Check if the "data" folder exists, if not create it
    CheckFolder("data")
    
    # Set the output path
    pathName = "./data/" + name[0] + ".h5ad"
    
    # Save the h5ad object to a file
    adata.write_h5ad(pathName)

    return pathName

def RandomKeys(keys_list):
    random_key_x = random.choice(keys_list)
    random_key_y = random.choice(keys_list)
    print("random_key_x -> " + random_key_x + ", random_key_y -> " + random_key_y)
    return random_key_x, random_key_y

def IntegrateIntoGraphHistology(gene, histology, xpixel = "null", ypixel = "null"):
    # Read the gene expression data from the gene file using Scanpy
    adata=sc.read(gene)
    
    # Read the histology image using OpenCV
    img=cv2.imread(histology)

    # Set the x_pixel and y_pixel values in the gene expression data based on the x4 and x5 values in the obs dataframe
    if xpixel == "null" and ypixel == "null":
        keys_list = list(adata.obs.keys())
        random_key_x, random_key_y = RandomKeys(keys_list)
        x_pixel=adata.obs[random_key_x].tolist()
        y_pixel=adata.obs[random_key_y].tolist()
    else:
        x_pixel=adata.obs[xpixel].tolist()
        y_pixel=adata.obs[ypixel].tolist()

    # Set a window around each coordinate on the image to test if the coordinates are accurate
    img_new=img.copy()
    for i in range(len(x_pixel)):
        x=x_pixel[i]
        y=y_pixel[i]
        img_new[int(x-20):int(x+20), int(y-20):int(y+20),:]=0

    # Save the image with the marked coordinates to a folder
    split = gene.split("/")
    name = split[len(split)-1].split(".")
    # Check if the "data" folder exists, if not create it
    CheckFolder("sample_results")
    pathName = "./sample_results/" + name[0] + ".jpg"
    cv2.imwrite(pathName, img_new)

    s=1
    b=49
    adj=spg.calculate_adj_matrix(x=x_pixel,y=y_pixel, x_pixel=x_pixel, y_pixel=y_pixel, image=img, beta=b, alpha=s, histology=True)
    
    # Save the adjacency matrix to a CSV file in a folder
    split = gene.split("/")
    name = split[len(split)-1].split(".")
    pathName = "./data/" + name[0] + ".csv"
    np.savetxt(pathName, adj, delimiter=',')

    return pathName

def ReadKeys(gene):
    adata=sc.read(gene)
    print(f"All keys > {str(adata.obs.keys())}")

def ReadSpecificKeys(gene, search_str):
    adata = sc.read(gene)
    keys_with_str = [key for key in adata.obs.keys() if search_str in key]
    print(f"Keys with '{search_str}' > {keys_with_str}")

def IntegrateIntoGraph(gene, xpixel = "null", ypixel = "null"):
    # Read the gene expression data from the gene file using Scanpy
    adata=sc.read(gene)

    # Extract x and y coordinates from the gene data
    if xpixel == "null" and ypixel == "null":
        keys_list = list(adata.obs.keys())
        random_key_x, random_key_y = RandomKeys(keys_list)
        x_pixel=adata.obs[random_key_x].tolist()
        y_pixel=adata.obs[random_key_y].tolist()
    else:
        x_pixel=adata.obs[xpixel].tolist()
        y_pixel=adata.obs[ypixel].tolist()

    # Calculate adjacency matrix using x and y coordinates
    adj=spg.calculate_adj_matrix(x=x_pixel,y=y_pixel, histology=False)
    
    # Save the adjacency matrix to the output file
    split = gene.split("/")
    name = split[len(split)-1].split(".")
    CheckFolder("data")
    pathName = "./data/" + name[0] + ".csv"
    np.savetxt(pathName, adj, delimiter=',')

    return pathName

def SpatialDomainsDetectionSpaGCN(gene, adjCsv, clusters=7, xpixel = "null", ypixel = "null", xarray = "null", yarray = "null", startL=0.01):
    # Read the gene expression data
    adata=sc.read(gene)

    # Add spatial information to the data
    if xpixel == "null" and ypixel == "null":
        keys_list = list(adata.obs.keys())
        random_key_x, random_key_y = RandomKeys(keys_list)
        adata.obs["x_pixel"]=adata.obs[random_key_x]
        adata.obs["y_pixel"]=adata.obs[random_key_y]
    else:
        adata.obs["x_pixel"]=adata.obs[xpixel]
        adata.obs["y_pixel"]=adata.obs[ypixel]

    x_pixel=adata.obs["x_pixel"].tolist()
    y_pixel=adata.obs["y_pixel"].tolist()

    if xarray == "null" and yarray == "null":
        keys_list = list(adata.obs.keys())
        random_key_x, random_key_y = RandomKeys(keys_list)
        adata.obs["x_array"]=adata.obs[random_key_x]
        adata.obs["y_array"]=adata.obs[random_key_y]
    else:
        adata.obs["x_array"]=adata.obs[xarray]
        adata.obs["y_array"]=adata.obs[yarray]

    x_array=adata.obs["x_array"].tolist()
    y_array=adata.obs["y_array"].tolist()

    # Load the adjacency matrix
    adj=np.loadtxt(adjCsv, delimiter=',')

    # Make sure gene names are unique
    adata.var_names_make_unique()

    # Perform gene filtering
    spg.prefilter_genes(adata,min_cells=3)

    # Perform filtering of special genes
    spg.prefilter_specialgenes(adata)
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)
    p=0.5
    # TODO: find better way for these parameters
    l=spg.search_l(p, adj, start=startL, end=1000, tol=0.01, max_run=100)

    #If the number of clusters known, we can use the spg.search_res() fnction to search for suitable resolution(optional)
    #For this toy data, we set the number of clusters=7 since this tissue has 7 layers
    n_clusters=clusters
    #Set seed
    r_seed=t_seed=n_seed=100

    # Search for optimal value of the resolution parameter
    # TODO: find better way for these parameters
    res=spg.search_res(adata, adj, l, n_clusters, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed) 

    # Train the model
    clf=spg.SpaGCN()
    clf.set_l(l)

    # Set random seeds for reproducibility
    random.seed(r_seed)
    torch.manual_seed(t_seed)
    np.random.seed(n_seed)

    # Train the model with the given parameters
    clf.train(adata,adj,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)

    # Make predictions and add to the data
    y_pred, prob=clf.predict()
    adata.obs["pred"]= y_pred
    adata.obs["pred"]=adata.obs["pred"].astype('category')
    #Do cluster refinement(optional)
    #shape="hexagon" for Visium data, "square" for ST data.
    adj_2d=spg.calculate_adj_matrix(x=x_array,y=y_array, histology=False)
    refined_pred=spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), dis=adj_2d, shape="hexagon")
    adata.obs["refined_pred"]=refined_pred
    adata.obs["refined_pred"]=adata.obs["refined_pred"].astype('category')

    split = gene.split("/")
    name = split[len(split)-1].split(".")
    CheckFolder("sample_results")
    pathNameResult = "./sample_results/" + name[0] + "_results.h5ad"
    #Save results
    adata.write_h5ad(pathNameResult)

    adata=sc.read(pathNameResult)
    
    #Set colors used
    plot_color=["#F56867","#FEB915","#C798EE","#59BE86","#7495D3","#D1D1D1","#6D1A9C","#15821E","#3A84E6","#997273","#787878","#DB4C6C","#9E7A7A","#554236","#AF5F3C","#93796C","#F9BD3F","#DAB370","#877F6C","#268785"]
    #Plot spatial domains
    domains="pred"
    num_celltype=len(adata.obs[domains].unique())
    adata.uns[f"{domains}_colors"] = list(plot_color[:num_celltype])
    ax=sc.pl.scatter(adata,alpha=1,x="y_pixel",y="x_pixel",color=domains,title=domains,color_map=plot_color,show=False,size=100000/adata.shape[0])
    ax.set_aspect('equal', 'box')
    ax.axes.invert_yaxis()
    pathNamePred = "./sample_results/" + name[0] + "_pred.png"
    plt.savefig(pathNamePred, dpi=600)
    plt.close()

    #Plot refined spatial domains
    domains="refined_pred"
    num_celltype=len(adata.obs[domains].unique())
    adata.uns[f"{domains}_colors"] = list(plot_color[:num_celltype])
    ax=sc.pl.scatter(adata,alpha=1,x="y_pixel",y="x_pixel",color=domains,title=domains,color_map=plot_color,show=False,size=100000/adata.shape[0])
    ax.set_aspect('equal', 'box')
    ax.axes.invert_yaxis()
    pathNameRefPred = "./sample_results/" + name[0] + "_refined_pred.png"
    plt.savefig(pathNameRefPred, dpi=600)
    plt.close()

    return pathNameResult + " " + pathNamePred + " " + pathNameRefPred

def IdentifySVG(gene, results, xarray = "null", yarray = "null", rawxpixel = "null", rawypixel = "null", rawxarray = "null", rawyarray = "null"):
    #Read in raw data
    adata=sc.read(results)

    # Assign x and y values to the adata object
    if xarray == "null" and yarray == "null":
        keys_list = list(adata.obs.keys())
        random_key_x, random_key_y = RandomKeys(keys_list)
        x_array=adata.obs[random_key_x].tolist()
        y_array=adata.obs[random_key_y].tolist()
    else:
        x_array=adata.obs[xarray].tolist()
        y_array=adata.obs[yarray].tolist()

    # Read in gene data and assign x, y, x_pixel, and y_pixel values to the gene data
    raw=sc.read(gene)
    raw.var_names_make_unique()
    raw.obs["pred"]=adata.obs["pred"].astype('category')

    if rawxpixel == "null" and rawypixel == "null":
        keys_list = list(raw.obs.keys())
        random_key_x, random_key_y = RandomKeys(keys_list)
        raw.obs["x_pixel"]=raw.obs[random_key_x]
        raw.obs["y_pixel"]=raw.obs[random_key_y]
    else:
        raw.obs["x_pixel"]=raw.obs[rawxpixel]
        raw.obs["y_pixel"]=raw.obs[rawypixel]

    if rawxarray == "null" and rawyarray == "null":
        keys_list = list(raw.obs.keys())
        random_key_x, random_key_y = RandomKeys(keys_list)
        raw.obs["x_array"]=raw.obs[random_key_x]
        raw.obs["y_array"]=raw.obs[random_key_y]
    else:
        raw.obs["x_array"]=raw.obs[rawxarray]
        raw.obs["y_array"]=raw.obs[rawyarray]
    
    #Convert sparse matrix to non-sparse
    raw.X=(raw.X.A if issparse(raw.X) else raw.X)
    raw.raw=raw
    sc.pp.log1p(raw)

    #Use domain 0 as an example
    target=0
    #Set filtering criterials
    min_in_group_fraction=0.8
    min_in_out_group_ratio=1
    min_fold_change=1.5
    #Search radius such that each spot in the target domain has approximately 10 neighbors on average
    adj_2d=spg.calculate_adj_matrix(x=x_array, y=y_array, histology=False)
    start, end= np.quantile(adj_2d[adj_2d!=0],q=0.001), np.quantile(adj_2d[adj_2d!=0],q=0.1)
    r=spg.search_radius(target_cluster=target, cell_id=adata.obs.index.tolist(), x=x_array, y=y_array, pred=adata.obs["pred"].tolist(), start=start, end=end, num_min=10, num_max=14,  max_run=100)
    #Detect neighboring domains
    nbr_domians=spg.find_neighbor_clusters(target_cluster=target,
                                    cell_id=raw.obs.index.tolist(), 
                                    x=raw.obs["x_array"].tolist(), 
                                    y=raw.obs["y_array"].tolist(), 
                                    pred=raw.obs["pred"].tolist(),
                                    radius=r,
                                    ratio=1/2)

    nbr_domians=nbr_domians[0:3]
    de_genes_info=spg.rank_genes_groups(input_adata=raw,
                                    target_cluster=target,
                                    nbr_list=nbr_domians, 
                                    label_col="pred", 
                                    adj_nbr=True, 
                                    log=True)
    #Filter genes
    de_genes_info=de_genes_info[(de_genes_info["pvals_adj"]<0.05)]
    filtered_info=de_genes_info
    filtered_info=filtered_info[(filtered_info["pvals_adj"]<0.05) &
                                (filtered_info["in_out_group_ratio"]>min_in_out_group_ratio) &
                                (filtered_info["in_group_fraction"]>min_in_group_fraction) &
                                (filtered_info["fold_change"]>min_fold_change)]
    filtered_info=filtered_info.sort_values(by="in_group_fraction", ascending=False)
    filtered_info["target_dmain"]=target
    filtered_info["neighbors"]=str(nbr_domians)
    print("SVGs for domain ", str(target),":", filtered_info["genes"].tolist())

    filtered_info

    pathNameList = ""
    CheckFolder("sample_results")

    #Plot refinedspatial domains
    color_self = clr.LinearSegmentedColormap.from_list('pink_green', ['#3AB370',"#EAE7CC","#FD1593"], N=256)
    for g in filtered_info["genes"].tolist():
        raw.obs["exp"]=raw.X[:,raw.var.index==g]
        ax=sc.pl.scatter(raw,alpha=1,x="y_pixel",y="x_pixel",color="exp",title=g,color_map=color_self,show=False,size=100000/raw.shape[0])
        ax.set_aspect('equal', 'box')
        ax.axes.invert_yaxis()
        pathName = "./sample_results/"+g+".png"
        pathNameList += pathName + " "
        plt.savefig(pathName, dpi=600)
        plt.close()
    
    return pathNameList
    
def IdentifyMetaGene(gene, results, xpixel = "null", ypixel = "null", xarray = "null", yarray = "null", rawxpixel = "null", rawypixel = "null", rawxarray = "null", rawyarray = "null"):
    #Read in raw data
    adata=sc.read(results)

    if xpixel == "null" and ypixel == "null":
        keys_list = list(adata.obs.keys())
        random_key_x, random_key_y = RandomKeys(keys_list)
        adata.obs["x_pixel"]=adata.obs[random_key_x]
        adata.obs["y_pixel"]=adata.obs[random_key_y]
    else:
        adata.obs["x_pixel"]=adata.obs[ypixel]
        adata.obs["y_pixel"]=adata.obs[xpixel]

    if xarray == "null" and yarray == "null":
        keys_list = list(adata.obs.keys())
        random_key_x, random_key_y = RandomKeys(keys_list)
        adata.obs["x_array"]=adata.obs[random_key_x]
        adata.obs["y_array"]=adata.obs[random_key_y]
    else:
        adata.obs["x_array"]=adata.obs[xarray]
        adata.obs["y_array"]=adata.obs[yarray]

    # Read in gene expression data
    raw=sc.read(gene)
    raw.var_names_make_unique()
    raw.obs["pred"]=adata.obs["pred"].astype('category')

    if rawxpixel == "null" and rawypixel == "null":
        keys_list = list(raw.obs.keys())
        random_key_x, random_key_y = RandomKeys(keys_list)
        raw.obs["x_pixel"]=raw.obs[random_key_x]
        raw.obs["y_pixel"]=raw.obs[random_key_y]
    else:
        raw.obs["x_pixel"]=raw.obs[rawxpixel]
        raw.obs["y_pixel"]=raw.obs[rawypixel]

    if rawxarray == "null" and rawyarray == "null":
        keys_list = list(raw.obs.keys())
        random_key_x, random_key_y = RandomKeys(keys_list)
        raw.obs["x_array"]=raw.obs[random_key_x]
        raw.obs["y_array"]=raw.obs[random_key_y]
    else:
        raw.obs["x_array"]=raw.obs[rawxarray]
        raw.obs["y_array"]=raw.obs[rawyarray]

    #Convert sparse matrix to non-sparse
    raw.X=(raw.X.A if issparse(raw.X) else raw.X)
    raw.raw=raw
    sc.pp.log1p(raw)

    #Use domain 2 as an example
    target=2

    # Find meta gene
    meta_name, meta_exp=spg.find_meta_gene(input_adata=raw,
                        pred=raw.obs["pred"].tolist(),
                        target_domain=target,
                        start_gene="GFAP",
                        mean_diff=0,
                        early_stop=True,
                        max_iter=3,
                        use_raw=False)

    # Add meta gene expression to raw data
    raw.obs["meta"]=meta_exp

    color_self = clr.LinearSegmentedColormap.from_list('pink_green', ['#3AB370',"#EAE7CC","#FD1593"], N=256)

    #Plot meta gene
    g="GFAP"
    raw.obs["exp"]=raw.X[:,raw.var.index==g]
    ax=sc.pl.scatter(raw,alpha=1,x="y_pixel",y="x_pixel",color="exp",title=g,color_map=color_self,show=False,size=100000/raw.shape[0])
    ax.set_aspect('equal', 'box')
    ax.axes.invert_yaxis()
    CheckFolder("sample_results")
    pathNameFirst = f"./sample_results/identify_meta_{g}.png"
    plt.savefig(pathNameFirst, dpi=600)
    plt.close()

    # Plot meta gene expression
    raw.obs["exp"]=raw.obs["meta"]
    ax=sc.pl.scatter(raw,alpha=1,x="y_pixel",y="x_pixel",color="exp",title=meta_name,color_map=color_self,show=False,size=100000/raw.shape[0])
    ax.set_aspect('equal', 'box')
    ax.axes.invert_yaxis()
    pathNameSecond = "./sample_results/meta_finally_gene.png"
    plt.savefig(pathNameSecond, dpi=600)
    plt.close()

    return pathNameFirst + " " + pathNameSecond

def MultipleTissue(firstTissue, secondTissue, firstHistology, secondHistology):
    # Load the first and second tissue datasets
    adata1=sc.read(firstTissue)
    adata2=sc.read(secondTissue)

    # Load the first and second histology images
    img1=cv2.imread(firstHistology)
    img2=cv2.imread(secondHistology)

    # Define beta and scale factors for extracting color and z values from histology images and extract color and z values for cells in the first tissue dataset
    b=49
    s=1
    x_pixel1=adata1.obs["x4"].tolist()
    y_pixel1=adata1.obs["x5"].tolist()
    adata1.obs["color"]=spg.extract_color(x_pixel=x_pixel1, y_pixel=y_pixel1, image=img1, beta=b)
    z_scale=np.max([np.std(x_pixel1), np.std(y_pixel1)])*s
    adata1.obs["z"]=(adata1.obs["color"]-np.mean(adata1.obs["color"]))/np.std(adata1.obs["color"])*z_scale

    # Extract color and z values for cells in the second tissue dataset
    x_pixel2=adata2.obs["x4"].tolist()
    y_pixel2=adata2.obs["x5"].tolist()
    adata2.obs["color"]=spg.extract_color(x_pixel=x_pixel2, y_pixel=y_pixel2, image=img2, beta=b)
    z_scale=np.max([np.std(x_pixel2), np.std(y_pixel2)])*s
    adata2.obs["z"]=(adata2.obs["color"]-np.mean(adata2.obs["color"]))/np.std(adata2.obs["color"])*z_scale
    
    # Delete the histology images to free up memory
    del img1, img2

    # Add x and y pixel coordinates to the first and second tissue datasets
    adata1.obs["x_pixel"]=x_pixel1
    adata1.obs["y_pixel"]=y_pixel1
    adata2.obs["x_pixel"]=x_pixel2-np.min(x_pixel2)+np.min(x_pixel1)
    adata2.obs["y_pixel"]=y_pixel2-np.min(y_pixel2)+np.max(y_pixel1)

    # Make variable names unique in both tissue datasets
    adata1.var_names_make_unique()
    adata2.var_names_make_unique()

    # Concatenate the two tissue datasets into a single AnnData object
    adata_all=AnnData.concatenate(adata1, adata2,join='inner',batch_key="dataset_batch",batch_categories=["0","1"])

    # Calculate pairwise distance between cells using their x, y, and z coordinates
    X=np.array([adata_all.obs["x_pixel"], adata_all.obs["y_pixel"], adata_all.obs["z"]]).T.astype(np.float32)
    adj=spg.pairwise_distance(X)

    # Normalize gene expression values by total counts and log-transform
    sc.pp.normalize_per_cell(adata_all, min_counts=0)
    sc.pp.log1p(adata_all)
    p=0.5 
    #Find the l value given p
    l=spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)

    res=1.0
    seed=100
    random.seed(seed)
    torch.manual_seed(seed)
    np.random.seed(seed)

    # train SpaGCN model with given parameters
    clf=spg.SpaGCN()
    clf.set_l(l)
    clf.train(adata_all,adj,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)
    y_pred, prob=clf.predict()
    adata_all.obs["pred"]= y_pred
    adata_all.obs["pred"]=adata_all.obs["pred"].astype('category')

    colors_use=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#bcbd22', '#17becf', '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#bec1d4', '#bb7784', '#0000ff', '#111010', '#FFFF00',   '#1f77b4', '#800080', '#959595', 
    '#7d87b9', '#bec1d4', '#d6bcc0', '#bb7784', '#8e063b', '#4a6fe3', '#8595e1', '#b5bbe3', '#e6afb9', '#e07b91', '#d33f6a', '#11c638', '#8dd593', '#c6dec7', '#ead3c6', '#f0b98d', '#ef9708', '#0fcfc0', '#9cded6', '#d5eae7', '#f3e1eb', '#f6c4e1', '#f79cd4']
    num_celltype=len(adata_all.obs["pred"].unique())
    adata_all.uns["pred_colors"]=list(colors_use[:num_celltype])
    
    # create scatter plot of predicted cell types
    ax=sc.pl.scatter(adata_all,alpha=1,x="y_pixel",y="x_pixel",color="pred",show=False,size=150000/adata_all.shape[0])
    ax.set_aspect('equal', 'box')
    ax.axes.invert_yaxis()
    
    # save scatter plot to file
    CheckFolder("sample_results")
    split = firstTissue.split("/")
    name = split[len(split)-1].split(".")
    pathName = "./sample_results/muti_sections_domains_" + name[0] + ".png"
    plt.savefig(pathName, dpi=600)
    plt.close()

    return pathName


def MultipleTissueMore(tissues, histologies): #send arrays
    adata = []
    img = []
    for i in tissues:
        # Load tissues datasets
        adata.append(sc.read(tissues[i]))
        # Load histology images
        img.append(cv2.imread(histologies[i]))

    # Define beta and scale factors for extracting color and z values from histology images
    b=49
    s=1
    x_pixels = []
    y_pixels = []

    # Extract color and z values for cells in the tissues dataset
    for i in adata:
        x_pixels.append(adata[i].obs["x4"].tolist())
        y_pixels.append(adata[i].obs["x5"].tolist())
        adata[i].obs["color"]=spg.extract_color(x_pixel=x_pixels[i], y_pixel=y_pixels[i], image=img[i], beta=b)
        z_scale=np.max([np.std(x_pixels[i]), np.std(y_pixels[i])])*s
        adata[i].obs["z"]=(adata[i].obs["color"]-np.mean(adata[i].obs["color"]))/np.std(adata[i].obs["color"])*z_scale

    # Delete the histology images to free up memory
    del img

    # Add x and y pixel coordinates to the first and second tissue datasets
    adata[0].obs["x_pixel"]=x_pixels[0]
    adata[0].obs["y_pixel"]=y_pixels[0]
    for i in range(1, len(adata)):
        adata[i].obs["x_pixel"]=x_pixels[i]-np.min(x_pixels[i])+np.min(x_pixels[i-1])
        adata[i].obs["y_pixel"]=y_pixels[i]-np.min(y_pixels[i])+np.max(y_pixels[i-1])

    # Make variable names unique in both tissue datasets
    for i in adata:
        adata[i].var_names_make_unique()

    # Concatenate the two tissue datasets into a single AnnData object
    for i in range(0, len(adata)-1, 2):
        adata_all=AnnData.concatenate(adata[i], adata[i+1],join='inner',batch_key="dataset_batch",batch_categories=["0","1"])

    
    # Calculate pairwise distance between cells using their x, y, and z coordinates
    X=np.array([adata_all.obs["x_pixel"], adata_all.obs["y_pixel"], adata_all.obs["z"]]).T.astype(np.float32)
    adj=spg.pairwise_distance(X)

    # Normalize gene expression values by total counts and log-transform
    sc.pp.normalize_per_cell(adata_all, min_counts=0)
    sc.pp.log1p(adata_all)
    p=0.5 
    #Find the l value given p
    l=spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)

    res=1.0
    seed=100
    random.seed(seed)
    torch.manual_seed(seed)
    np.random.seed(seed)

    # train SpaGCN model with given parameters
    clf=spg.SpaGCN()
    clf.set_l(l)
    clf.train(adata_all,adj,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)
    y_pred, prob=clf.predict()
    adata_all.obs["pred"]= y_pred
    adata_all.obs["pred"]=adata_all.obs["pred"].astype('category')

    colors_use=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#bcbd22', '#17becf', '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#bec1d4', '#bb7784', '#0000ff', '#111010', '#FFFF00',   '#1f77b4', '#800080', '#959595', 
    '#7d87b9', '#bec1d4', '#d6bcc0', '#bb7784', '#8e063b', '#4a6fe3', '#8595e1', '#b5bbe3', '#e6afb9', '#e07b91', '#d33f6a', '#11c638', '#8dd593', '#c6dec7', '#ead3c6', '#f0b98d', '#ef9708', '#0fcfc0', '#9cded6', '#d5eae7', '#f3e1eb', '#f6c4e1', '#f79cd4']
    num_celltype=len(adata_all.obs["pred"].unique())
    adata_all.uns["pred_colors"]=list(colors_use[:num_celltype])
    
    # create scatter plot of predicted cell types
    ax=sc.pl.scatter(adata_all,alpha=1,x="y_pixel",y="x_pixel",color="pred",show=False,size=150000/adata_all.shape[0])
    ax.set_aspect('equal', 'box')
    ax.axes.invert_yaxis()
    
    # save scatter plot to file
    CheckFolder("sample_results")
    split = tissues[0].split("/")
    name = split[len(split)-1].split(".")
    pathName = "./sample_results/muti_sections_domains_" + name[0] + ".png"
    plt.savefig(pathName, dpi=600)
    plt.close()

    return pathName
