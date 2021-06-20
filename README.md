# SpaGCN

## SpaGCN: Integrating gene expression, spatial location and histology to identify spatial domains and spatially variable genes by graph convolutional network


### Jian Hu*, Xiangjie Li, Kyle Coleman, Amelia Schroeder, Nan Ma, David J. Irwin, Edward B. Lee, Russell T. Shinohara, Mingyao Li*

SpaGCN is a graph convolutional network to integrate gene expression and histology to identify spatial domains and spatially variable genes. To jointly model all spots in a tissue slide, SpaGCN integrates information from gene expression, spatial locations and histological pixel intensities across spots into an undirected weighted graph. Each vertex in the graph contains gene expression information of a spot and the edge weight between two vertices quantifies their expression similarity that is driven by spatial dependency of their coordinates and the corresponding histology. To aggregate gene expression of each spot from its neighboring spots, SpaGCN utilizes a convolutional layer based on edge weights specified by the graph. The aggregated gene expression is then fed into a deep embedding clustering algorithm to cluster the spots into different spatial domains. After spatial domains are identified, genes that are enriched in each spatial domain can be detected by differential expression analysis between domains. SpaGCN is applicable to both in-situ transcriptomics with single-cell resolution (seqFISH, seqFISH+, MERFISH, STARmap, and FISSEQ) and spatial barcoding based transcriptomics (Spatial Transcriptomics , SLIDE-seq, SLIDE-seqV2, HDST, 10x Visium, DBiT-seq, Stero-seq, and PIXEL-seq) data. 

![SpaGCN workflow](docs/asserts/images/workflow.jpg)
<br>
For thorough details, see the preprint: [Biorxiv](https://www.biorxiv.org/content/10.1101/2020.11.30.405118v1)
<br>

## Usage

The [**SpaGCN**](https://github.com/jianhuupenn/SpaGCN) package is an implementation of a garph convolutional network for spatial transcriptomics. With SpaGCN, you can:

- Preprocess spatial transcriptomics data from various formats.
- Build a graph convolutional network with deep iterative clustering algorithm to identify spatial domains
- identify spatially variable genes for each spatial domain.
- Create mete genes to mark each spatial domains

<br>
For tutorial, please refer to: https://github.com/jianhuupenn/SpaGCN/blob/master/tutorial/tutorial.md
<br>
A Jupyter notebook of the tutorial is accessible from : https://github.com/jianhuupenn/SpaGCN/blob/master/tutorial/tutorial.ipynb
<br>
Please install jupyter in order to open this notebook.
<br>
Toy data and results can be downloaded at: https://drive.google.com/drive/folders/1zten54vkjorp26T4iD0ApQGa9ut5eY42?usp=sharing

## System Requirements
Python support packages: igraph, torch, pandas, numpy, scipy, scanpy > 1.5, anndata, louvain, sklearn.

## Versions the software has been tested on
Environment 1:
- System: Mac OS 10.13.6
- Python: 3.7.0
- Python packages: pandas = 1.1.3, numpy = 1.18.1, python-igraph=0.7.1,torch=1.5.1,louvain=0.6.1,scipy = 1.4.1, scanpy = 1.5.1, anndata = 0.6.22.post1, natsort = 7.0.1, sklearn = 0.22.1

Environment 2:
- System: Anaconda
- Python: 3.7.9
- Python packages: pandas = 1.1.3, numpy = 1.20.2, python-igraph=0.8.3, torch=1.6.0,louvain=0.7.0, scipy = 1.5.2, scanpy = 1.6.0, anndata = 0.7.4, natsort = 7.0.1, sklearn = 0.23.3

Environment 3:
- System: Anaconda
- Python: 3.8.8
- Python packages: pandas = 1.2.4, numpy = 1.19.1, python-igraph=0.9.1, torch=1.8.1, louvain=0.7.0, scipy = 1.6.3, scanpy = 1.7.2, anndata = 0.7.6, natsort = 7.1.1, sklearn = 0.24.2


## Contributing

Souce code: [Github](https://github.com/jianhuupenn/SpaGCN)  

We are continuing adding new features. Bug reports or feature requests are welcome. 

Last update: 06/19/2021, version 1.2.0



## References

Please consider citing the following reference:

- https://www.biorxiv.org/content/10.1101/2020.11.30.405118v1

<br>
