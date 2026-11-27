# SpaGCN

## SpaGCN: Integrating gene expression and histology to identify spatial domains and spatially variable genes using graph convolutional networks

### Jian Hu,*, Xiangjie Li, Kyle Coleman, Amelia Schroeder, Nan Ma, David J. Irwin, Edward B. Lee, Russell T. Shinohara, Mingyao Li,*

SpaGCN is a graph convolutional network to integrate gene expression and histology to identify spatial domains and spatially variable genes. To jointly model all spots in a tissue slide, SpaGCN integrates information from gene expression, spatial locations and histological pixel intensities across spots into an undirected weighted graph. Each vertex in the graph contains gene expression information of a spot and the edge weight between two vertices quantifies their expression similarity that is driven by spatial dependency of their coordinates and the corresponding histology. To aggregate gene expression of each spot from its neighboring spots, SpaGCN utilizes a convolutional layer based on edge weights specified by the graph. The aggregated gene expression is then fed into a deep embedding clustering algorithm to cluster the spots into different spatial domains. After spatial domains are identified, genes that are enriched in each spatial domain can be detected by differential expression analysis between domains. SpaGCN is applicable to both in-situ transcriptomics with single-cell resolution (seqFISH, seqFISH+, MERFISH, STARmap, and FISSEQ) and spatial barcoding based transcriptomics (Spatial Transcriptomics , SLIDE-seq, SLIDE-seqV2, HDST, 10x Visium, DBiT-seq, Stero-seq, and PIXEL-seq) data. 
For more info, please go to: 
https://github.com/jianhuupenn/SpaGCN

### CUDA training and p sweeps

Single-slice clustering defaults to CPU. With a CUDA-enabled PyTorch 2.x
installation, set `device="cuda"` on `SpaGCN.train`, `search_res`, or
`detect_spatial_domains_ez_mode`. PCA and Louvain initialization run on CPU.

To cluster one slice at several p values, pass the preprocessed AnnData and its
dense distance matrix to `sweep_p`. The rows and columns of `adj` must follow
the AnnData spot order.

```python
import SpaGCN as spg

results = spg.sweep_p(
    adata, adj, p_values=[0.1, 0.3, 0.5, 0.7, 0.9],
    n_clusters=7, device="cuda",
)
labels = results[0.5].labels
```

`sweep_p` reuses one PCA embedding across p values. It searches for the length
scale and resolution and fits a model for each p. Each result provides `l`,
`resolution`, `labels`, and `probabilities`.
