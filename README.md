# SpaGCN

## SpaGCN: Integrating gene expression and histology to identify spatial domains and spatially variable genes using graph convolutional networks

### Jian Hu, Xiangjie Li, Kyle Coleman, Edward B. Lee, Russell T. Shinohara, Mingyao Li*

SpaGCN is a graph convolutional network to integrate gene expression and histology to identify spatial domains and spatially variable genes. To jointly model all spots in a tissue slide, SpaGCN integrates information from gene expression, spatial locations and histological pixel intensities across spots into an undirected weighted graph. Each vertex in the graph contains gene expression information of a spot and the edge weight between two vertices quantifies their expression similarity that is driven by spatial dependency of their coordinates and the corresponding histology. To aggregate gene expression of each spot from its neighboring spots, SpaGCN utilizes a convolutional layer based on edge weights specified by the graph. The aggregated gene expression is then fed into a deep embedding clustering algorithm to cluster the spots into different spatial domains. After spatial domains are identified, genes that are enriched in each spatial domain can be detected by differential expression analysis between domains. SpaGCN is applicable to both fluorescence in situ hybridization (FISH)-based (e.g. seqFISH, seqFISH+, or MERFISH) and sequencing-based (e.g. 10X Visium spatial gene expression) data. 
![SpaGCN workflow](docs/asserts/images/workflow.jpg)
Figure created by Nan Ma.
For thorough details, see the preprint: [Bioxiv]()
<br>

## Usage

The [**SpaGCN**](https://github.com/jianhuupenn/SpaGCN) package is an implementation of a garph convolutional network for spatial transcriptomics. With SpaGCN, you can:

- Preprocess spatial transcriptomics data from various formats.
- Build a graph convolutional network with deep iterative clustering algorithm to identify spatial domains
- identify spatially variable genes for each spatial domain.
- Create mete genes to mark each spatial domains

<br>
For tutorial, please refer to: https://github.com/jianhuupenn/SpaGCN/blob/master/tutorial/tutorial.md

## Contributing

Souce code: [Github](https://github.com/jianhuupenn/SpaGCN)  

We are continuing adding new features. Bug reports or feature requests are welcome.

<br>


## References

Please consider citing the following reference:

- 
<br>
