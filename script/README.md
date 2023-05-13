# Computational-Genomics-SpaGCN

## Table of Contents

- [Computational-Genomics-SpaGCN](#computational-genomics-spagcn)
  - [Short intro about SpaGCN](#short-intro-about-spagcn)
  - [SpaGCN presentation](#spagcn-presentation)
  - [SpaGCN script](#spagcn-script)
    - [Introduction](#introduction)
    - [Commands](#commands)
      - [-convert_h5](#-convert_h5)
      - [-integrate_gene](#-integrate_gene)
        - [-histology](#-histology)
        - [-pixels](#-pixels)
      - [-spatial_domains](#-spatial_domains)
        - [-arrays](#-arrays)
        - [-clusters](#-clusters)
        - [-start](#-start)
      - [-identify_csv](#-identify_csv)
        - [-raws](#-raws)
      - [-identify_meta](#-identify_meta)
      - [-multiple_tissue](#-multiple_tissue)
      - [-read_keys](#-read_keys)
      - [-read_specific_keys](#-read_specific_keys)
  - [Dockerfile](#dockerfile)
  - [File-E9.5_E1S1.MOSTA.h5ad](#file-e95_e1s1mostah5ad)
  - [File-Dorsal_midbrain_cell_bin.h5ad](#file-dorsal_midbrain_cell_binh5ad)

### Short intro about SpaGCN

SpaGCN is a graph convolutional network to integrate gene expression and histology to identify spatial domains and spatially variable genes. To jointly model all spots in a tissue slide, SpaGCN integrates information from gene expression, spatial locations, and histological pixel intensities across spots into an undirected weighted graph.

------------


### SpaGCN presentation
A short presentation was made for this project following the [SpaGCN paper](https://www.biorxiv.org/content/10.1101/2020.11.30.405118v1.full "SpaGCN paper").

The presentation is located in the following repository -> [PRESENTATION LINK](https://github.com/Master-Computational-Genomics-SpaGCN/spaGCN-presentation "PRESENTATION LINK")

----
### SpaGCN script
#### Introduction

A script has been created that performs the functionality shown in the [tutorial](https://github.com/jianhuupenn/SpaGCN/tree/master/tutorial "tutorial") using the [argparse](https://docs.python.org/3/library/argparse.html "argparse") library.

The script is run as follows:
`python3 spaGCN.py [command]`

You need python, SpaGCN and cv2 to even run the command.
You can install python from this link [python](https://www.python.org/downloads/ "python"). Python version 3.10.6 was used.

To install SpaGCN and cv2 you need pip.
`sudo apt-get -y install python3-pip`

After that you can install the necessary libraries:
`pip3 install SpaGCN`
`pip3 install opencv-python`

> Using Python on Windows can be problematic, to avoid wasting time install [Ubuntu 22.04.2 LTS](https://apps.microsoft.com/store/detail/ubuntu-22042-lts/9PN20MSR04DW?hl=en-rs&gl=rs&rtc=1 "Ubuntu 22.04.2 LTS") from the Windows Store. After that you install the [WSL](https://code.visualstudio.com/docs/remote/wsl "WSL") extension on VSCode and connect to ubuntu where everything will work flawlessly for you.

------------


#### Commands

By typing the following command, we get a list of commands that we can use
`python3 spaGCN.py -h`

> Computational Genomics Project SpaGCN [-h] [-convert_h5 CONVERT_H5 CONVERT_H5] [-integrate_gene INTEGRATE_GENE] [-histology HISTOLOGY]
                                             [-spatial_domains SPATIAL_DOMAINS SPATIAL_DOMAINS] [-clusters CLUSTERS] [-start START]
                                             [-identify_csv IDENTIFY_CSV IDENTIFY_CSV] [-identify_meta IDENTIFY_META IDENTIFY_META]
                                             [-multiple_tissue MULTIPLE_TISSUE MULTIPLE_TISSUE MULTIPLE_TISSUE MULTIPLE_TISSUE] [-pixels PIXELS PIXELS]
                                             [-arrays ARRAYS ARRAYS] [-raws RAWS RAWS RAWS RAWS] [-read_keys READ_KEYS]
                                             [-read_specific_keys READ_SPECIFIC_KEYS READ_SPECIFIC_KEYS] [--version]

> This program does multiple things, it can convert h5 file to h5ad. It does Integrate gene expression and histology into a Graph, Spatial domain detection using
SpaGCN, Identify SVGs, Identify Meta Gene and Multiple tissue sections analysis.

> options:
  -h, --help            show this help message and exit
  -convert_h5 CONVERT_H5 CONVERT_H5
                        Read original 10x_h5 data and save it to h5ad. The path to the h5 file is required.
  -integrate_gene INTEGRATE_GENE
                        Integrate gene expression into a Graph. The path to the h5ad file is required
  -histology HISTOLOGY  Integrate and histology into a Graph. The path to the .tif file is required. It is used in combination with command: integrate_gene.
  -spatial_domains SPATIAL_DOMAINS SPATIAL_DOMAINS
                        Spatial domain detection using SpaGCN. Paths to h5ad and csv files are required.
  -clusters CLUSTERS    Number of clusters for spatial domains. It is used in combination with command: spatial_domains.
  -start START          Start value for search l. It is used in combination with command: spatial_domains.
  -identify_csv IDENTIFY_CSV IDENTIFY_CSV
                        Identify SVGs. Paths to h5ad gene matrix and h5ad results files are required.
  -identify_meta IDENTIFY_META IDENTIFY_META
                        Identify Meta Gene. Paths to h5ad gene matrix and h5ad results files are required.
  -multiple_tissue MULTIPLE_TISSUE MULTIPLE_TISSUE MULTIPLE_TISSUE MULTIPLE_TISSUE
                        Multiple tissue sections analysis.Paths to h5ad first tissue, h5ad second tissue, tif first tissue and tif second tissue files are
                        required.
  -pixels PIXELS PIXELS
                        The x and y coordinates for pixels are typed here. It is used in combination with commands: integrate_gene, spatial_domains, identify_meta.
  -arrays ARRAYS ARRAYS
                        The x and y coordinates for arrays are typed here. It is used in combination with commands: spatial_domains, identify_csv, identify_meta.
  -raws RAWS RAWS RAWS RAWS
                        The x array, y array, x pixel and y pixel coordinates for rows are typed here. It is used in combination with commands: identify_csv,
                        identify_meta.
  -read_keys READ_KEYS  Print all keys found in your file.
  -read_specific_keys READ_SPECIFIC_KEYS READ_SPECIFIC_KEYS
                        Print all specific keys found in your file.
  --version             show program's version number and exit

------------


##### -convert_h5
This command reads original 10x_h5 data and save it to h5ad.
This command accepts two parameters, the first parameter must be a file of type .h5 and the second parameter must be a file of type .txt with positions.
`python3 spaGCN.py -convert_h5 ./path/to/file/h5 ./path/to/file/txt`
Example:
`python3 spaGCN.py -convert_h5 ./data/151673/expression_matrix.h5 ./data/151673/positions.txt`
As a response we get the path to the h5ad file.
> Done, your converted h5 file are located here -> ./data/expression_matrix.h5ad

![-convert_h5](https://media.discordapp.net/attachments/962669952362496003/1106639329717125240/image.png?width=1101&height=60 "-convert_h5")

------------


##### -integrate_gene
This command  integrate gene expression into a Graph. 
This command accepts one parameter, parameter must be a file of type .h5ad.
`python3 spaGCN.py -integrate_gene ./path/to/file/h5ad`
Example:
`python3 spaGCN.py -integrate_gene ./data/151673/sample_data.h5ad`
As a response we get the path to the csv file.
> Done, your csv file is located here -> ./data/sample_data.csv

![-integrate_gene](https://media.discordapp.net/attachments/962669952362496003/1106640045449957498/image.png?width=1101&height=85 "-integrate_gene")

###### -histology
This command is only used in combination with -integrate_gene and accepts one parameter, parameter must be a file of type .tif or .png or .jpg.
`python3 spaGCN.py -integrate_gene ./path/to/file/h5ad -histology ./path/to/histology/file`
Example
`python3 spaGCN.py -integrate_gene ./data/151673/sample_data.h5ad -histology ./data/151673/histology.tif`
As a response we get the path to the csv file.
> Done, your csv file is located here -> ./data/sample_data.csv

![histology](https://media.discordapp.net/attachments/962669952362496003/1106641138695291071/image.png?width=1101&height=95 "histology")

###### -pixels
This command is only used in combination with commands: integrate_gene, spatial_domains, identify_meta. And they accept two parameters representing the coordinates. If this command is not used, then the program randomly selects the coordinates and informs the user (you can see it in the previous two pictures).
`python3 spaGCN.py -integrate_gene ./path/to/file/h5ad -histology ./path/to/histology/file -pixels x y`
Example:
`python3 spaGCN.py -integrate_gene ./data/151673/sample_data.h5ad -histology ./data/151673/histology.tif -pixels x4 x5`

As a response we get the path to the csv file.
> Done, your csv file is located here -> ./data/sample_data.csv

![pixels](https://media.discordapp.net/attachments/962669952362496003/1106642719285198909/image.png?width=1101&height=75 "pixels")

------------

##### -spatial_domains
This command detects the spatial domain using SpaGCN. And accepts two parameters, the first parameter must be a file of type .h5ad and the second parameter must be a file of type .csv(The file we got as a result of executing the -integrate_gene command).
`python3 spaGCN.py -spatial_domains ./path/to/file/h5ad ./path/to/file/csv` + 
`-pixels x y`
###### -arrays
This command is only used in combination with commands: spatial_domains, identify_csv, identify_meta. And they accept two parameters representing the coordinates. If this command is not used, then the program randomly selects the coordinates and informs the user (you can see it in the first two pictures).
`python3 spaGCN.py -spatial_domains ./path/to/file/h5ad ./path/to/file/csv` + 
`-pixels x y` + `-arrays x y`
###### -clusters
If you know how many clusters are needed for your file then this command is also used, it accepts one parameter of type int(number). The default value is 7.
`python3 spaGCN.py -spatial_domains ./path/to/file/h5ad ./path/to/file/csv` + 
`-pixels x y` + `-arrays x y` + `-clusters 10

###### -start
An int type value is added here, where we determine the starting point for determining the value of l. The default value is 0.01.
`python3 spaGCN.py -spatial_domains ./path/to/file/h5ad ./path/to/file/csv` + 
`-pixels x y` + `-arrays x y` + `-clusters 10` + `-start 0.01`

Example:
`python3 spaGCN.py -spatial_domains ./data/151673/sample_data.h5ad ./data/sample_data.csv -pixels x4 x5 -arrays x2 x3 -clusters 7`
As a result, the result.h5ad file and images of Spatial Domains and Refined Spatial Domains are obtained.
>Done, your result file and pictures are located here -> ./sample_results/sample_data_results.h5ad ./sample_results/sample_data_pred.png ./sample_results/sample_data_refined_pred.png

![spatial_domains](https://media.discordapp.net/attachments/962669952362496003/1106645677271949353/image.png?width=966&height=535 "spatial_domains")

Results:
![spatial_domains_img](https://media.discordapp.net/attachments/962669952362496003/1106646230433546250/image.png?width=926&height=535 "spatial_domains_img")

------------

##### -identify_csv
This command identifies csv based on results from -spatial_domains. And accepts two parameters, the first parameter must be a original file of type .h5ad and the second parameter must be a file of type .h5ad(The file we got as a result of executing the -spatial_domains command).
`python3 spaGCN.py -identify_csv ./path/to/file/h5ad ./path/to/file/h5ad` + `-arrays x y`

###### -raws
The x array, y array, x pixel and y pixel coordinates for rows are typed here for result.h5ad. It is used in combination with commands: identify_csv, identify_meta.
`python3 spaGCN.py -identify_csv ./path/to/file/h5ad ./path/to/file/h5ad` + `-arrays x y` + `-raws x1 y1 x2 y2`
Example:
`python3 spaGCN.py -identify_csv ./data/151673/sample_data.h5ad ./sample_results/sample_data_results.h5ad -arrays x2 x3 -raws x2 x3 x4 x5`

As a result we get svg images:
> Done, your pictures are located here -> ./sample_results/TMSB10.png ./sample_results/PCP4.png 

![identify_csv](https://media.discordapp.net/attachments/962669952362496003/1106648046932066344/image.png?width=1101&height=313 "identify_csv")

Results:
![identify_csv_img](https://media.discordapp.net/attachments/962669952362496003/1106648702174634044/image.png?width=1005&height=535 "identify_csv_img")

> Observed differences in obtaining results:
Based on sample_data.h5ad and sample_data_results.h5ad only managed to find two svg('TMSB10', 'PCP4'). While 5 svg('CAMK2N1', 'ENC1', 'GPM6A', 'ARPP19', 'HPCAL1') were found in the tutorial description. Among the files for downloading the tutorial there is also a differently generated result.h5ad through which I get the same results as in the tutorial, so I come to the conclusion that some part of the algorithm has been changed in the meantime because it generates a different result.h5ad that leads to a different result. Below is a proof image that I get the same result as in the tutorial using their result.h5ad file.
![identify_csv_problem](https://media.discordapp.net/attachments/962669952362496003/1106650686705369108/image.png?width=1101&height=330 "identify_csv_problem")
![identify_csv_problem_img](https://media.discordapp.net/attachments/962669952362496003/1106651242194808944/image.png?width=883&height=535 "identify_csv_problem_img")

------------

##### -identify_meta
It is completely identical to the previous command. And accepts two parameters, the first parameter must be a original file of type .h5ad and the second parameter must be a file of type .h5ad(The file we got as a result of executing the -spatial_domains command). There is only one difference, here the command pixels is added.
`python3 spaGCN.py -identify_csv ./path/to/file/h5ad ./path/to/file/h5ad` + 
`-arrays x y` + `-pixels x y` + `-raws x1 y1 x2 y2`

Example:
`python3 spaGCN.py -identify_meta ./data/151673/sample_data.h5ad ./sample_results/sample_data_results.h5ad -arrays x2 x3 -pixels x4 x5 -raws x2 x3 x4 x5`

As a result we get images:
> Done, your pictures are located here -> ./sample_results/identify_meta_GFAP.png ./sample_results/meta_finally_gene.png

![identify_meta](https://media.discordapp.net/attachments/962669952362496003/1106652916846178335/image.png?width=1101&height=342 "identify_meta")

Results:
![identify_meta_img](https://media.discordapp.net/attachments/962669952362496003/1106653290743210014/image.png?width=1025&height=535 "identify_meta_img")

> What do the results look like due to the problem mentioned above due to a different generated result. I get the same results here as in the tutorial.
![identify_meta_problem_img](https://media.discordapp.net/attachments/962669952362496003/1106653979984809984/image.png?width=1016&height=535 "identify_meta_problem_img")

------------

##### -multiple_tissue

This command analyzes multiple tissue sections. It receives 4 parameters, the first parameter is the h5ad of the first tissue, the second parameter is the h5ad of the second tissue, the third parameter is the tif of the first tissue and the fourth parameter is the tif of the second tissue.
`python3 spaGCN.py -multiple_tissue ./path/to/file/h5ad ./path/to/file/h5ad ./path/to/file/tif ./path/to/file/tif`

Example:
`python3 spaGCN.py -multiple_tissue ./data/Mouse_brain/MA1.h5ad ./data/Mouse_brain/MP1.h5ad ./data/Mouse_brain/MA1_histology.tif ./data/Mouse_brain/MP1_histology.tif`

As a result we get pictures:
> Done, your picture are located here -> ./sample_results/muti_sections_domains_MA1.png

![multiple_tissue](https://media.discordapp.net/attachments/962669952362496003/1106658992442327052/image.png?width=1101&height=285 "multiple_tissue")

Results:

![multiple_tissue_img](https://media.discordapp.net/attachments/962669952362496003/1106659149338661025/image.png?width=911&height=533 "multiple_tissue_img")

------------

##### -read_keys
If you don't know the keys that are in your h5ad file, you can see them using this command. It accepts one parameter, which is a file of type h5ad. That way you can choose coordinates for pixels, arrays and raws.


------------

##### -read_specific_keys
It does the same as the previous command except that here goes another parameter which is of type string where we try to narrow down the selection of keys to select.

------------

### Dockerfile

Dockerfiles with requirement.txt can be found in Docker folder: [Docker & requirements](https://github.com/Master-Computational-Genomics-SpaGCN/SpaGCN/tree/master/script/Docker "Docker & requirements") 

With this dockerfile you can clone this repository and use script on docker container.
You need [docker](https://www.docker.com/ "docker") to run the following commands to create a docker image for you and run that docker image.
`docker build -t spagcn .`
And
`docker run --rm -it spagcn:latest`
add h5ad file/files to your docker container
`docker cp -r /path/to/my-folder/my-file my-container:/app/my-folder`
![docker](https://media.discordapp.net/attachments/962669952362496003/1106664137032994968/image.png?width=853&height=535 "docker")

The link to the DockerHub image is here ->[DockerHub SpaGCN](https://hub.docker.com/r/aleksa1902/spagcn/tags "DockerHub SpaGCN")

------------

###  File-E9.5_E1S1.MOSTA.h5ad

> The problem is the lack of information that I have on the site of those files and tutorials, I don't know in detail what each of the parameters means and what the values should be, mostly spagcn is adapted to their files, all parameters are adapted to their spample_data(start, step, tol, lr...) and that's why I don't get the expected results, but I get these.

![e95_fist](https://media.discordapp.net/attachments/962669952362496003/1106683839948001290/image.png?width=881&height=296 "e95_fist")
![E95_second](https://media.discordapp.net/attachments/962669952362496003/1106684132379066509/image.png?width=751&height=428 "E95_second")
![E95_third](https://media.discordapp.net/attachments/962669952362496003/1106684540434530324/image.png?width=881&height=365 "E95_third")
![E95_results](https://media.discordapp.net/attachments/962669952362496003/1106684667576451223/image.png?width=770&height=428 "E95_results")

Second chance with different keys

![E95_secondChance](https://media.discordapp.net/attachments/962669952362496003/1106728156636917800/image.png?width=765&height=428 "E95_secondChance")
![E95_secondChance_results](https://media.discordapp.net/attachments/962669952362496003/1106728356856221796/image.png?width=768&height=428 "E95_secondChance_results")
![E95_secondChance_csv](https://media.discordapp.net/attachments/962669952362496003/1106729815144411226/image.png?width=881&height=289 "E95_secondChance_csv")
![E95_secondChance_csv_results](https://media.discordapp.net/attachments/962669952362496003/1106730026323431474/image.png?width=688&height=428 "E95_secondChance_csv_results")

### File-Dorsal_midbrain_cell_bin.h5ad

> The file is huge, which leads to a much longer execution time, the command integrate_gene succeeds after a few minutes but the spatial_domains break, the error is MemoryError.
