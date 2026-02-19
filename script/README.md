# Computational-Genomics-SpaGCN

## Table of Contents

- [Computational-Genomics-SpaGCN](#computational-genomics-spagcn)
  - [Short intro about SpaGCN](#short-intro-about-spagcn)
  - [SpaGCN presentation](#spagcn-presentation)
  - [SpaGCN script](#spagcn-script)
    - [Introduction](#introduction)
    - [Commands](#commands)
      - [--mode](#--mode)
        - [--histology_image](#--histology_image)
        - [--pixels](#--pixels)
        - [--arrays](#--arrays)
        - [--start](#--start)
        - [--raws](#--raws)
      - [--verbose](#--verbose)
      - [--read_obs](#--read_obs)
      - [--read_specific_obs](#--read_specific_obs)
  - [Dockerfile](#dockerfile)
  - [File-E9.5_E1S1.MOSTA.h5ad](#file-e95_e1s1mostah5ad)
  - [File-Dorsal_midbrain_cell_bin.h5ad](#file-dorsal_midbrain_cell_binh5ad)
  - [Presentation-and-youtube](#presentation-and-youtube)

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
```
usage: Computational Genomics Project SpaGCN [-h] [--mode MODE [MODE ...]] [--histology_image HISTOLOGY_IMAGE] [--start START] [--pixels PIXELS PIXELS] [--arrays ARRAYS ARRAYS] [--raws RAWS RAWS RAWS RAWS]
                                             [--read_obs READ_OBS] [--read_specific_obs READ_SPECIFIC_OBS READ_SPECIFIC_OBS] [--verbose {DEBUG,INFO,WARNING,ERROR,CRITICAL}]

This program does multiple things, it can convert h5 file to h5ad. It does Integrate gene expression and histology into a Graph, Spatial domain detection using SpaGCN, Identify SVGs, Identify Meta Gene and Multiple tissue
sections analysis.

options:
  -h, --help            show this help message and exit
  --mode MODE [MODE ...]
                        Modes on which these programs work.(Mode 1: integrate gene & spatial domains(default), Mode 2: identify svg & identify meta, Mode 3: multiple tissue)
  --histology_image HISTOLOGY_IMAGE
                        Integrate and histology into a Graph. The path to the .tif file is required. It is used in combination with command: mode 1 or mode without number(default).
  --start START         Start value for search l. It is used in combination with command: mode 1 or mode.
  --pixels PIXELS PIXELS
                        The x and y coordinates for pixels are typed here. It is used in combination with commands: mode, mode 1 and mode 2
  --arrays ARRAYS ARRAYS
                        The x and y coordinates for arrays are typed here. It is used in combination with commands: mode, mode 1
  --raws RAWS RAWS RAWS RAWS
                        The x array, y array, x pixel and y pixel coordinates for rows are typed here. It is used in combination with commands: mode 2
  --read_obs READ_OBS   Print all keys found in your file.
  --read_specific_obs READ_SPECIFIC_OBS READ_SPECIFIC_OBS
                        Print all specific keys found in your file.
  --verbose {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Set the logging level
```
  
------------

##### --mode

This command represents the modes of our project, we have 3 modes:

------------

##### mode 1(or mode without number) spatial domains and(this mode is also the default, if no mode number is defined)

This command accepts one parameter, parameter must be a file of type .h5ad.
`python3 spaGCN.py -mode 1 ./path/to/file/h5ad`

Example:
`python3 spaGCN.py --mode 1 ./data/151673/sample_data.h5ad`

As a response we get the path to the csv file.
> INFO:root:Done, your csv file is located here -> ./data/sample_data.csv

###### --histology_image
This command is only used in combination with -mode 1 or mode and accepts one parameter, parameter must be a file of type .tif or .png or .jpg.
`python3 spaGCN.py --mode 1 ./path/to/file/h5ad --histology_image ./path/to/histology/file`
Example
`python3 spaGCN.py --mode 1 ./data/151673/sample_data.h5ad --histology_image ./data/151673/histology.tif`
As a response we get the path to the csv file.
> INFO:root:Done, your csv file is located here -> ./data/sample_data.csv

If we have a histology, we get the result and a picture like in this location -> [PHOTO](https://github.com/Master-Computational-Genomics-SpaGCN/spaGCN-presentation/blob/main/resultImg/sample_data.jpg "PHOTO")

###### --pixels
This command is only used in combination with commands: mode, mode 1, mode 2. And they accept two parameters representing the coordinates. If this command is not used, then the program randomly selects the coordinates and informs the user.
`python3 spaGCN.py --mode 1/or2 ./path/to/file/h5ad --histology_image ./path/to/histology/file --pixels x y`
Example:
`python3 spaGCN.py --mode 1 ./data/151673/sample_data.h5ad --histology_image ./data/151673/histology.tif --pixels x4 x5`

###### --arrays
This command is only used in combination with commands: mode, mode 1, mode 2. And they accept two parameters representing the coordinates. If this command is not used, then the program randomly selects the coordinates and informs the user (you can see it in the first two pictures).
`python3 spaGCN.py -mode 1 ./path/to/file/h5ad --histology_image ./path/to/histology/file` +  `--pixels x y` + `--arrays x y`

###### --start
An int type value is added here, where we determine the starting point for determining the value of l. The default value is 0.01.
`python3 spaGCN.py -mode 1/or2 ./path/to/file/h5ad --histology_image ./path/to/histology/file` +  `--pixels x y` + `--arrays x y` + `--start 0.01`

As a result, the result.h5ad file and images of Spatial Domains and Refined Spatial Domains are obtained.
> INFO:root:Done, your result file and pictures are located here -> ./sample_results/sample_data_results.h5ad ./sample_results/sample_data_pred.png ./sample_results/sample_data_refined_pred.png

![mode_1](https://media.discordapp.net/attachments/962669952362496003/1113437577433853962/image.png?width=1173&height=662 "mode_1")

Results:
![mode_1_res](https://media.discordapp.net/attachments/962669952362496003/1113438021774221352/image.png?width=1138&height=662 "mode_1_res")

------------

##### mode 2 -  identify svg and identify meta are executed
This command identifies svg based on results from --mode 1. And accepts two parameters, the first parameter must be a original file of type .h5ad and the second parameter must be a file of type .h5ad(The file we got as a result of executing the --mode 1 command).
`python3 spaGCN.py -mode 2 ./path/to/file/h5ad ./path/to/file/h5ad` +  `--pixels x y` + `--arrays x y`

###### --raws
The x array, y array, x pixel and y pixel coordinates for rows are typed here for result.h5ad. It is used in combination with commands: --mode 2.
`python3 spaGCN.py --mode 2 ./path/to/file/h5ad ./path/to/file/h5ad` +  `--pixels x y` + `--arrays x y` + `--raws x1 y1 x2 y2`
Example:
`python3 spaGCN.py --mode 2 ./data/151673/sample_data.h5ad ./sample_results/sample_data_results.h5ad --pixels x y --arrays x x --raws x2 x3 x4 x5`

As a result we get svg images:
> INFO:root:Done, your pictures are located here -> ./sample_results/TMSB10.png ./sample_results/PCP4.png 

![mode_2](https://media.discordapp.net/attachments/962669952362496003/1113439282120622141/image.png?width=1280&height=662 "mode_2")

Results:
![identify_csv_img](https://media.discordapp.net/attachments/962669952362496003/1106648702174634044/image.png?width=1005&height=535 "identify_csv_img")

> Observed differences in obtaining results:
Based on sample_data.h5ad and sample_data_results.h5ad only managed to find two svg('TMSB10', 'PCP4'). While 5 svg('CAMK2N1', 'ENC1', 'GPM6A', 'ARPP19', 'HPCAL1') were found in the tutorial description. Among the files for downloading the tutorial there is also a differently generated result.h5ad through which I get the same results as in the tutorial, so I come to the conclusion that some part of the algorithm has been changed in the meantime because it generates a different result.h5ad that leads to a different result. Below is a proof image that I get the same result as in the tutorial using their result.h5ad file.
![identify_csv_problem](https://media.discordapp.net/attachments/962669952362496003/1113441263233007629/image.png?width=1300&height=662 "identify_csv_problem")
![identify_csv_problem_img](https://media.discordapp.net/attachments/962669952362496003/1106651242194808944/image.png?width=883&height=535 "identify_csv_problem_img")

For identify meta ->
As a result we get images:
> INFO:root:Done, your pictures are located here -> ./sample_results/identify_meta_GFAP.png ./sample_results/meta_finally_gene.png

Results:
![identify_meta_img](https://media.discordapp.net/attachments/962669952362496003/1106653290743210014/image.png?width=1025&height=535 "identify_meta_img")

> What do the results look like due to the problem mentioned above due to a different generated result. I get the same results here as in the tutorial.
![identify_meta_problem_img](https://media.discordapp.net/attachments/962669952362496003/1106653979984809984/image.png?width=1016&height=535 "identify_meta_problem_img")

------------

##### mode 3 - multiple tissue are executed
This command analyzes multiple tissue sections. It receives 4 parameters, the first parameter is the h5ad of the first tissue, the second parameter is the h5ad of the second tissue, the third parameter is the tif of the first tissue and the fourth parameter is the tif of the second tissue. This command can also accept more than 2 tissues.
`python3 spaGCN.py -mode 3 ./paths/to/files/h5ad ./paths/to/files/tif...`

Example:
`python3 spaGCN.py --multiple_tissue ./data/Mouse_brain/MA1.h5ad ./data/Mouse_brain/MP1.h5ad ./data/Mouse_brain/MA1_histology.tif ./data/Mouse_brain/MP1_histology.tif`

As a result we get pictures:
> INFO:root:Done, your picture are located here -> ./sample_results/muti_sections_domains_MA1.png

![mode_3](https://media.discordapp.net/attachments/962669952362496003/1113441979217485924/image.png?width=1440&height=318 "mode_3")

Results:

![multiple_tissue_img](https://media.discordapp.net/attachments/962669952362496003/1106659149338661025/image.png?width=911&height=533 "multiple_tissue_img")

------------

##### --verbose

The --verbose command is to provide a command-line option for setting the logging level of the application, allowing the user to control the verbosity of log messages. The --verbose argument is expected to have one of the following choices: 'DEBUG', 'INFO', 'WARNING', 'ERROR', or 'CRITICAL'. It is used to set the logging level of the application. If --verbose is not provided, the default logging level is used. The log messages are printed to the standard output (sys.stdout).

------------

##### --read_obs
If you don't know the keys that are in your h5ad file, you can see them using this command. It accepts one parameter, which is a file of type h5ad. That way you can choose coordinates for pixels, arrays and raws.


------------

##### --read_specific_obs
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
![docker](https://media.discordapp.net/attachments/962669952362496003/1113449620564082858/image.png?width=1185&height=662 "docker")

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

------------

### Presentation-and-youtube

Link to final presentation(in Serbian): [FINAL PRESENTATION](https://github.com/Master-Computational-Genomics-SpaGCN/spaGCN-presentation/blob/main/finalPresentation/SpaGCN_rezultati.pptx "FINAL PRESENTATION")

Link to youtube video: [YOUTUBE LINK](https://youtu.be/e4Ofv_IijKw "YOUTUBE LINK")
