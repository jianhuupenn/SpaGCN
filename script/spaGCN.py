import argparse
import logging
import os, sys
from functions import *


def PrintError(parser):
    logging.error("Your files do not exist or are of the wrong type.")
    parser.print_help()

# create parser object with program name and description
parser = argparse.ArgumentParser(prog='Computational Genomics Project SpaGCN', description='This program does multiple things, it can convert h5 file to h5ad. It does Integrate gene expression and histology into a Graph, Spatial domain detection using SpaGCN, Identify SVGs, Identify Meta Gene and Multiple tissue sections analysis.')

# add arguments to parser
parser.add_argument("--convert_h5", nargs=2, help="Read original 10x_h5 data and save it to h5ad. The path to the h5 file is required.", default='')
parser.add_argument("--integrate_gene", help="Integrate gene expression into a Graph. The path to the h5ad file is required", default='')
parser.add_argument("--histology_image", help="Integrate and histology into a Graph. The path to the .tif file is required. It is used in combination with command: integrate_gene.", default='')
parser.add_argument("--spatial_domains", nargs=2, help="Spatial domain detection using SpaGCN. Paths to h5ad and csv files are required.", default='')
parser.add_argument("--start", help="Start value for search l. It is used in combination with command: spatial_domains.", type=float, default=0.01)
parser.add_argument("--identify_svg", nargs=2, help="Identify SVGs. Paths to h5ad gene matrix and h5ad results files are required.", default='')
parser.add_argument("--identify_meta", nargs=2, help="Identify Meta Gene. Paths to h5ad gene matrix and h5ad results files are required.", default='')
parser.add_argument("--multiple_tissue", nargs='+', help="Multiple tissue sections analysis.Paths to h5ad first tissue, h5ad second tissue, tif first tissue and tif second tissue files are required.", default='')
parser.add_argument("--pixels", nargs=2, help="The x and y coordinates for pixels are typed here. It is used in combination with commands: integrate_gene, spatial_domains, identify_meta.", default='')
parser.add_argument("--arrays", nargs=2, help="The x and y coordinates for arrays are typed here. It is used in combination with commands: spatial_domains, identify_svg, identify_meta.", default='')
parser.add_argument("--raws", nargs=4, help="The x array, y array, x pixel and y pixel coordinates for rows are typed here. It is used in combination with commands: identify_svg, identify_meta.", default='')
parser.add_argument("--read_keys", help="Print all keys found in your file.", default='')
parser.add_argument("--read_specific_keys", nargs=2, help="Print all specific keys found in your file.", default='')
parser.add_argument("--mode", nargs="+", help="Modes on which these programs work.(Mode 1: spatial domains & identify svg, Mode 2: identify meta, Mode 3: multiple tissue)", default='')
parser.add_argument("--verbose", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], help="Set the logging level")


if __name__ == '__main__':

    # parse command line arguments
    parsed_args = parser.parse_args()

    if parsed_args.verbose:
        logging_level = getattr(logging, parsed_args.verbose)
        logging.basicConfig(level=logging_level, stream=sys.stdout)
    else:
        logging.basicConfig(level=logging.INFO, stream=sys.stdout)

    if parsed_args.mode:
        if parsed_args.mode[0] == '1' and os.path.exists(parsed_args.mode[1]) and ".h5ad" in parsed_args.mode[1] and os.path.exists(parsed_args.mode[2]) and ".csv" in parsed_args.mode[2]:
            if len(parsed_args.mode) == 8:
                pathName = SpatialDomainsDetectionSpaGCN(parsed_args.mode[1], parsed_args.mode[2], xarray=parsed_args.mode[3], yarray=parsed_args.mode[4], xpixel=parsed_args.mode[5], ypixel=parsed_args.mode[6], startL=parsed_args.mode[7])
            elif len(parsed_args.mode) == 7:
                pathName = SpatialDomainsDetectionSpaGCN(parsed_args.mode[1], parsed_args.mode[2], xarray=parsed_args.mode[3], yarray=parsed_args.mode[4], xpixel=parsed_args.mode[5], ypixel=parsed_args.mode[6])
            else:
                pathName = SpatialDomainsDetectionSpaGCN(parsed_args.mode[1], parsed_args.mode[2])

            logging.info(f"Done, your result file and pictures are located here -> {pathName}")
            split_string = pathName.split(" ")
            pathName = IdentifySVG(gene=parsed_args.mode[1], results=split_string[0], xarray=parsed_args.mode[3], yarray=parsed_args.mode[4], rawxarray=parsed_args.mode[3], rawyarray=parsed_args.mode[4], rawxpixel=parsed_args.mode[5], rawypixel=parsed_args.mode[6])
            logging.info(f"Done, your pictures are located here -> {pathName}")
        elif parsed_args.mode[0] == '2' and os.path.exists(parsed_args.mode[1]) and ".h5ad" in parsed_args.mode[1] and os.path.exists(parsed_args.mode[2]) and ".h5ad" in parsed_args.mode[2]:
            if len(parsed_args.mode) == 7:
                pathName = IdentifyMetaGene(gene=parsed_args.mode[1], results=parsed_args.mode[2], xpixel=parsed_args.mode[3], ypixel=parsed_args.mode[4], xarray=parsed_args.mode[5], yarray=parsed_args.mode[6], rawxarray=parsed_args.mode[3], rawyarray=parsed_args.mode[4], rawxpixel=parsed_args.mode[5], rawypixel=parsed_args.mode[6])
            else:
                pathName = IdentifyMetaGene(gene=parsed_args.mode[1], results=parsed_args.mode[2])
            logging.info(f"Done, your pictures are located here -> {pathName}")
        elif parsed_args.mode[0] == '3':
            if len(parsed_args.mode) == 5 and os.path.exists(parsed_args.mode[1]) and ".h5ad" in parsed_args.mode[1] and os.path.exists(parsed_args.mode[2]) and ".h5ad" in parsed_args.mode[2] and os.path.exists(parsed_args.mode[3]) and (".tif" in parsed_args.mode[3] or ".png" in parsed_args.mode[3] or ".jpeg" in parsed_args.mode[3]) and os.path.exists(parsed_args.mode[4]) and (".tif" in parsed_args.mode[4] or ".png" in parsed_args.mode[4] or ".jpeg" in parsed_args.mode[4]):
                pathName = MultipleTissue(parsed_args.mode[1], parsed_args.mode[2], parsed_args.mode[3], parsed_args.mode[4])
                logging.info(f"Done, your picture are located here -> {pathName}")
            elif len(parsed_args.mode) > 5:
                given_tissues = int(len(parsed_args.mode)/2 + 1)
                num_tissue = 0
                num_histology = 0
                tissues = []
                histology = []
                for i in range(1, given_tissues): #tissues
                    if os.path.exists(parsed_args.mode[i]) and ".h5ad" in parsed_args.mode[i]:
                        tissues.append(parsed_args.mode[i])
                        num_tissue=num_tissue+1
                
                for j in range(given_tissues, len(parsed_args.mode)): #histology_images
                    if os.path.exists(parsed_args.mode[j]) and (".tif" in parsed_args.mode[j] or ".png" in parsed_args.mode[j] or ".jpeg" in parsed_args.mode[j]):
                        histology.append(parsed_args.mode[j])
                        num_histology = num_histology + 1
                
                if num_histology==num_tissue and num_tissue==given_tissues and num_histology==given_tissues:
                    pathName = MultipleTissueMore(tissues, histology)
                    logging.info(f"Done, your picture are located here -> {pathName}")
                else:
                    # print error message and help if arguments are incorrect
                    PrintError(parser)
            else:
                # print error message and help if arguments are incorrect
                PrintError(parser)
        else:
            if os.path.exists(parsed_args.mode[0]) and ".h5ad" in parsed_args.mode[0] and os.path.exists(parsed_args.mode[1]) and ".csv" in parsed_args.mode[1]:
                if len(parsed_args.mode) == 7:
                    pathName = SpatialDomainsDetectionSpaGCN(parsed_args.mode[0], parsed_args.mode[1], xarray=parsed_args.mode[2], yarray=parsed_args.mode[3], xpixel=parsed_args.mode[4], ypixel=parsed_args.mode[5], startL=parsed_args.mode[6])
                elif len(parsed_args.mode) == 6:
                    pathName = SpatialDomainsDetectionSpaGCN(parsed_args.mode[0], parsed_args.mode[1], xarray=parsed_args.mode[2], yarray=parsed_args.mode[3], xpixel=parsed_args.mode[4], ypixel=parsed_args.mode[5])
                else:
                    pathName = SpatialDomainsDetectionSpaGCN(parsed_args.mode[0], parsed_args.mode[1])

                logging.info(f"Done, your result file and pictures are located here -> {pathName}")
                split_string = pathName.split(" ")
                pathName = IdentifySVG(gene=parsed_args.mode[0], results=split_string[0], xarray=parsed_args.mode[2], yarray=parsed_args.mode[3], rawxarray=parsed_args.mode[2], rawyarray=parsed_args.mode[3], rawxpixel=parsed_args.mode[4], rawypixel=parsed_args.mode[5])
                logging.info(f"Done, your pictures are located here -> {pathName}")
            else:
                # print error message and help if arguments are incorrect
                PrintError(parser)

    # check which argument was specified and perform corresponding function
    elif parsed_args.convert_h5:
        # check if there are two arguments and the first is an existing h5 file and the second is an existing txt file
        if len(parsed_args.convert_h5) == 2 and os.path.exists(parsed_args.convert_h5[0]) and ".h5" in parsed_args.convert_h5[0] and os.path.exists(parsed_args.convert_h5[1]) and ".txt" in parsed_args.convert_h5[1]:
            # perform conversion and print message with resulting file location
            pathName = ConvertH5File(parsed_args.convert_h5[0], parsed_args.convert_h5[1])
            logging.info(f"Done, your converted h5 file are located here -> {pathName}")
        else:
            # print error message and help if arguments are incorrect
            PrintError(parser)

    elif parsed_args.integrate_gene:
        # check if argument is an existing h5ad file
        if os.path.exists(parsed_args.integrate_gene) and ".h5ad" in parsed_args.integrate_gene:
            # if histology argument is also specified, check if it is an existing tif, png, or jpeg file
            if parsed_args.histology_image:
                if os.path.exists(parsed_args.histology_image) and (".tif" in parsed_args.histology_image or ".png" in parsed_args.histology_image or ".jpeg" in parsed_args.histology_image):
                    # integrate gene expression and histology into graph and print message with resulting csv file location
                    if parsed_args.pixels:
                        pathName = IntegrateIntoGraphHistology(parsed_args.integrate_gene, parsed_args.histology_image, parsed_args.pixels[0], parsed_args.pixels[1])
                    else:
                        pathName = IntegrateIntoGraphHistology(parsed_args.integrate_gene, parsed_args.histology_image)
                    logging.info(f"Done, your csv file is located here -> {pathName}")
                else:
                    # print error message and help if arguments are incorrect
                    PrintError(parser)
            else:
                # integrate gene expression into graph and print message with resulting csv file location
                if parsed_args.pixels:
                    pathName = IntegrateIntoGraph(parsed_args.integrate_gene, parsed_args.pixels[0], parsed_args.pixels[1])
                else:
                    pathName = IntegrateIntoGraph(parsed_args.integrate_gene)
                logging.info(f"Done, your csv file is located here -> {pathName}")    
        else:
            # print error message and help if arguments are incorrect
            PrintError(parser)

    elif parsed_args.spatial_domains:
        # check if arguments are an existing h5ad and csv files
        if len(parsed_args.spatial_domains) == 2 and os.path.exists(parsed_args.spatial_domains[0]) and ".h5ad" in parsed_args.spatial_domains[0] and os.path.exists(parsed_args.spatial_domains[1]) and ".csv" in parsed_args.spatial_domains[1]:
            gene = parsed_args.spatial_domains[0]
            adjCsv = parsed_args.spatial_domains[1]
            xpixel = parsed_args.pixels[0] if parsed_args.pixels else None
            ypixel = parsed_args.pixels[1] if parsed_args.pixels else None
            xarray = parsed_args.arrays[0] if parsed_args.arrays else None
            yarray = parsed_args.arrays[1] if parsed_args.arrays else None
            startL = parsed_args.start if parsed_args.start else None

            pathName = SpatialDomainsDetectionSpaGCN(gene=gene, adjCsv=adjCsv, xpixel=xpixel, ypixel=ypixel, xarray=xarray, yarray=yarray, startL=startL)
            logging.info(f"Done, your result file and pictures are located here -> {pathName}")
        else:
            # print error message and help if arguments are incorrect
            PrintError(parser)

    elif parsed_args.identify_svg:
        # check if arguments are an existing h5ad files
        if len(parsed_args.identify_svg) == 2 and os.path.exists(parsed_args.identify_svg[0]) and ".h5ad" in parsed_args.identify_svg[0] and os.path.exists(parsed_args.identify_svg[1]) and ".h5ad" in parsed_args.identify_svg[1]:
            # Identify cvg and print a message with the resulting png file location
            gene = parsed_args.identify_svg[0]
            results = parsed_args.identify_svg[1]
            xarray = parsed_args.arrays[0] if parsed_args.arrays else None
            yarray = parsed_args.arrays[1] if parsed_args.arrays else None
            rawxarray = parsed_args.raws[0] if parsed_args.raws else None
            rawyarray = parsed_args.raws[1] if parsed_args.raws else None
            rawxpixel = parsed_args.raws[2] if parsed_args.raws else None
            rawypixel = parsed_args.raws[3] if parsed_args.raws else None

            pathName = IdentifySVG(gene=gene, results=results, xarray=xarray, yarray=yarray, rawxarray=rawxarray, rawyarray=rawyarray, rawxpixel=rawxpixel, rawypixel=rawypixel)
            logging.info(f"Done, your pictures are located here -> {pathName}")
        else:
            # print error message and help if arguments are incorrect
            PrintError(parser)

    elif parsed_args.identify_meta:
        # check if arguments are an existing h5ad files
        if len(parsed_args.identify_meta) == 2 and os.path.exists(parsed_args.identify_meta[0]) and ".h5ad" in parsed_args.identify_meta[0] and os.path.exists(parsed_args.identify_meta[1]) and ".h5ad" in parsed_args.identify_meta[1]:
            # Identify meta gene and print a message with the resulting png file location
            gene = parsed_args.identify_meta[0]
            results = parsed_args.identify_meta[1]
            xpixel = parsed_args.pixels[0] if parsed_args.pixels else None
            ypixel = parsed_args.pixels[1] if parsed_args.pixels else None
            xarray = parsed_args.arrays[0] if parsed_args.arrays else None
            yarray = parsed_args.arrays[1] if parsed_args.arrays else None
            rawxarray = parsed_args.raws[0] if parsed_args.raws else None
            rawyarray = parsed_args.raws[1] if parsed_args.raws else None
            rawxpixel = parsed_args.raws[2] if parsed_args.raws else None
            rawypixel = parsed_args.raws[3] if parsed_args.raws else None

            pathName = IdentifyMetaGene(gene=gene, results=results, xpixel=xpixel, ypixel=ypixel, xarray=xarray, yarray=yarray, rawxarray=rawxarray, rawyarray=rawyarray, rawxpixel=rawxpixel, rawypixel=rawypixel)
            logging.info(f"Done, your pictures are located here -> {pathName}")
        else:
            # print error message and help if arguments are incorrect
            PrintError(parser)

    elif parsed_args.multiple_tissue:
        # check if arguments are an existing h5ad files and tif/jpeg/jpg files
        if len(parsed_args.multiple_tissue) == 4 and os.path.exists(parsed_args.multiple_tissue[0]) and ".h5ad" in parsed_args.multiple_tissue[0] and os.path.exists(parsed_args.multiple_tissue[1]) and ".h5ad" in parsed_args.multiple_tissue[1] and os.path.exists(parsed_args.multiple_tissue[2]) and (".tif" in parsed_args.multiple_tissue[2] or ".png" in parsed_args.multiple_tissue[2] or ".jpeg" in parsed_args.multiple_tissue[2]) and os.path.exists(parsed_args.multiple_tissue[3]) and (".tif" in parsed_args.multiple_tissue[3] or ".png" in parsed_args.multiple_tissue[3] or ".jpeg" in parsed_args.multiple_tissue[3]):
            pathName = MultipleTissue(parsed_args.multiple_tissue[0], parsed_args.multiple_tissue[1], parsed_args.multiple_tissue[2], parsed_args.multiple_tissue[3])
            logging.info(f"Done, your picture are located here -> {pathName}")
        elif len(parsed_args.multiple_tissue) > 4:
            given_tissues = int(len(parsed_args.multiple_tissue)/2)
            num_tissue = 0
            num_histology = 0
            tissues = []
            histology = []
            for i in range(0, given_tissues): #tissues
                if os.path.exists(parsed_args.multiple_tissue[i]) and ".h5ad" in parsed_args.multiple_tissue[i]:
                    tissues.append(parsed_args.multiple_tissue[i])
                    num_tissue=num_tissue+1
            
            for j in range(given_tissues, len(parsed_args.multiple_tissue)): #histology_images
                if os.path.exists(parsed_args.multiple_tissue[j]) and (".tif" in parsed_args.multiple_tissue[j] or ".png" in parsed_args.multiple_tissue[j] or ".jpeg" in parsed_args.multiple_tissue[j]):
                    histology.append(parsed_args.multiple_tissue[j])
                    num_histology = num_histology + 1
            
            if num_histology==num_tissue and num_tissue==given_tissues and num_histology==given_tissues:
                pathName = MultipleTissueMore(tissues, histology)
                logging.info(f"Done, your picture are located here -> {pathName}")
            else:
                # print error message and help if arguments are incorrect
                PrintError(parser)
        else:
            # print error message and help if arguments are incorrect
            PrintError(parser)

    elif parsed_args.read_keys:
        if os.path.exists(parsed_args.read_keys) and ".h5ad" in parsed_args.read_keys:
            ReadKeys(parsed_args.read_keys)
        else:
            # print error message and help if arguments are incorrect
            PrintError(parser)

    elif parsed_args.read_specific_keys:
        if len(parsed_args.read_specific_keys) == 2 and os.path.exists(parsed_args.read_specific_keys[0]) and ".h5ad" in parsed_args.read_specific_keys[0]:
            ReadSpecificKeys(parsed_args.read_specific_keys[0], parsed_args.read_specific_keys[1])
        else:
            # print error message and help if arguments are incorrect
            PrintError(parser)