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
parser.add_argument("--mode", nargs="+", help="Modes on which these programs work.(Mode 1: integrate gene & spatial domains(default), Mode 2: identify svg & identify meta, Mode 3: multiple tissue)", default='')
#parser.add_argument("--convert_h5", nargs=2, help="Read original 10x_h5 data and save it to h5ad. The path to the h5 file is required.", default='')
parser.add_argument("--start", help="Start value for search l. It is used in combination with command: spatial_domains.", type=float, default=0.01)
parser.add_argument("--pixels", nargs=2, help="The x and y coordinates for pixels are typed here. It is used in combination with commands: integrate_gene, spatial_domains, identify_meta.", default='')
parser.add_argument("--arrays", nargs=2, help="The x and y coordinates for arrays are typed here. It is used in combination with commands: spatial_domains, identify_svg, identify_meta.", default='')
parser.add_argument("--raws", nargs=4, help="The x array, y array, x pixel and y pixel coordinates for rows are typed here. It is used in combination with commands: identify_svg, identify_meta.", default='')
parser.add_argument("--read_obs", help="Print all keys found in your file.", default='')
parser.add_argument("--read_specific_obs", nargs=2, help="Print all specific keys found in your file.", default='')
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
        if parsed_args.mode[0] == '1':
            gene = parsed_args.mode[1]
            histology = parsed_args.mode[2]
            xpixel = parsed_args.pixels[0] if parsed_args.pixels else "null"
            ypixel = parsed_args.pixels[1] if parsed_args.pixels else "null"
            xarray = parsed_args.arrays[0] if parsed_args.arrays else "null"
            yarray = parsed_args.arrays[1] if parsed_args.arrays else "null"
            startL = parsed_args.start if parsed_args.start else 0.01

            if os.path.exists(gene) and ".h5ad" in gene and os.path.exists(gene) and (".tif" in histology or ".png" in histology or ".jpeg" in histology):
                csvPathName = IntegrateIntoGraphHistology(gene=gene, histology=histology, xpixel=xpixel, ypixel=ypixel)
            elif os.path.exists(gene) and ".h5ad" in gene and os.path.exists(gene):
                csvPathName = IntegrateIntoGraph(gene=gene, xpixel=xpixel, ypixel=ypixel)
            else:
                PrintError(parser)
                sys.exit(1)

            logging.info(f"Done, your csv file is located here -> {csvPathName}")

            pathName = SpatialDomainsDetectionSpaGCN(gene=gene, adjCsv=csvPathName, xarray=xarray, yarray=yarray, xpixel=xpixel, ypixel=ypixel, startL=startL)

            logging.info(f"Done, your result file and pictures are located here -> {pathName}")
        elif parsed_args.mode[0] == '2' and os.path.exists(parsed_args.mode[1]) and ".h5ad" in parsed_args.mode[1] and os.path.exists(parsed_args.mode[2]) and ".h5ad" in parsed_args.mode[2]:
            gene = parsed_args.mode[1]
            results = parsed_args.mode[2]
            xpixel = parsed_args.pixels[0] if parsed_args.pixels else "null"
            ypixel = parsed_args.pixels[1] if parsed_args.pixels else "null"
            xarray = parsed_args.arrays[0] if parsed_args.arrays else "null"
            yarray = parsed_args.arrays[1] if parsed_args.arrays else "null"
            rawxarray = parsed_args.raws[0] if parsed_args.raws else "null"
            rawyarray = parsed_args.raws[1] if parsed_args.raws else "null"
            rawxpixel = parsed_args.raws[2] if parsed_args.raws else "null"
            rawypixel = parsed_args.raws[3] if parsed_args.raws else "null"

            if os.path.exists(gene) and ".h5ad" in gene and os.path.exists(results) and ".h5ad" in results:
                pathName = IdentifySVG(gene=gene, results=results, xarray=xarray, yarray=yarray, rawxarray=rawxarray, rawyarray=rawyarray, rawxpixel=rawxpixel, rawypixel=rawypixel)
                logging.info(f"Done, your pictures(identify svg) are located here -> {pathName}")

                pathName = IdentifyMetaGene(gene=gene, results=results, xpixel=xpixel, ypixel=ypixel, xarray=xarray, yarray=yarray, rawxarray=rawxarray, rawyarray=rawyarray, rawxpixel=rawxpixel, rawypixel=rawypixel)
                logging.info(f"Done, your pictures(identify meta) are located here -> {pathName}")
            else:
                PrintError(parser)
                sys.exit(1)

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
            gene = parsed_args.mode[0]
            histology = parsed_args.mode[1]
            xpixel = parsed_args.pixels[0] if parsed_args.pixels else "null"
            ypixel = parsed_args.pixels[1] if parsed_args.pixels else "null"
            xarray = parsed_args.arrays[0] if parsed_args.arrays else "null"
            yarray = parsed_args.arrays[1] if parsed_args.arrays else "null"
            startL = parsed_args.start if parsed_args.start else 0.01

            if os.path.exists(gene) and ".h5ad" in gene and os.path.exists(gene) and (".tif" in histology or ".png" in histology or ".jpeg" in histology):
                csvPathName = IntegrateIntoGraphHistology(gene=gene, histology=histology, xpixel=xpixel, ypixel=ypixel)
            elif os.path.exists(gene) and ".h5ad" in gene and os.path.exists(gene):
                csvPathName = IntegrateIntoGraph(gene=gene, xpixel=xpixel, ypixel=ypixel)
            else:
                PrintError(parser)
                sys.exit(1)

            logging.info(f"Done, your csv file is located here -> {csvPathName}")

            pathName = SpatialDomainsDetectionSpaGCN(gene=gene, adjCsv=csvPathName, xarray=xarray, yarray=yarray, xpixel=xpixel, ypixel=ypixel, startL=startL)

            logging.info(f"Done, your result file and pictures are located here -> {pathName}")
    
    # check which argument was specified and perform corresponding function
    #elif parsed_args.convert_h5:
        # check if there are two arguments and the first is an existing h5 file and the second is an existing txt file
    #    if len(parsed_args.convert_h5) == 2 and os.path.exists(parsed_args.convert_h5[0]) and ".h5" in parsed_args.convert_h5[0] and os.path.exists(parsed_args.convert_h5[1]) and ".txt" in parsed_args.convert_h5[1]:
            # perform conversion and print message with resulting file location
    #        pathName = ConvertH5File(parsed_args.convert_h5[0], parsed_args.convert_h5[1])
    #        logging.info(f"Done, your converted h5 file are located here -> {pathName}")
    #    else:
            # print error message and help if arguments are incorrect
    #        PrintError(parser)

    elif parsed_args.read_obs:
        if os.path.exists(parsed_args.read_obs) and ".h5ad" in parsed_args.read_obs:
            ReadKeys(parsed_args.read_obs)
        else:
            # print error message and help if arguments are incorrect
            PrintError(parser)

    elif parsed_args.read_specific_obs:
        if len(parsed_args.read_specific_obs) == 2 and os.path.exists(parsed_args.read_specific_obs[0]) and ".h5ad" in parsed_args.read_specific_obs[0]:
            ReadSpecificKeys(parsed_args.read_specific_obs[0], parsed_args.read_specific_obs[1])
        else:
            # print error message and help if arguments are incorrect
            PrintError(parser)