import argparse
import os
from functions import *

parser = argparse.ArgumentParser(prog='Computational Genomics Project SpaGCN', description='This program does multiple things, it can convert h5 file to h5ad. It does Integrate gene expression and histology into a Graph, Spatial domain detection using SpaGCN, Identify SVGs, Identify Meta Gene and Multiple tissue sections analysis.')

parser.add_argument("-convert_h5", nargs=2, help="Read original 10x_h5 data and save it to h5ad. The path to the h5 file is required.", default='')
parser.add_argument("-integrate_gene", help="Integrate gene expression into a Graph. The path to the h5ad file is required", default='')
parser.add_argument("-histology", help="Integrate and histology into a Graph. The path to the .tif file is required", default='')
parser.add_argument("-spatial_domains", nargs=2, help="Spatial domain detection using SpaGCN. Paths to h5ad and csv files are required.", default='')
parser.add_argument("-identify_csv", nargs=2, help="Identify SVGs. Paths to h5ad gene matrix and h5ad results files are required.", default='')
parser.add_argument("-identify_meta", nargs=2, help="Identify Meta Gene. Paths to h5ad gene matrix and h5ad results files are required.", default='')
parser.add_argument("-multiple_tissue", nargs=4, help="Multiple tissue sections analysis.Paths to h5ad first tissue, h5ad second tissue, tif first tissue and tif second tissue files are required.", default='')
parser.add_argument('--version', action='version', version='%(prog)s v0.9')
parsed_args = parser.parse_args()

if parsed_args.convert_h5:
    if len(parsed_args.convert_h5) == 2 and os.path.exists(parsed_args.convert_h5[0]) and ".h5" in parsed_args.convert_h5[0] and os.path.exists(parsed_args.convert_h5[1]) and ".txt" in parsed_args.convert_h5[1]:
        pathName = ConvertH5File(parsed_args.convert_h5[0], parsed_args.convert_h5[1])
        print("Done, your converted h5 file are located here -> " + pathName)
    else:
        print("Your files do not exist or are of the wrong type.")
        parser.print_help()

elif parsed_args.integrate_gene:
    if os.path.exists(parsed_args.integrate_gene) and ".h5ad" in parsed_args.integrate_gene:
        if parsed_args.histology:
            if os.path.exists(parsed_args.histology) and (".tif" in parsed_args.histology or ".png" in parsed_args.histology or ".jpeg" in parsed_args.histology):
                pathName = IntegrateIntoGraphHistology(parsed_args.integrate_gene, parsed_args.histology)
            else:
                print("Your histology file does not exist or is of the wrong type.")
        else:
            pathName = IntegrateIntoGraph(parsed_args.integrate_gene)
        print("Done, your csv file is located here -> " + pathName)
    else:
        print("Your gene expression file does not exist or is of the wrong type.")
        parser.print_help()
        
elif parsed_args.spatial_domains:
    if len(parsed_args.spatial_domains) == 2 and os.path.exists(parsed_args.spatial_domains[0]) and ".h5ad" in parsed_args.spatial_domains[0] and os.path.exists(parsed_args.spatial_domains[1]) and ".csv" in parsed_args.spatial_domains[1]:
        pathName = SpatialDomainsDetectionSpaGCN(parsed_args.spatial_domains[0], parsed_args.spatial_domains[1])
        print("Done, your result file and pictures are located here -> " + pathName)
    else:
        print("Your files do not exist or are of the wrong type.")
        parser.print_help()
        
elif parsed_args.identify_csv:
    if len(parsed_args.identify_csv) == 2 and os.path.exists(parsed_args.identify_csv[0]) and ".h5ad" in parsed_args.identify_csv[0] and os.path.exists(parsed_args.identify_csv[1]) and ".h5ad" in parsed_args.identify_csv[1]:
        pathName = IdentifyCSV(parsed_args.identify_csv[0], parsed_args.identify_csv[1])
        print("Done, your pictures are located here -> " + pathName)
    else:
        print("Your files do not exist or are of the wrong type.")
        parser.print_help()
        
elif parsed_args.identify_meta:
    if len(parsed_args.identify_meta) == 2 and os.path.exists(parsed_args.identify_meta[0]) and ".h5ad" in parsed_args.identify_meta[0] and os.path.exists(parsed_args.identify_meta[1]) and ".h5ad" in parsed_args.identify_meta[1]:
        pathName = IdentifyMetaGene(parsed_args.identify_meta[0], parsed_args.identify_meta[1])
        print("Done, your pictures are located here -> " + pathName)
    else:
        print("Your files do not exist or are of the wrong type.")
        parser.print_help()
        
elif parsed_args.multiple_tissue:
    if len(parsed_args.multiple_tissue) == 4 and os.path.exists(parsed_args.multiple_tissue[0]) and ".h5ad" in parsed_args.multiple_tissue[0] and os.path.exists(parsed_args.multiple_tissue[1]) and ".h5ad" in parsed_args.multiple_tissue[1] and os.path.exists(parsed_args.multiple_tissue[2]) and ".tif" in parsed_args.multiple_tissue[2] and os.path.exists(parsed_args.multiple_tissue[3]) and ".tif" in parsed_args.multiple_tissue[3]:
        pathName = MultipleTissue(parsed_args.multiple_tissue[0], parsed_args.multiple_tissue[1], parsed_args.multiple_tissue[2], parsed_args.multiple_tissue[3])
        print("Done, your picture are located here -> " + pathName)
    else:
        print("Your files do not exist or are of the wrong type.")
        parser.print_help()