# -*- coding: utf-8 -*-

# Author: Pedro Furió Tarí - then fixed for python3 by ALB
# Data: 30/05/2013
#

# The input should have the following format:
#############################################
# chr6	47788069	47788180	RNY1	.	-
# chr6	47781624	47781725	RNY3	.	-

# The output will be named the following way:
# chr6	RNY1	featureName	47788069	47788180	.	-	.	feature_id "RNY1_chr6_47788069_47788180";
# chr6	RNY3	featureName	47781624	47781725	.	-	.	feature_id "RNY3_chr6_47781624_47781725";

import getopt, sys, os.path


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:f:o:b:", ["help", "input=","featureName=","output=","biotype="])
    except getopt.GetoptError as err:
        print(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    infile = None
    outfile = None
    featName = "gene"
    genetype = "protein_coding"

    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("-i", "--input"):
            if os.path.isfile(arg):
                infile = arg
        elif opt in ("-o", "--output"):
            outfile = arg
        elif opt in ("-f", "--featureName"):
            featName = arg
        elif opt in ("-b", "--biotype"):
            genetype = arg
        else:
            assert False, "Unhandled option"

    if infile is not None and outfile is not None:
        run(infile, outfile,featName,genetype)
    else:
        usage()


def usage():
    print("\nUsage: python bed2gtf [options] <mandatory>")
    print("Options:")
    print("\t-h, --help:\n\t\t show this help message and exit")
    print("Mandatory:")
    print("\t-i, --input:\n\t\t File with the regions in bed format")
    print("\t-o, --output:\n\t\t Name of the gtf file output file. Directory where the file will be created should exist!")
    print("Optional:")
    print("\t-f, --featureName:\n\t\t This is what you're going to call the feature,default gene (exon, gene, peak, potato whatever you like)")
    print("\t-b, --biotype:\n\t\t This is what you're going to call the gene biotype")

    print("\n30/05/2013. Pedro Furió Tarí. - then modified by ALB 28/01/2020\n")

def run(infile, outfile, featureNameUser, biotype):
    metaFeatureName = "gene_id "

    inf  = open(infile, 'r')
    outf = open(outfile,'w')
    featureName = featureNameUser
    cont = 1
    for linea in inf:
        linea_split = linea.split()
        chrom = linea_split[0]
        ini_pos = int(linea_split[1])
        fin_pos = int(linea_split[2])
        name = str(linea_split[3])
        strand = str(linea_split[5])
        #peak = linea_split[3]

        #outf.write(chrom + "\tpfurio\tpeak\t" + str(ini_pos) + "\t" + str(fin_pos) + '\t.\t+\t.\tpeak_id "' + peak + '";\n')
        outf.write(chrom + "\t" + name + "\t" + featureNameUser + "\t" + str(ini_pos) + "\t" + str(fin_pos) + '\t.\t' + strand + "\t.\t" + metaFeatureName + '"' + name + "_" + chrom + "_" + str(ini_pos) + "_" + str(fin_pos) + "; gene_biotype " + '"' + str(biotype) +  '"' + ";\n")

        cont += 1

    inf.close()
    outf.close()


if __name__ == "__main__":
    main()
