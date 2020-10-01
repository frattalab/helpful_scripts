#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import argparse

def cmdline_args():
        # Make parser object
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument("-f", "--folder", help="folder containing splicejunctions")

    p.add_argument("-s", "--suffix", help="suffix of the splicejunctions",
                   default=".SJ.out.tab")

    p.add_argument("-c", "--column", help="column to sum use 4 for beds,6 for STar",
                   default=4,type=int)

    p.add_argument("-o", "--output", help="output file name",
                   default="splice_junction_sums.txt")


    return(p.parse_args())

def main():
    """ Main program """
    args = cmdline_args()
    #collect the names of the files in the given folder
    files_to_count = args.folder
    #take only the files which contain the suffix
    onlysuffix = [f for f in os.listdir(files_to_count) if args.suffix in f]
    #loop through all the files and get the sum of column 4(counting frm 0)
    count_column = args.column

    output_file = open(args.output, 'w')
    for file in onlysuffix:
        file_name = file.replace(args.suffix,'')
        with open(os.path.join(files_to_count,file)) as infile:
            try:
                tot = sum(int(line.split("\t")[count_column]) for line in infile)
            except ValueError:
                print("%s does not have a column %d" % (file_name,count_column))
        print([file_name,tot])
        output_file.write("%s,%d\n" % (file_name,tot))
    output_file.close()

    return(0)

if __name__ == "__main__":
    main()
