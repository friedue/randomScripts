#!/usr/bin/env python

'''This script generates a list of CpG pairs based on pileOMeth`s bedGraph
 file of methylation calls. 
 The output is a tab-separated list of CpG pairs with 
 one line per pair: chr, start, end.
'''

import sys
import argparse
import getopt
import gzip

def get_args():
    parser=argparse.ArgumentParser(description='Generate a list of adjacent CpGs based on pileOMeth output.')
    parser.add_argument('--CpGfile', '-i', type=str, required=True, help="The bedGraph output of pileOMeth with chromosome, CpG position, DNA me percent, reads with C, reads with T. The program expects the output in an _unstranded_ fashion, i.e., pileOMeth must have been used with --mergeContext.")
    parser.add_argument('--outfile', '-o', type=str, required=True, help = 'Prefix for the output file (list of CpG pairs)')
    parser.add_argument('--tabSeparated', '-tab', action='store_true', default=False, help = 'Set this if the file is tab, rather than space separated.')
    parser.add_argument('--maxDist', '-d', type=int, required=False, default = 200, help = 'Maximum distance between adjacent CpG; usually the max. read length.')
    
    args=parser.parse_args()
    return args
    
def main():
    args = get_args()

    if ".gz" in args.CpGfile:
        infile = gzip.open(args.CpGfile, "r")
    else:
        infile = open(args.CpGfile, "r")
        
    out = open(args.outfile + '.bed', 'wb')
    
    max_dist = args.maxDist
    
    prev_chrom = 0
    prev_cpg = max_dist * -2

    for Line in infile.readlines():
        if args.tabSeparated:
            Line = Line.strip("\n")
            Line = Line.split("\t")
        else:
            Line = Line.split()
            
        chrom, cpg1, cpg2, coverage = Line[0], int(Line[1]), int(Line[2]), int(Line[4]) + int(Line[5])
        ## could add a filter based on coverage here
            
        if chrom == prev_chrom:
                
            dist = cpg1 - prev_cpg
               
            if 1 < dist <= max_dist:
                out.write("%s\t%d\t%d\n" % (chrom, prev_cpg, cpg2))
            
        prev_chrom = chrom
        prev_cpg = cpg1

    out.close()


if __name__ == '__main__':
    main()
