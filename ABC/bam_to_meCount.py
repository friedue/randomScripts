#!/usr/bin/env python

'''This program should take a BAM file produced by Bismark and return a table where for each read the following values are indicated: # methylated C, # unmethylated C, # methylated other C, # unmethylated other C, # not a C, read length. '''

import sys
import argparse
import getopt
import pysam
import os
import doctest
import gzip


def get_args():
    parser=argparse.ArgumentParser(description='Count numbers of methylated and unmethylated Cs per read.'
    'The output table contains: #meCpgG, #unmeCpG, #meC-otherContext, #unmeC-otherContext,read length')
    parser.add_argument('--BAMfile', '-in', type=str, required=True, help="BAM file")
    parser.add_argument('--outfile', '-o', type=str, required=True, help = 'Prefix for the output file of counts per reads')
    parser.add_argument('--trimStart', '-strim', type=int, required=False, default = 0, help = 'Number indicating how many bp should be ignored at the 5 prime end of the read.')
    parser.add_argument('--trimEnd', '-etrim', type=int, required=False, default = 0, help = 'Number indicating how many bp should be ignored at the 3 prime end of the read.')
    
    args=parser.parse_args()
    return args

def get_read_seq(Read, strim, etrim):
    '''returns the sequence part of Bismark's XM tag; if necessary,
    the sequence is trimmed'''
    
    Seq = [item for item in Read.tags if item[0] == 'XM'][0][1]
    
    beg = 0 + strim
    if etrim > 0:
        end = -etrim
    else:
        end = None
        
    return Seq[beg:end]



def get_mecounts(Read_seq):
    '''counts the (un)methylated residues in a read'''
    Counts = [0,0,0,0,0]

    for bp in Read_seq:
        if bp == "Z": 
            Counts[0] += 1 # methylated CpG
        elif bp == "z":
            Counts[1] += 1 # unmethylated CpG
        elif bp.isupper():
            Counts[2] += 1 # methylated other C
        elif bp.islower():
            Counts[3] += 1 # unmethylated other C
        elif bp == ".":
            Counts[4] += 1 # not a C
    
    if sum(Counts) != len(Read_seq):
        raise StandardError("Individual counts (%d) do not match the length of the (trimmed) read sequence (%d)" % (sum(Counts), len(Read_seq)))
    
    return Counts

def get_mepercent(C):
    
    if C[1] > 0:
        meperc = float(C[0])/float(C[1] + C[0]) * 100
    elif C[0] == 0:
        meperc = 'NA'
    else:
        meperc = 100
        
    return meperc

def main():
    args = get_args()
    
    infile = pysam.Samfile(args.BAMfile, "rb")

    fo = gzip.open(args.outfile + '.tab.gz', 'wb')
        
    for DNAread in infile:
        me_seq = get_read_seq(DNAread, args.trimStart, args.trimEnd)
        c = get_mecounts(me_seq)
      #  me_percent = get_mepercent(c)

        fo.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(DNAread.query_name, c[0], c[1], c[2], c[3], c[4], sum(c)))

    fo.close() 
    
if __name__ == '__main__':
    main()
    
