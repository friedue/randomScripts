#!/usr/bin/env python

'''This program should take a BAM file produced by Bismark and return a table where for each read the following values are indicated: %meC, # methylated C, # unmethylated C, read length. '''

import sys
import argparse
import getopt
import pysam
import os
import doctest
import gzip


def get_args():
    parser=argparse.ArgumentParser(description='Split suspender BAM files')
    parser.add_argument('--BAMfile', '-in', type=str, required=True, help="BAM file")
    parser.add_argument('--outfile', '-o', type=str, required=True, help = 'Prefix for the output file of counts per reads')
    parser.add_argument('--trimStart', '-strim', type=int, required=False, default = 0, help = 'Number indicating how many bp should be ignored at the 5 prime end of the read.')
    parser.add_argument('--trimEnd', '-etrim', type=int, required=False, default = 0, help = 'Number indicating how many bp should be ignored at the 3 prime end of the read.')
    
    args=parser.parse_args()
    return args

def get_read_seq(Read):
    '''returns the sequence part of Bismark's XM tag'''
    Seq = [item for item in Read.tags if item[0] == 'XM'][0][1]
    
    return Seq



def get_mecounts(Read_seq, strim, etrim):
    '''counts the (un)methylated residues in a read'''
    Counts = [0,0,0]
    
    beg = 0 + strim
    if etrim > 0:
        end = -etrim
    else:
        end = None
        
    Read = Read_seq[beg:end]
   
    for bp in Read:
        if bp.isupper():
            Counts[0] += 1 # methylated Cs
        if bp.islower():
            Counts[1] += 1 # unmethylated Cs
        elif bp == ".": 
            Counts[2] += 1 # not a C
    
    if sum(Counts) != len(Read):
            raise NameError("Individual counts (%d) do not match the length of the (trimmed) read sequence (%d)" % (sum(Counts), len(Read)))
    
    return Counts

def get_mepercent(C):
    
    if C[1] > 0:
        meperc = float(C[0])/float(C[1]) * 100
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
        me_seq = get_read_seq(DNAread)
        c = get_mecounts(me_seq, args.trimStart, args.trimEnd)
        me_percent = get_mepercent(c)

        fo.write('{}\t{}\t{}\t{}\t{}\n'.format(DNAread.query_name, me_percent, c[0], c[1], sum(c)))

    fo.close() 
    
if __name__ == '__main__':
    main()
    
