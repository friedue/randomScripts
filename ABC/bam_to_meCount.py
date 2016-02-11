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
    parser=argparse.ArgumentParser(description='Count numbers of methylated and unmethylated Cs per read')
    parser.add_argument('--BAMfile', '-in', type=str, required=True, help="BAM file")
    parser.add_argument('--outfile', '-o', type=str, required=True, help = 'Prefix for the output file of counts per reads')
    parser.add_argument('--trimStart', '-strim', type=int, required=False, default = 0, help = 'Number indicating how many bp should be ignored at the 5 prime end of the read.')
    parser.add_argument('--trimEnd', '-etrim', type=int, required=False, default = 0, help = 'Number indicating how many bp should be ignored at the 3 prime end of the read.')
    parser.add_argument('--CpGonly', action='store_true', default=False, help = 'If set, only cytosines in CpG contexts are taken into account.')
    
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



def get_mecounts(Read_seq, CpG_only):
    '''counts the (un)methylated residues in a read'''
    Counts = [0,0,0]
    
    if not CpG_only:
        for bp in Read_seq:
            if bp.isupper():
                Counts[0] += 1 # methylated Cs
            elif bp.islower():
                Counts[1] += 1 # unmethylated Cs
            elif bp == ".": 
                Counts[2] += 1 # not a C
    else:
        for bp in Read_seq:
            if bp == "Z":
                Counts[0] += 1
            elif bp == "z":
                Counts[1] += 1
            elif bp == "." or bp.isupper() or bp.islower:
                Counts[2] += 1
    
    if sum(Counts) != len(Read_seq):
            raise NameError("Individual counts (%d) do not match the length of the (trimmed) read sequence (%d)" % (sum(Counts), len(Read_seq)))
    
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
        c = get_mecounts(me_seq, CpG_only = args.CpGonly)
        me_percent = get_mepercent(c)

        fo.write('{}\t{}\t{}\t{}\t{}\n'.format(DNAread.query_name, me_percent, c[0], c[1], sum(c)))

    fo.close() 
    
if __name__ == '__main__':
    main()
    
