#!/usr/bin/env python

'''This program should take a BAM file produced by Bismark and return a table suitable for histogram generation where the columns represent the counts per bp: # methylated C, # unmethylated C, no C. 
There will be two output files: one for reads of 100 bp length, the other for reads of 200 bp length.

usage: $ bam_to_meCountHist.py -in test.bam -o testHistograms 
'''

import sys
import argparse
import getopt
import pysam
import os
import doctest
import gzip
import numpy as np


def get_args():
    parser=argparse.ArgumentParser(description='Count number of Cs etc. per read position')
    parser.add_argument('--BAMfile', '-in', type=str, required=True, help="BAM file")
    parser.add_argument('--outfile', '-o', type=str, required=True, help = 'Prefix for the output files of methylation distribution')
    parser.add_argument('--CpGonly', action='store_true', default=False, help = 'If set, only cytosines in CpG contexts are taken into account.')
    parser.add_argument('--trimStart', '-strim', type=int, required=False, default = 0, help = 'Number indicating how many bp should be ignored at the 5 prime end of the read.')
    parser.add_argument('--trimEnd', '-etrim', type=int, required=False, default = 0, help = 'Number indicating how many bp should be ignored at the 3 prime end of the read.')
    parser.add_argument('--checkStrand', action='store_true', default=False, help = 'If set, only reads mapping to the forward strand are taken into account.')
    
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


def counting(Read_seq, count_object, CpG_only):
        
    if count_object.ndim != 2 or len(count_object) != 3:
        raise NameError("The count object does not have the correct number of dimensions.")
    
    if len(count_object[1,]) != len(Read_seq):
        raise NameError("The numpy array count object does not have the appropriate number of slots for the length of the read.")
    
    i = 0
    if not CpG_only:
        for bp in Read_seq:
            if bp.isupper():
                count_object[0,i] += 1 # methylated Cs
            elif bp.islower():
                count_object[1,i] += 1 # unmethylated Cs
            elif bp == ".": 
                count_object[2,i] += 1 # not a C
            i +=1
    else:
        for bp in Read_seq:
            if bp == "Z":
                count_object[0,i] += 1
            elif bp == "z":
                count_object[1,i] += 1
            elif bp == "." or bp.isupper() or bp.islower:
                count_object[2,i] += 1
            i +=1
    
    return(count_object)
    



def main():
    args = get_args()
    
    infile = pysam.Samfile(args.BAMfile, "rb")

    np100 = np.zeros([3, 100], float)
    np200 = np.zeros([3, 200], float)
    
    i100 = 0
    i200 = 0
        
    for DNAread in infile:
        
        # ignore reverse strand reads if the checkStrand option is set
        if args.checkStrand and DNAread.is_reverse:
            continue            
            
        me_seq = get_read_seq(DNAread, args.trimStart, args.trimEnd) 
        
        if len(me_seq) == 100:
            np100 = counting(me_seq, np100, CpG_only = args.CpGonly)
            i100 += 1
        elif len(me_seq) == 200:
            np200 = counting(me_seq, np200, CpG_only = args.CpGonly)
            i200 += 1
    
    CountOut="Reads of length 100bp: %d. Reads of length 200bp: %d" % (i100, i200)
    print CountOut
    
    np.savetxt(args.outfile + '_100bpReads.csv', np100.transpose(), delimiter=',', fmt="%12.6G")
    np.savetxt(args.outfile + '_200bpReads.csv', np200.transpose(), delimiter=',', fmt="%12.6G")
    
    
if __name__ == '__main__':
    main()
    
