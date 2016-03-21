#!/usr/bin/env python

'''This program should take a BAM file produced by Bismark and a list of adjacent CpG sites for which the methylation statuses will be reported in a pairwise manner. '''

import sys
import argparse
import getopt
import pysam
import os
import doctest
import gzip
import itertools


def get_args():
    parser=argparse.ArgumentParser(description='Count numbers of methylated and unmethylated Cs per read.'
    'The output table contains: #meCpgG, #unmeCpG, #meC-otherContext, #unmeC-otherContext,read length')
    parser.add_argument('--BAMfile', '-bam', type = str, required=True, help="BAM file")
    parser.add_argument('--CpGpairs', '-pairs', type = str, required=True, help="The bedGraph output of pileOMeth with chromosome, CpG position, %DNA me, reads with C, reads with T.")
    parser.add_argument('--minMapQual', '-mmq', type = int, default=0, help="Min. mapping quality accepted for the reads that will be used to count the methylation state ocurrences. See http://bit.ly/25glGcI for information about the different aligners' MapQ calculations.")
    parser.add_argument('--minCoverage', '-mc', type = int, default = 0, help = 'Indicate the minimum number of reads that must overlap with each adjacent pair.')
    parser.add_argument('--outfile', '-out', type = str, required=True, help = 'Prefix for the output file')
    
    args=parser.parse_args()
    return args


def get_CIGARbased_sequence(readinfo):
    
    ori_seq = [item for item in readinfo.tags if item[0] == 'XM'][0][1]
    ct = readinfo.cigartuples
    new_seq = ''
    start_pos = 0
    
    '''see http://pysam.readthedocs.org/en/latest/api.html#pysam.AlignedSegment.cigartuples for information about how the integers are mapped to the CIGAR operations'''
    for c in ct:
        # match/mismatch
        if c[0] == 0 or c[0]== 7 or c[0] == 8: 
            end_pos = start_pos + c[1]
            new_seq = new_seq + ori_seq[start_pos:end_pos]
            start_pos = end_pos
            
        # deletion in read --> insertion necessary to make it fit with the reference
        elif c[0] == 2 or c[0] == 3 or c[0] == 5: 
            new_seq = new_seq + '_'*c[1]
            
        # insertion in read compared to reference sequence --> must be deleted from the read sequence
        elif c[0] == 1: 
            start_pos =+ c[1]
    
    return new_seq



def main():
    
    args = get_args()
    
    # read in CpG pair information    
    pairsfile = open(args.CpGpairs, "r")
    pairs_coords = []
    for Line in pairsfile.readlines(): #assuming tab-separated output from get_pairs.py
        Line = Line.strip("\n")
        Line = Line.split("\t")
        pairs_coords.append( (Line[0], int(Line[1]), int(Line[2])) )
    
    bc = []
    for i in range(0, len(pairs_coords)):
        bc.append([0,0,0,0])
    
    pairs_dict = dict(itertools.izip(pairs_coords, bc ))
    
    # get the reads
    bamfile = pysam.Samfile(args.BAMfile, "rb")
    
    # for each read, check which adjacent pairs it overlaps and extract
    # the methylation state information for those pairs
    for Read in bamfile:
        chrom = Read.reference_name
        r_start = int(Read.reference_start)
        r_end = int(Read.reference_end)
        
        # get pairs within the read's range
        elig_pairs = [Pair for Pair in pairs_coords if Pair[0] == chrom and Pair[1] >= r_start and Pair[2] <= r_end]
        
        if len(elig_pairs) > 0 and Read.mapping_quality >= args.minMapQual:
            # first, check CIGAR string for coordinate-altering operations
            cs = Read.cigarstring
            if 'D' in cs or 'I' in cs or 'N' in cs or 'H' in cs or 'P' in cs:
                bs_seq = get_CIGARbased_sequence(Read)
            else:
                bs_seq = [item for item in Read.tags if item[0] == 'XM'][0][1]
            
            for p in elig_pairs:
                if not Read.is_reverse:
                    b1 = p[1] - r_start
                    b2 = p[2] - 1 - r_start
                    
                else:
                    b1 = p[1] + 1 - r_start
                    b2 = p[2] - r_start
                    
                state = bs_seq[b1] + bs_seq[b2-1]
                
                if '.' in state:
                    raise StandardError("Did not find a z or Z at the expected position (%d, %d) within read %s" % (p[1], p[2], Read.query_name))
                
                # record state in temporary dictionary
                sdc = dict(itertools.izip(['ZZ','Zz','zZ','zz'], [0,0,0,0]))
                sdc[state] += 1
                
                # update dictionary of pairs
                pairs_dict[p][0] += sdc['ZZ'] 
                pairs_dict[p][1] += sdc['Zz']
                pairs_dict[p][2] += sdc['zZ']
                pairs_dict[p][3] += sdc['zz']
    
    # save output
    out = open(args.outfile + 'CpG_pair_states.txt', 'wb')
    header = ['chr','cpg1','cpg2','ZZ','Zz','zZ','zz']
    out.write('\t'.join(header) + '\n')
    
    for i in pairs_dict:
        if sum(pairs_dict[i]) >= args.minCoverage:
            chrom, cpg1, cpg2 = i[0], i[1], i[2]
            ZZ, Zz, zZ, zz = pairs_dict[i][0], pairs_dict[i][1], pairs_dict[i][2], pairs_dict[i][3]
            out.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\n" % (chrom, cpg1, cpg2, ZZ, Zz, zZ, zz))
    
    
    out.close()

if __name__ == '__main__':
    main()
    
    
    
