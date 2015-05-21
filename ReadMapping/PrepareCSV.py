#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
import time
from Bio import SeqIO

def usage():
    test="name"
    message='''
python PrepareCSV.py --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/MSU7.Chr4.mPing


    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


def readtable(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[0]] = unit[1]
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-g', '--genome')
    parser.add_argument('-r', '--repeat')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.genome:
        args.genome = '/rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/MSU7.Chr4.fa'
  
    if not args.repeat:
        args.repeat = '/rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/mping.fa'

    #-t ../input/mping_UNK.fa -g /rhome/cjinfeng/HEG4_cjinfeng/seqlib/MSU_r7.fa -d ../input/FC52_7 -e HEG4 -o mPing_HEG4_UNK -r 1 -p 1 -a 1   
    Reference= args.genome
    Repeat   = args.repeat 

    project = os.path.split(args.input)[1]
    if not os.path.exists(project):
        os.mkdir(project)
    #print project
    s = re.compile(r'(.*_(\d+))_1.fq')
    reps = glob.glob('%s/*.*_reads' %(args.input))
    for rep in reps:
        #print rep
        rep_read = os.path.split(rep)[1]
        rep_read_dir = './%s/%s' %(project, rep_read)
        if not os.path.exists(rep_read_dir):
            os.mkdir(rep_read_dir)
        #print rep_read
        reads = glob.glob('%s/*.*_1.fq' %(rep))
        for read in reads:
            m = s.search(read)
            if m:
                read_pre = m.groups(0)[0]
                lib_size = m.groups(0)[1]
                lib_sd   = float(lib_size) * 0.2
                read_idx = '%s_?.fq' %(read_pre)
                #MSU7.Chr4.mPing.rep1_reads_2X_100_500.bam
                bam = '%s/%s.bam' %(rep_read_dir, os.path.split(read_pre)[1])
                if not os.path.isfile(bam):
                    print '%s\t%s\t%s\t%s' %(read_idx, lib_size, lib_size, lib_sd)
                

if __name__ == '__main__':
    main()

