#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
from Bio import SeqIO

def usage():
    test="name"
    message='''
python RedistributeBam.py --input .

Read in all the bam files in specified dir and redistributed into structured directory.
mv MSU7.Chr4.mPing.rep1_reads_10X_100_500.bam ./MSU7.Chr4.mPing/MSU7.Chr4.mPing.rep1_reads/
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
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.output:
        args.output = args.input

    #./MSU7.Chr4.mPing.rep1_reads_10X_100_500.bam
    s = re.compile(r'.*/((.*)\.rep\d+\_reads).*')
    bams = glob.glob('%s/*.bam' %(args.input))
    for bam in bams:
        m = s.search(bam)
        if m:
            outdir = '%s/%s/%s' %(args.output, m.groups(0)[1], m.groups(0)[0])
            clear  = '%s.clear.sh*' %(os.path.splitext(bam)[0])
            dupli  = '%s.dupli' %(os.path.splitext(bam)[0])
            mv     = 'mv %s %s' %(bam, outdir)
            mv2    = 'mv %s.bai %s' %(bam, outdir)
            clean  = 'rm -R %s %s' %(clear, dupli)
            os.system(mv)
            os.system(mv2)
            os.system(clean)
            #print bam, outdir
            #print clear, dupli

if __name__ == '__main__':
    main()

