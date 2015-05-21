#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO

def usage():
    test="name"
    message='''
python Move.py --input MSU7.Chr4.mPing 

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

    #MSU7.Chr4.mPing.rep3_reads_9X_100_500_ITIS
    for i in [1,2,3]:
        for j in [1,2,3,4,5,6,7,8,9,10,15,20,40]:
            itis = '%s.rep%s_reads_%sX_100_500_ITIS*' %(args.input, i, j)
            out  = '%s/%s.rep%s_reads' %(args.input, args.input, i)
            #print itis, out
            os.system('mv %s %s' %(itis, out))

if __name__ == '__main__':
    main()

