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
python RunITIS.py --input ../ReadMapping/MSU7.Chr4.mPing

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
        args.genome = '/rhome/cjinfeng/BigData/00.RD/RelocaTE_i/ITIS/reference/MSU7.Chr4.fa'
  
    if not args.repeat:
        args.repeat = '/rhome/cjinfeng/BigData/00.RD/RelocaTE_i/ITIS/reference/mping.fa'

    ITIS = '/rhome/cjinfeng/software/tools/SVcaller/ITIS_v1'
    Reference= args.genome
    Repeat   = args.repeat 
    lib_range= 0.1

    header = '''#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=2gb
#PBS -l walltime=100:00:00
#PBS -V

cd $PBS_O_WORKDIR

module load bwa/0.7.9
module load samtools/0.1.19
module load ncbi-blast/2.2.30+

'''

    #../ReadMapping/MSU7.Chr4.mPing
    project = os.path.split(args.input)[1]
    if not os.path.exists(project):
        os.mkdir(project)
    #print project
    s = re.compile(r'(.*)_1.fq')
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
                read1 = '%s_1.fq' %(read_pre)
                read2 = '%s_2.fq' %(read_pre)
                itis_run = 'perl %s/itis.pl -g %s -t %s -l %s -N %s -1 %s -2 %s -e Y -T %s_ITIS' %(ITIS, Reference, Repeat, 500, os.path.split(read_pre)[1], read1, read2, os.path.split(read_pre)[1])
                #print itis_run 
                outdir = '%s_ITIS' %(read_pre)
                tmpdir = os.path.split(outdir)[-1]
                #print rep_read_dir, tmpdir
                if not os.path.exists(tmpdir) and not os.path.exists('%s/%s' %(rep_read_dir, tmpdir)):
                    ofile = open ('%s.sh' %(tmpdir), 'w')
                    print >> ofile, header
                    print >> ofile, itis_run
                    print >> ofile, '\n\necho "Done"'
                    ofile.close()
                    os.system('qsub %s.sh' %(tmpdir))
                #    print '%s/%s' %(rep_read_dir, tempdir)
                #    if not os.path.exists(tempdir):
                #        os.mkdir(tempdir)
                #    start_time = time.time()
                #    os.system(cmd)
                #    end_time   = time.time()
                #    ofile = open('%s.run.log' %(read_pre), 'w')
                #    print >> ofile, 'Run Time: %s seconds' %(end_time - start_time)
                #    ofile.close()
                #    os.system(gff)
                #    os.system(rm)
                #    os.system(mv)
                #    os.system(mv2)
                

if __name__ == '__main__':
    main()

