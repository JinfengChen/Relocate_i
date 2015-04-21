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
python RunTEMP.py --input ../ReadMapping/MSU7.Chr4.mPing

Can not run TEMP jobs same time in one directory. Some files have same name in process?

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
        args.repeat = '/rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/mping_pogo.fa'

    TEMP = '/rhome/cjinfeng/software/tools/TEMP'
    BED  = '%s.RepeatMasker.out.bed' %(args.genome)
    Reference= args.genome
    Repeat   = args.repeat 
    lib_range= 0.1

    #../ReadMapping/MSU7.Chr4.mPing
    project = os.path.split(args.input)[1]
    if not os.path.exists(project):
        os.mkdir(project)
    #print project
    #MSU7.Chr4.mPing.rep1_reads_20X_100_500.bam
    s = re.compile(r'(.*_(\d+)).bam')
    reps = glob.glob('%s/*.*_reads' %(args.input))
    for rep in reps:
        #print rep
        rep_read = os.path.split(rep)[1]
        rep_read_dir = './%s/%s' %(project, rep_read)
        if not os.path.exists(rep_read_dir):
            os.mkdir(rep_read_dir)
        #print rep_read
        bams = glob.glob('%s/*.bam' %(rep))
        for bam in bams:
            m = s.search(bam)
            if m:
                read_pre = os.path.split(m.groups(0)[0])[1]
                lib_size = m.groups(0)[1]
                lib_sd   = float(lib_size) * float(lib_range)
                #print '%s\t%s\t%s\t%s' %(bam, lib_size, lib_size, lib_sd)
                #print read_pre
                #bash $TEMP/scripts/TEMP_Insertion.sh -i $bam -s $TEMP/scripts -r $TE_Fasta -t $TE_BED -m 3 -f 400 -c 8
                tempout = '%s.insertion.refined.bp.summary' %(read_pre)
                tempdir = '%s_TEMP' %(read_pre)
                cmd = 'bash %s/scripts/TEMP_Insertion.sh -i %s -s %s/scripts -r %s -t %s -m 3 -f %s -c 1' %(TEMP, bam, TEMP, Repeat, BED, lib_size)
                gff = 'python /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/TEMP/TEMP2GFF.py --insertion %s' %(tempout)
                rm  = 'rm %s.bam %s.bam.bai' %(read_pre, read_pre)
                mv  = 'mv %s.* %s' %(read_pre, tempdir)
                mv2 = 'mv %s %s' %(tempdir, rep_read_dir)
                #print '%s\n%s\n%s'%(cmd, gff, mv)
                print '%s/%s' %(rep_read_dir, tempdir)
                if not os.path.isfile(tempout) and not os.path.exists('%s/%s' %(rep_read_dir, tempdir)):
                    print '%s/%s' %(rep_read_dir, tempdir)
                    if not os.path.exists(tempdir):
                        os.mkdir(tempdir)
                    start_time = time.time()
                    os.system(cmd)
                    end_time   = time.time()
                    ofile = open('%s.run.log' %(read_pre), 'w')
                    print >> ofile, 'Run Time: %s seconds' %(end_time - start_time)
                    ofile.close()
                    os.system(gff)
                    os.system(rm)
                    os.system(mv)
                    os.system(mv2)
                

if __name__ == '__main__':
    main()

