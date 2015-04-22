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
python RunRelocaTE.py --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/MSU7.Chr4.mPing

Can be run at same time for multi-project
    '''
    print message


def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --maxjob 20 --lines %s --interval 120 --resource nodes=1:ppn=1,walltime=100:00:00,mem=2G --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)


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
    #RelocaTE = 'python /rhome/cjinfeng/software/tools/RelocaTE_1.0.3_i/RelocaTE/scripts/relocaTE.py'
    RelocaTE = 'python /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/scripts_dev/relocaTE.py'
    Reference= args.genome
    Repeat   = args.repeat
    project = os.path.split(args.input)[1]
    if not os.path.exists(project):
        os.mkdir(project)
    #print project
    s = re.compile(r'(.*)_1.fq')
    reps = glob.glob('%s/*.*_reads' %(args.input))
    ofile = open('%s.run.sh' %(project), 'w')
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
                read_dir = './%s/%s' %(rep_read_dir, os.path.split(read_pre)[1])
                if not os.path.exists(read_dir):
                    os.mkdir(read_dir)
                read1 = '%s_1.fq' %(read_pre)
                read2 = '%s_2.fq' %(read_pre)
                ln1   = 'ln -s %s %s' %(read1, read_dir)
                ln2   = 'ln -s %s %s' %(read2, read_dir)
                if not os.path.exists('%s_1.fq' %('./%s/%s_1.fq' %(read_dir, os.path.split(read_pre)[1]))):
                    os.system(ln1)
                    os.system(ln2)
                outdir = '%s_RelocaTEi' %(os.path.split(read_pre)[1])
                #existingTE = os.path.dirname(read_pre)
                #existingTE = re.sub(r'_reads', '.RepeatMasker.out', existingTE)
                existingTE  = '%s.RepeatMasker.out' %(Reference)
                #print read_pre, existingTE
                # relocate will not run if there is result exists
                if not os.path.exists('./%s/%s' %(rep_read_dir, outdir)):
                    relocaTE = '%s --te_fasta %s --genome_fasta %s --fq_dir %s --outdir %s --reference_ins %s' %(RelocaTE, Repeat, Reference, os.path.abspath(read_dir), os.path.abspath(outdir), existingTE)
                    shell    = 'bash %s/run_these_jobs.sh > %s/run.log 2>%s/run.log2' %(os.path.abspath(outdir), os.path.abspath(outdir), os.path.abspath(outdir))
                    mv       = 'mv %s %s' %(os.path.abspath(outdir), os.path.abspath(rep_read_dir))
                    #os.system(relocaTE)
                    #start_time = time.time()
                    #os.system(shell)
                    #end_time   = time.time()
                    #ofile = open('./%s/run.log' %(outdir), 'a')
                    #print >> ofile, 'Run Time: %s seconds' %(end_time - start_time)
                    #ofile.close()
                    #os.system(mv)
                    #print relocaTE
                    print >> ofile, relocaTE
                    print >> ofile, shell
                    print >> ofile, mv
    ofile.close()
    runjob('%s.run.sh' %(project), 9)

if __name__ == '__main__':
    main()

