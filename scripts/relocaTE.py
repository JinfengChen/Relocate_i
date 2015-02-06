#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO
import glob

def usage():
    test="name"
    message='''
python CircosConf.py --input circos.config --output pipe.conf

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

def createdir(dirname):
    if not os.path.exists(dirname):
        os.mkdir(dirname)

def writefile(outfile, lines):
    ofile = open(outfile, 'w')
    print >> ofile, lines
    ofile.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bam')
    parser.add_argument('-t', '--te_fasta')
    parser.add_argument('-d', '--fq_dir')
    parser.add_argument('-g', '--genome_fasta')
    parser.add_argument('-r', '--reference_ins')
    parser.add_argument('-f', '--fastmap')
    parser.add_argument('-o', '--outdir')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    
    try:
        os.path.isfile(args.bam) and os.path.isfile(args.te_fasta) and os.path.isfile(args.genome_fasta)
    except:
        pass
    else:
        try:
            os.path.exists(args.fq_dir) and os.path.isfile(args.te_fasta) and os.path.isfile(args.genome_fasta)
        except:
            usage()
            sys.exit(2)

    #absolute path for all files, directory and scripts
    RelocaTE_bin = os.path.split(os.path.abspath(__file__))[0]
    reference    = os.path.abspath(args.genome_fasta)
    te_fasta     = os.path.abspath(args.te_fasta)
    mode         = 'bam'
    bam          = ''
    fastq_dir    = ''
    
    try: 
        if os.path.isfile(args.bam):
            bam  = os.path.abspath(args.bam)
    except:
        try:
            if os.path.abspath(args.fq_dir):
                fastq_dir = os.path.abspath(args.fq_dir)
                mode      = 'fastq'
        except:
            usage()
            exit(2)

    print fastq_dir
    print mode

    #Prepare directory and script
    if args.outdir is None:
        args.outdir = '%s/RelocaTE_output' %(os.getcwd())
        createdir(args.outdir)
    else:
        args.outdir = os.path.abspath(args.outdir)
        createdir(args.outdir)

    writefile('%s/regex.txt' %(args.outdir), '_1\t_2\t.unPaired\t...')
    createdir('%s/shellscripts' %(args.outdir))
    createdir('%s/repeat' %(args.outdir))
    createdir('%s/repeat/blat_output' %(args.outdir))
    createdir('%s/repeat/flanking_seq' %(args.outdir))
    createdir('%s/repeat/te_containing_fq' %(args.outdir))
    createdir('%s/repeat/te_only_read_portions_fa' %(args.outdir))

    shells = []
    #step0 existing TE blat
    reference_ins_flag = 'NONE'
    if args.reference_ins is None or args.reference_ins == '0':
        step0_file = '%s/shellscripts/step_0_do_not_call_reference_insertions' %(args.outdir)
        writefile(step0_file, '')
    elif os.path.isfile(args.reference_ins):
        step0_file = '%s/shellscripts/step_0_te_annotation_provided' %(args.outdir)
        reference_ins_flag = args.reference_ins
        writefile(step0_file, '')
    elif args.reference_ins == '1':
        createdir('%s/shellscripts/step_0' %(args.outdir))
        step0_file = '%s/shellscripts/step_0/step_0.existingTE_blat.sh' %(args.outdir)
        shells.append('sh %s' %(step0_file))
        existingTE_blat = 'blat %s %s %s/existingTE.blatout 1> %s/existingTE.blat.stdout' %(reference, te_fasta, args.outdir, args.outdir)
        reference_ins_flag = '%s/existingTE.blatout' %(args.outdir)
        writefile(step0_file, existingTE_blat)

    #step1 format reference genome

    #step2 fastq to fasta
    fastas = defaultdict(lambda : str)
    if mode == 'fastq':
        fastqs = glob.glob('%s/*.f*q' %(fastq_dir))
        step2_flag = 0
        step2_count= 0
        for fq in fastqs:
            fa    = '%s.fa' %(os.path.splitext(fq)[0])
            fastas[fa] = fq
            if not os.path.isfile(fa):
                createdir('%s/shellscripts/step_2' %(args.outdir))
                fq2fa = '%s/relocaTE_fq2fa.pl %s %s' %(RelocaTE_bin, fq, fa)
                step2_file = '%s/shellscripts/step_2/%s.fq2fq.sh' %(args.outdir, step2_count)
                shells.append('sh %s' %(step2_file))
                writefile(step2_file, fq2fa)
                step2_flag == 1
                step2_count += 1
        if step2_flag == 0:
            step2_file = '%s/shellscripts/step_2_not_needed_fq_already_converted_2_fa' %(args.outdir)
            writefile(step2_file, '')        
    elif mode == 'bam':
        print 'Add module of obtaining reads from bam then prepare as fa files'
    
    #step3 blat fasta to repeat
    step3_count = 0
    for fa in sorted(fastas.keys()):
        createdir('%s/shellscripts/step_3' %(args.outdir))
        fq      = fastas[fa]
        #fq      = '%s.fq' %(os.path.splitext(fa)[0]) if os.path.isfile('%s.fq' %(os.path.splitext(fa)[0])) else '%s.fastq' %(os.path.splitext(fa)[0])
        fa_prefix = os.path.split(os.path.splitext(fa)[0])[1]
        blatout = '%s/repeat/blat_output/%s.te_repeat.blatout' %(args.outdir, fa_prefix)
        blatstd = '%s/repeat/blat_output/blat.out' %(args.outdir)
        blat = 'blat -minScore=10 -tileSize=7 %s %s %s 1>> %s' %(te_fasta, fa, blatout, blatstd)
        if args.fastmap is not None:
            blat = 'blat -minScore=10 -tileSize=7 -fastMap %s %s %s 1>> %s' %(te_fasta, fa, blatout, blatstd)
        flank= '%s/repeat/flanking_seq/%s.te_repeat.flankingReads.fq' %(args.outdir, fa_prefix)
        trim = 'perl %s/relocaTE_trim.py %s %s 10 0 > %s' %(RelocaTE_bin, blatout, fq, flank)
        step3_file = '%s/shellscripts/step_3/%s.te_repeat.blat.sh' %(args.outdir, step3_count)
        shells.append('sh %s' %(step3_file))
        step3_cmds = '%s\n%s' %(blat, trim)
        writefile(step3_file, step3_cmds)
        step3_count += 1
            

    #step4 align TE trimed reads to genome
    ref = os.path.split(os.path.splitext(reference)[0])[1]
    createdir('%s/shellscripts/step_4' %(args.outdir))
    step4_file= '%s/shellscripts/step_4/step_4.%s.repeat.align.sh' %(args.outdir, ref)
    shells.append('sh %s' %(step4_file))
    step4_cmd = 'perl %s/relocaTE_align.py %s %s/repeat %s %s %s/regex.txt repeat not.given 0' %(RelocaTE_bin, RelocaTE_bin, args.outdir, reference, fastq_dir, args.outdir)
    writefile(step4_file, step4_cmd)
    
    #step5 find insertions
    ids = fasta_id(reference)
    createdir('%s/shellscripts/step_5' %(args.outdir))
    step5_count = 0
    for chrs in ids:
        step5_cmd = 'perl %s/relocaTE_insertionFinder.py %s/repeat/bwa_aln/%s.repeat.bwa.sorted.bam %s %s repeat %s/regex.txt not.give 100 %s 0 0' %(RelocaTE_bin, args.outdir, ref, chrs, reference, args.outdir, reference_ins_flag)
        step5_file= '%s/shellscripts/step_5/%s.repeat.findSites.sh' %(args.outdir, step5_count)
        shells.append('sh %s' %(step5_file))
        writefile(step5_file, step5_cmd)
        step5_count +=1
    
    
    #write script
    writefile('%s/run_these_jobs.sh' %(args.outdir), '\n'.join(shells))

if __name__ == '__main__':
    main()

