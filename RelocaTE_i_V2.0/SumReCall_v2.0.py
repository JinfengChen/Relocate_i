#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
import glob
import time
from Bio import SeqIO

def usage():
    test="name"
    message='''
python SumReCall.py --call TEMP --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/MSU7.Chr4.mPing
Summary RelocatTE call in current direcory using simulation in input directory
--call: RelocaTE, TEMP or other
--input: dir of simulation, where we can find insertion simulated gff file "MSU7.Chr4.mPing.rep1.gff"

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


def parse_ref_gff(infile):
    data = defaultdict(int)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'#'):
                #print line 
                unit = re.split(r'\t',line)
                te   = '%s_%s_%s' %(unit[0], unit[3], unit[4])
                data[te] = 1
    return data


#Chr3    MSU7    Gaijin  472029  472030  .       -       .       ID=Gaijin.Chr3.472029;TSD=CTT;
def parse_overlap_reloacte(infile, ref_te, call_te):
    #Total number of call, true call, call with breakpoint near TSD, No call, False call
    data = defaultdict(lambda : int)
    dupli= defaultdict(lambda : int())
    true = 0
    tsd  = 0
    r = re.compile(r'TSD=(\w+);')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                #print line
                if dupli.has_key('%s_%s' %(unit[12], unit[13])):
                    continue
                else:
                    dupli['%s_%s' %(unit[12], unit[13])] == 1
                #print 'pass'
                tsd_wd = r.search(unit[8]).groups(0)[0] if r.search(unit[8]) else 'NA'
                tsd_s  = int(unit[3]) - len(tsd_wd) + 1
                tsd_e  = int(unit[3])
                if int(unit[12]) == int(tsd_s) and int(unit[13]) == int(tsd_e):
                    tsd += 1
                #pos  = map(int, [unit[3], unit[4], unit[12], unit[13]])
                #dist_min = min([abs(pos[2]-pos[0]), abs(pos[2]-pos[1]), abs(pos[3]-pos[0]), abs(pos[3]-pos[1])])
                #dist_max = max([abs(pos[2]-pos[0]), abs(pos[2]-pos[1]), abs(pos[3]-pos[0]), abs(pos[3]-pos[1])])
                true += 1
                #if dist_min <= 10 and dist_min <= 10:
                    #tsd += 1
    call = len(call_te.keys())
    ref  = len(ref_te.keys())
    data = [call, true, tsd, ref-true, call-true]
    return data

def parse_overlap_reloacte1(infile, ref_te, call_te):
    #Total number of call, true call, call with breakpoint near TSD, No call, False call
    data = defaultdict(lambda : int)
    dupli= defaultdict(lambda : int())
    true = 0
    tsd  = 0
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                #print line
                if dupli.has_key('%s_%s' %(unit[12], unit[13])):
                    continue
                else:
                    dupli['%s_%s' %(unit[12], unit[13])] == 1
                #print 'pass'
                pos  = map(int, [unit[3], unit[4], unit[12], unit[13]])
                dist_min = min([abs(pos[2]-pos[0]), abs(pos[2]-pos[1]), abs(pos[3]-pos[0]), abs(pos[3]-pos[0])])
                dist_max = max([abs(pos[2]-pos[0]), abs(pos[2]-pos[1]), abs(pos[3]-pos[0]), abs(pos[3]-pos[0])])
                true += 1
                if dist_min <= 5 and dist_max <= 5:
                    tsd += 1
    call = len(call_te.keys())
    ref  = len(ref_te.keys())
    data = [call, true, tsd, ref-true, call-true]
    return data




def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-c', '--call')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    data = defaultdict(lambda : defaultdict(lambda : list)) 
    project = os.path.split(args.input)[1]
    reps = glob.glob('%s/*.*_reads' %(args.input))
    for rep in reps:
        rep_read = os.path.split(rep)[1]
        rep_read_dir = './%s/%s' %(project, rep_read)
        #print rep_read
        #print rep_read_dir
        calls = glob.glob('%s/*.*_%s' %(rep_read_dir, args.call))
        for call in calls: 
            #MSU7.Chr4.mPing.rep1_reads_10X_100_500_RelocaTE
            unit = re.split(r'_', os.path.split(call)[1])
            head = re.split(r'\.', unit[0])
            genome = head[0]
            chro   = head[1]
            repeat = head[2]
            replicate = head[3]
            depth     = unit[2].replace('X', '')
            readlen   = unit[3]
            libsize   = unit[4]
            #print '%s\t%s\t%s\t%s' %(genome, chro, repeat, replicate)
            #
            ref_gff = '%s/%s.gff' %(args.input, unit[0])
            ref_te  = parse_ref_gff(ref_gff)
            eval_list = []
            if args.call == 'RelocaTE':
                if 0: #use only inseretions with two breakpoint
                    #MSU7.Chr4.mPing.rep1_reads_10X_100_500_RelocaTE/mPing/results/HEG4.mPing.all_inserts.gff
                    result_gff = '%s/%s/results/HEG4.%s.all_inserts.gff' %(call, repeat, repeat)
                    #print result_gff
                    non_ref_gff= 'grep "Non-reference" %s > temp.relocate.gff' %(result_gff)
                    win_overlap= 'bedtools window -w 100 -a %s -b temp.relocate.gff > temp.relocate.overlap' %(ref_gff)
                    #print win_overlap
                    os.system(non_ref_gff)
                    os.system(win_overlap)
                    call_te = parse_ref_gff('temp.relocate.gff')
                    eval_list = parse_overlap_reloacte('temp.relocate.overlap', ref_te, call_te)
                    index = '_'.join([depth, genome, chro, repeat, readlen, libsize, args.call])
                    data[index][replicate] = eval_list
                    #print index
                    #print '%s\t%s\t%s\t%s\t%s' %(eval_list[0], eval_list[1], eval_list[2], eval_list[3], eval_list[4])
                else:
                    non_ref_txt = '%s/%s/results/HEG4.%s.all_nonref.txt' %(call, repeat, repeat)
                    non_ref_gff = '%s/%s/results/HEG4.%s.all_nonref.gff' %(call, repeat, repeat)
                    convertgff  = 'perl /rhome/cjinfeng/software/bin/relocate2gff.pl --mping %s' %(non_ref_txt)
                    os.system(convertgff)
                    os.system('cp %s temp.relocate.gff' %(non_ref_gff))
                    win_overlap= 'bedtools window -w 100 -a %s -b temp.relocate.gff > temp.relocate.overlap' %(ref_gff)
                    os.system(win_overlap)
                    call_te = parse_ref_gff('temp.relocate.gff')
                    eval_list = parse_overlap_reloacte('temp.relocate.overlap', ref_te, call_te)
                    index = '_'.join([depth, genome, chro, repeat, readlen, libsize, args.call])
                    data[index][replicate] = eval_list
                    #print index
                    #print '%s\t%s\t%s\t%s\t%s' %(eval_list[0], eval_list[1], eval_list[2], eval_list[3], eval_list[4])
            elif args.call == 'RelocaTEi':
                non_ref_txt = '%s/repeat/results/%s.repeat.all_nonref_insert.txt' %(call, chro)
                non_ref_gff = '%s/repeat/results/%s.repeat.all_nonref_insert.gff' %(call, chro)
                non_ref_gff_char = '%s/repeat/results/ALL.all_nonref_insert.characTErized.gff' %(call)
                non_ref_gff_all  = '%s/repeat/results/ALL.all_nonref_insert.gff' %(call)
                convertgff  = 'perl /rhome/cjinfeng/software/bin/relocate2gff.pl --mping %s' %(non_ref_txt)
                if not os.path.isfile(non_ref_gff):
                    os.system(convertgff)
                #os.system('cp %s temp.relocate.gff' %(non_ref_gff))
                os.system('cp %s temp.relocate.gff' %(non_ref_gff_all))
                win_overlap= 'bedtools window -w 100 -a %s -b temp.relocate.gff > temp.relocate.overlap' %(ref_gff)
                os.system(win_overlap)
                call_te = parse_ref_gff('temp.relocate.gff')
                eval_list = parse_overlap_reloacte('temp.relocate.overlap', ref_te, call_te)
                index = '_'.join([depth, genome, chro, repeat, readlen, libsize, args.call])
                data[index][replicate] = eval_list
            elif args.call == 'TEMP':
                #MSU7.Chr4.mPing.rep1_reads_10X_100_500.insertion.refined.bp.summary.gff
                gff_pre = '_'.join([unit[0], unit[1], unit[2], unit[3], unit[4]])
                summary    = '%s/%s.insertion.refined.bp.summary' %(call, gff_pre)
                result_gff = '%s/%s.insertion.refined.bp.summary.gff' %(call, gff_pre)
                os.system('python /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/TEMP/TEMP2GFF.py --insertion %s' %(summary))
                win_overlap= 'bedtools window -w 100 -a %s -b %s > temp.temp.overlap' %(ref_gff, result_gff)
                os.system(win_overlap)
                call_te = parse_ref_gff(result_gff)
                eval_list = parse_overlap_reloacte('temp.temp.overlap', ref_te, call_te)
                index = '_'.join([depth, genome, chro, repeat, readlen, libsize, args.call])
                data[index][replicate] = eval_list

    ofile = open('%s.%s.summary' %(project, args.call), 'a')
    print >> ofile, 'Depth\tGenome\tChr\tElement\tReadLength\tLibrarySize\tTools\t#Call\tSTD\t#TRUE\tSTD\t#TSD\tSTD\t#Fail_Call\tSTD\t#False_Call\tSTD'
    for idx in sorted(data.keys()):
        eval1 = data[idx]['rep1']
        eval2 = data[idx]['rep2']
        eval3 = data[idx]['rep3']
        #print '%s\t%s\t%s\t%s\t%s' %(eval2[0], eval2[1], eval2[2], eval2[3], eval2[4])
        call_avg = mean([eval1[0], eval2[0], eval3[0]])
        call_sd  = std([eval1[0], eval2[0], eval3[0]])
        true_avg = mean([eval1[1], eval2[1], eval3[1]])
        true_sd  = std([eval1[1], eval2[1], eval3[1]])
        tsd_avg  = mean([eval1[2], eval2[2], eval3[2]])
        tsd_sd   = std([eval1[2], eval2[2], eval3[2]])
        Non_avg  = mean([eval1[3], eval2[3], eval3[3]])
        Non_sd   = std([eval1[3], eval2[3], eval3[3]])
        False_avg  = mean([eval1[4], eval2[4], eval3[4]])
        False_sd   = std([eval1[4], eval2[4], eval3[4]])
        head = '\t'.join(re.split(r'_', idx))
        eval0 = map(str, [call_avg, call_sd, true_avg, true_sd, tsd_avg, tsd_sd, Non_avg, Non_sd, False_avg, False_sd])
        print >> ofile, '%s\t%s' %(head, '\t'.join(eval0))
    ofile.close()
  
if __name__ == '__main__':
    main()

