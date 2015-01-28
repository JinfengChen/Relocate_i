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
python relocaTE_trim.py te_repeat.blatout 500_1.fq 10 0

Parse blat or bam file. Write TE matched reads and their pairs into files.
*.ContainingReads.fq: all reads have matches to TE
*.five_prime.fa/*.three_prime.fa: TE proportion of reads that matched to TE
*.potential_tandemInserts_containing_reads.list.txt: reads have two or more matches on TE and on same strand. Tandem insertions and Highly false positive.
    '''
    print message

def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(seq)
    for i in range(len(bases)):
        bases[i] = complement[bases[i]] if complement.has_key(bases[i]) else bases[i]
    return ''.join(bases)

def reverse_complement(seq):
    return complement(seq[::-1])

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


def parse_align_blat(infile, tandem):
    coord = defaultdict(lambda : defaultdict(lambda : str))

    ##align_file
    ofile = open(tandem, 'w')
    with open (infile, 'r') as filehd:
        for i in range(5):
            next(filehd)
        for line in filehd:
            line = line.rstrip()
            unit = re.split(r'\t',line)
            match    = int(unit[0])
            mismatch = int(unit[1])
            strand   = unit[8]
            qName    = unit[9]
            qLen     = int(unit[10])
            qStart   = int(unit[11])
            qEnd     = int(unit[12]) - 1
            tLen     = int(unit[14])
            tStart   = int(unit[15])
            #get all values into 1st base = 0 postion notation
            tEnd     = int(unit[16]) - 1 
            addRecord = 0

            if coord.has_key(qName):
                if strand == coord[qName]['strand']:
                    ##if there is are tandem insertions, theses reads
                    ##will call many false insertions events.
                    ##if on the same strand, the matches are not likely to overlap,
                    ##the two regions would be separated by mismatches or gaps if they were not
                    ##tandem insetions
                    print >> ofile, qName
                    qStart = qStart if qStart < int(coord[qName]['start']) else int(coord[qName]['start'])
                    qEnd   = qEnd if qEnd > int(coord[qName]['end']) else int(coord[qName]['end'])
                    tStart = tStart if tStart < int(coord[qName]['tStart']) else int(coord[qName]['tStart'])
                    tEnd   = tEnd if tEnd > int(coord[qName]['tEnd']) else int(coord[qName]['tEnd'])
                    match  = qLen if (int(coord[qName]['match']) + int(match)) > int(qLen) else int(coord[qName]['match']) + int(match)
                    addRecord = 1
                else:
                    ##keep the best match to TE
                    if int(coord[qName]['match']) >= match:
                        addRecord = 0
                    else:
                        addRecord = 1
            else:
                addRecord = 1

            if addRecord == 1:
                coord[qName]['match']    = match
                coord[qName]['len']      = qLen
                coord[qName]['start']    = qStart
                coord[qName]['end']      = qEnd
                coord[qName]['tLen']     = tLen
                coord[qName]['mismatch'] = mismatch
                coord[qName]['strand']   = strand
                coord[qName]['tStart']   = tStart
                coord[qName]['tEnd']     = tEnd
    ofile.close()
    return coord
 

def main():
    if not len(sys.argv) == 5:
        usage()
        sys.exit(2)

    align_file = sys.argv[1]
    fq_file1   = sys.argv[2]
    len_cutoff = int(sys.argv[3])
    mismatch_allowance = float(sys.argv[4])
    align_type = 'bam' if os.path.splitext(align_file)[1].replace(r'.', '') == 'bam' else 'blat'
  
    #set output directories and files
    align_path = os.path.abspath(os.path.split(align_file)[0])
    file_name  = os.path.split(align_file)[1]
    te_path    = re.split(r'/', align_path)
    del te_path[-1]
    out_fq_path= '%s/te_containing_fq' %('/'.join(te_path))
    out_fa_path= '%s/te_only_read_portions_fa' %('/'.join(te_path))
    tandem_file= '%s/%s.potential_tandemInserts_containing_reads.list.txt' %(out_fq_path, file_name)
    
    #parse align
    coord = defaultdict(lambda : defaultdict(lambda : str))
    if align_type == 'blat':
        coord = parse_align_blat(align_file, tandem_file)
    
    #outfiles
    TE = 'unspecified'
    FA = 'unspecified' 
    s = re.compile(r'(\S+)\.te_(\S+)\.blatout')
    m = s.search(file_name)
    if m:
        FA = m.groups(0)[0]
        TE = m.groups(0)[1]
    outfq  = 0
    outte5 = 0
    outte3 = 0
    
    if os.path.isfile('%s/%s.te_%s.ContainingReads.fq' %(out_fq_path, FA, TE)) and int(os.path.getsize('%s/%s.te_%s.ContainingReads.fq' %(out_fq_path, FA, TE))) > 0:
        outfq = 1        
    if os.path.isfile('%s/%s.te_%s.five_prime.fa' %(out_fa_path, FA, TE)) and int(os.path.getsize('%s/%s.te_%s.five_prime.fa' %(out_fa_path, FA, TE))) > 0:
        outte5 = 1
    if os.path.isfile('%s/%s.te_%s.three_prime.fa' %(out_fa_path, FA, TE)) and int(os.path.getsize('%s/%s.te_%s.five_prime.fa' %(out_fa_path, FA, TE))) > 0:
        outte3 = 1
    #print 'check1: %s\t%s\t%s' %(str(outfq), str(outte5), str(outte3))
    if outfq + outte5 + outte3 == 3:
        exit(0)
    else:
        if os.path.isfile('%s/%s.te_%s.ContainingReads.fq' %(out_fq_path, FA, TE)):
            cmd1 = 'rm %s/%s.te_%s.ContainingReads.fq' %(out_fq_path, FA, TE)
            os.system(cmd1)
        if os.path.isfile('%s/%s.te_%s.five_prime.fa' %(out_fa_path, FA, TE)):
            cmd2 = 'rm %s/%s.te_%s.five_prime.fa' %(out_fa_path, FA, TE)
            os.system(cmd2)
        if os.path.isfile('%s/%s.te_%s.three_prime.fa' %(out_fa_path, FA, TE)):
            cmd3 = 'rm %s/%s.te_%s.three_prime.fa' %(out_fa_path, FA, TE)
            os.system(cmd3)
        ofile_fq = open('%s/%s.te_%s.ContainingReads.fq' %(out_fq_path, FA, TE), 'w')
        ofile_te5= open('%s/%s.te_%s.five_prime.fa' %(out_fa_path, FA, TE), 'w')
        ofile_te3= open('%s/%s.te_%s.three_prime.fa' %(out_fa_path, FA, TE), 'w')
        ##go thru each fq record in the fq files. if the name of the seq is in the blat file
        ##trim the seq
        with open (fq_file1, 'r') as filehd:
            for line in filehd:
                line = line.rstrip()
                header = line[1:]
                header = re.split(r'\s+', header)[0]
                seq    = filehd.next().rstrip()
                qualh  = filehd.next().rstrip()
                qual   = filehd.next().rstrip()
                
                if coord.has_key(header):
                    start    = int(coord[header]['start'])
                    length   = int(coord[header]['len'])
                    end      = int(coord[header]['end'])
                    tStart   = int(coord[header]['tStart'])
                    tEnd     = int(coord[header]['tEnd'])
                    tLen     = int(coord[header]['tLen'])
                    mismatch = int(coord[header]['mismatch'])
                    match    = int(coord[header]['match'])
                    strand   = coord[header]['strand']
                   
                    #want to cut and keep anything not matching to database TE
                    trimmed_seq  = ''
                    trimmed_qual = ''
                    #print header, seq, qual
                    #print 'check2: %s\t%s\t%s'  %(str(tStart), str((length - (match + mismatch))), str((mismatch/(match + mismatch))))
                    ##query read overlaps 5' end of database TE & trimmed seq > cutoff
                    if tStart == 0 and (length - (match + mismatch)) > len_cutoff and (float(mismatch)/(float(match) + float(mismatch))) <= mismatch_allowance:
                        tS = int(tStart) + 1
                        tE = int(tEnd) + 1
                        qS = int(start) + 1
                        qE = int(end) + 1
                        ## te_subseq = portion of the seq that matches to TE
                        ## trimmed_seq = portion of the seq that does not match to T
                        #end1 = end + 1
                        te_subseq = seq[start:end+1]
                        if strand == '-':
                            te_subseq = reverse_complement(te_subseq)
                            trimmed_seq = seq[end:]
                            trimmed_qual= qual[end:]
                            seq_id   = header
                            seq_desc = ''
                            seq_id   = '%s:start:5' %(seq_id)
                            header = '%s%s' %(seq_id, seq_desc)
                        else:
                            trimmed_seq = seq[0:start]
                            trimmed_qual= qual[0:start]
                            seq_id   = header
                            seq_desc = ''
                            seq_id   = '%s:end:5' %(seq_id)
                            header = '%s%s' %(seq_id, seq_desc)
                        #print 'trimmed: %s %s' %(trimmed_seq, str(end))
                        if len(trimmed_seq) >= len_cutoff:
                            print >> ofile_te5, '>%s %s..%s matches %s:%s..%s mismatches:%s\n%s' %(header, qS, qE, TE, tS, tE, mismatch, te_subseq)
                    elif tEnd == tLen - 1 and (length - (match + mismatch)) > len_cutoff and (float(mismatch)/(float(match) + float(mismatch))) <= mismatch_allowance:
                        tS = int(tStart) + 1
                        tE = int(tEnd) + 1
                        qS = int(start) + 1
                        qE = int(end) + 1
                        te_subseq = seq[start:end+1]
                        if strand == '-':
                            te_subseq = reverse_complement(te_subseq)
                            trimmed_seq = seq[0:start]
                            trimmed_qual= qual[0:start]
                            seq_id   = header
                            seq_desc = ''
                            seq_id   = '%s:end:3' %(seq_id)
                            header = '%s%s' %(seq_id, seq_desc)
                        else:
                            trimmed_seq = seq[end:]
                            trimmed_qual= qual[end:]
                            seq_id   = header
                            seq_desc = ''
                            seq_id   = '%s:start:3' %(seq_id)
                            header = '%s%s' %(seq_id, seq_desc)
                        if len(trimmed_seq) >= len_cutoff:
                            print >> ofile_te3, '>%s %s..%s matches %s:%s..%s mismatches:%s\n%s' %(header, qS, qE, TE, tS, tE, mismatch, te_subseq)
                    ##trimmed reads
                    if len(trimmed_seq) >= len_cutoff:
                        print '@%s\n%s\n%s\n%s' %(header, trimmed_seq, qualh, trimmed_qual)
                    ##any read that was in the blat file is written here
                    print >> ofile_fq, '@%s\n%s\n%s\n%s' %(header, seq, qualh, qual)

if __name__ == '__main__':
    main()

