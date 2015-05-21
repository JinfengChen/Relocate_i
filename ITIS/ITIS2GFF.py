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
python ITIS2GFF.py --insertion all.filterd.bed

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

'''
Chr     Start   End     TransposonName  TransposonDirection     Class   VariantSupport  Frequency       Junction1       Junction1Support        Junction2       Junction2Support        5'_Support      3'_Support
chr2L   2003846 2003873 FBgn0003122_pogo        sense   1p1     39      1.0000  2003871 2       2003871 5       18      14

Chr1    HEG4_2  transposable_element_insertion_site     23783719        23783723        .       -       .       ID=Copia2_LTR.te_insertion_site.Chr1.23783723;Note=Non-reference
'''
def insertion2gff(infile, title):
    data = defaultdict(str)
    gff = infile + '.gff'
    ofile = open(gff, 'w')
    count = 0
    with open (infile, 'r') as filehd:
        filehd.readline()
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                count += 1 
                unit = re.split(r'\t',line)
                chrs = unit[0]
                start= unit[1]
                end  = unit[2]
                name = unit[3]
                ID   = 'TE_Insertion_%s' %(count) 
                strand = '+' if unit[4] is 'sense' else '-'
                if int(unit[9]) > 0 and int(unit[11]) > 0:
                    gffline = '%s\t%s\tTEMP\t%s\t%s\t.\t%s\t.\tID=%s;Note=Non-reference;Name=%s;Class=%s;VariantSupport=%s;Frequency=%s;Junction1Support=%s;Junction2Support=%s;Five_Support=%s;Three_Support=%s;' %(chrs, title, unit[8], unit[10], strand, ID, name, unit[5], unit[6], unit[7], unit[9], unit[11], unit[12], unit[13])
                    print >> ofile, gffline
                else:
                    gffline = '%s\t%s\tTEMP\t%s\t%s\t.\t%s\t.\tID=%s;Note=Non-reference;Name=%s;Class=%s;VariantSupport=%s;Frequency=%s;Junction1Support=%s;Junction2Support=%s;Five_Support=%s;Three_Support=%s;' %(chrs, title, start, end, strand, ID, name, unit[5], unit[6], unit[7], unit[9], unit[11], unit[12], unit[13])
                    print >> ofile, gffline
    return data
    ofile.close()

'''
Chr1    HEG4_2  transposable_element_insertion_site     23783719        23783723        .       -       .       ID=Copia2_LTR.te_insertion_site.Chr1.23783723;Note=Non-reference
'''
def absence2gff(infile):
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
    parser.add_argument('-i', '--insertion')
    parser.add_argument('-a', '--absence')
    parser.add_argument('-p', '--project')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
 
    if (args.project is None):
        args.project = 'HEG4'

    try:
        len(args.insertion) > 0
    except:
        usage()
        sys.exit(2)

    result_bed = args.insertion
    result_gff = re.sub(r'bed', r'gff', result_bed)
    print result_bed
    print result_gff
    bed2gff    = 'awk \'{print $1"\\tHEG4\\tITIS\\t"$2"\\t"$3"\\t.\\t"$6"\\t.\\tID="$1"_"$2";"$4}\' %s > %s' %(result_bed, result_gff) 
    os.system(bed2gff)

if __name__ == '__main__':
    main()

