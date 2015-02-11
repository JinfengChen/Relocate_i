#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from collections import OrderedDict
import re
import os
import argparse
import pysam
import glob

def usage():
    test="name"
    message='''
python CircosConf.py --input circos.config --output pipe.conf

    '''
    print message

#[name, seq, start, strand]
def Supporting_count(event, tsd_start, teSupportingReads):
    total = 0
    right = 0
    left  = 0
    read1 = []
    read2 = []
    if teSupportingReads.has_key(event):
        for read in teSupportingReads[event]:
            name   = read[0]
            seq    = read[1]
            start  = read[2]
            strand = read[3]
            #print name, seq, start, strand
            if int(start) + len(seq) <= int(tsd_start) and strand == '+':
                total += 1
                left  += 1
                read1.append(name)
            elif int(start) >= int(tsd_start) and strand == '-':
                total += 1
                right += 1
                read2.append(name)
        return (total, left, right, ','.join(read1), ','.join(read2))
    else:
        return (0,0,0,'','') 

def read_repeat_name(infiles):
    data = defaultdict(str)
    for infile in infiles:
        with open (infile, 'r') as filehd:
            for line in filehd:
                line = line.rstrip()
                if len(line) > 2: 
                    unit = re.split(r'\t',line)
                    data[unit[0]] = unit[1]
    return data

#s1= re.compile(r'(\S+)\.[rf]')
#s2= re.compile(r'(\S+)\/[12]')
## define insertion as repeat family according reads <-> repeats relation from blat or bam files
def insertion_family_supporting(reads, read_repeat):
    repeat_family = defaultdict(lambda : int())
    for read in re.split(r',', reads):
        read_name  = read
        read_name1 = '%s/1' %(read)
        read_name2 = '%s/2' %(read)
        read_name3 = '%s.f' %(read)
        read_name4 = '%s.r' %(read)
        if read_repeat.has_key(read_name):
            repeat_family[read_repeat[read_name]] += 1
        elif read_repeat.has_key(read_name1):
            repeat_family[read_repeat[read_name1]] += 1
            #print '%s,%s,%s,' %(read_name, read_name1, read_repeat[read_name1])
        elif read_repeat.has_key(read_name2):
            repeat_family[read_repeat[read_name2]] += 1
            #print '%s,%s,%s,' %(read_name, read_name2, read_repeat[read_name2])
        elif read_repeat.has_key(read_name3):
            repeat_family[read_repeat[read_name3]] += 1
        elif read_repeat.has_key(read_name4):
            repeat_family[read_repeat[read_name4]] += 1
    if len(repeat_family.keys()) == 1:
        #return first element if only have one repeat
        return repeat_family.keys()[0]
    elif len(repeat_family.keys()) > 1:
        #return the one have largest value
        sorted_by_value = OrderedDict(sorted(repeat_family.items(), key=lambda x: x[1]))
        return sorted_by_value.keys()[-1]
    else:
        return 'NA'


## define insertion as repeat family according reads <-> repeats relation from blat or bam files
def insertion_family(reads, read_repeat):
    repeat_family = defaultdict(lambda : int())
    r = re.compile(r'(.*):(start|end):(5|3)')
    for read in re.split(r',', reads):
        m = r.search(read)
        read_name = m.groups(0)[0] if m else 'NA'
        if read_name != 'NA' and read_repeat.has_key(read_name):
            repeat_family[read_repeat[read_name]] += 1
    if len(repeat_family.keys()) == 1:
        #return first element if only have one repeat
        return repeat_family.keys()[0]
    elif len(repeat_family.keys()) > 1:
        #return the one have largest value
        sorted_by_value = OrderedDict(sorted(repeat_family.items(), key=lambda x: x[1]))
        return sorted_by_value.keys()[-1]
    else:
        return ''

def write_output(result, read_repeat, usr_target, exper, TE, required_reads, required_left_reads, required_right_reads, teInsertions, teInsertions_reads, teSupportingReads):
    createdir(result)
    NONREF = open ('%s/%s.%s.all_nonref_insert.txt' %(result, usr_target, TE), 'w')
    READS  = open ('%s/%s.%s.reads.list' %(result, usr_target, TE), 'w')
    #teInsertions[event][TSD_seq][TSD_start]['count']   += 1   ## total junction reads
    #teInsertions[event][TSD_seq][TSD_start][pos]       += 1   ## right/left junction reads
    #teInsertions[event][TSD_seq][TSD_start][TE_orient] += 1   ## plus/reverse insertions 
    for event in sorted(teInsertions.keys(), key=int):
        #print 'event: %s' %(event)
        cluster_collection = []
        start_collection = []
        start_both_junction = 0
        for start in sorted(teInsertions[event].keys(), key=int):
            #####supporting reads
            total_supporting, left_supporting, right_supporting = [0, 0, 0]
            left_reads, right_reads = ['', '']
            total_supporting, left_supporting, right_supporting, left_reads, right_reads = Supporting_count(event, start, teSupportingReads)
            repeat_supporting = insertion_family_supporting('%s,%s' %(left_reads, right_reads), read_repeat)
            
            #####tsd within one start position into one start position
            tsd_count         = {}
            TE_orient, repeat_junction = ['','']
            TE_orient_foward, TE_orient_reverse    = [0,0]
            total_count, left_count, right_count   = [0,0,0]
            reads             = []
            #if two or more TSD at this loci, some might be caused by sequencing error, we add the count
            #and use high frequent one to get te orientation and family
            #if only one TSD found, we acturally use only this one
            for foundTSD in sorted(teInsertions[event][start].keys()):
                total_count += teInsertions[event][start][foundTSD]['count']
                left_count  += teInsertions[event][start][foundTSD]['left']
                right_count += teInsertions[event][start][foundTSD]['right']
                TE_orient_foward  += teInsertions[event][start][foundTSD]['+']
                TE_orient_reverse += teInsertions[event][start][foundTSD]['-']
                tsd_count[foundTSD] = teInsertions[event][start][foundTSD]['count']
                reads.extend(teInsertions_reads[event][start][foundTSD]['read'])
            tsd_top = OrderedDict(sorted(tsd_count.items(), key=lambda x: x[1])).keys()[-1] 
            TE_orient         = '+' if int(TE_orient_foward) > int(TE_orient_reverse) else '-'
            repeat_junction   = insertion_family(','.join(reads), read_repeat)

            #####information on each start position
            start_both_junction = start_both_junction + 1 if (int(left_count) > 0 and int(right_count) > 0) else start_both_junction
            start_collection.append([start, foundTSD, total_count, left_count, right_count, repeat_junction, ','.join(reads), total_supporting, left_supporting, right_supporting, left_reads, right_reads, repeat_supporting, TE_orient])

        if len(teInsertions[event].keys()) > 1:
            #two or more start found in one read cluster
            if int(start_both_junction) == len(teInsertions[event].keys()):
                #we keep all if they all supported by junction at both direction
                cluster_collection.extend(start_collection)
            else:
                #not all supported at both end of junction
                if int(start_both_junction) > 0:
                    #we have tsd supported by junction at both direction
                    #we only keep these supported by junction at both direction
                    #even some supported by more than 1 junction reads, we thought due to sequence context
                    #these insertions should have junction at both direction if they are true
                    for start_it in start_collection:
                        if int(start_it[3]) > 0 and int(start_it[4] > 0):
                            cluster_collection.append(start_it)
                else:
                    #do not have tsd supported by junction at both direction
                    #we keep the top one with higher number of junction reads
                    top_start_it = []
                    top_start_junction = 0
                    for start_it in start_collection:
                        if int(top_start_junction) < int(start_it[2]):
                                top_start_junction = int(start_it[2])
                                top_start_it       = start_it
                    cluster_collection.append(top_start_it)
        else:
            #only one tsd found, just add to start_collection
            cluster_collection.extend(start_collection)
         
        #[start, foundTSD, total_count, left_count, right_count, repeat_junction, ','.join(reads),
        #total_supporting, left_supporting, right_supporting, left_reads, right_reads, repeat_supporting, TE_orient]
        for insertion in cluster_collection:
            #####insertion te name
            r_supporting = insertion[12]
            r_junction   = insertion[5]
            repeat_family     = 'NA'
            if r_junction == r_supporting and r_junction != 'NA':
                repeat_family = r_junction
            elif r_junction != r_supporting and r_junction != 'NA' and r_supporting != 'NA':
                repeat_family = '%s/%s' %(r_junction, r_supporting)
            elif r_junction != 'NA':
                repeat_family = r_junction
            elif r_supporting != 'NA':
                repeat_family = r_supporting
            
            #####write to file
            t_count = insertion[2]
            l_count = insertion[3]
            r_count = insertion[4]
            i_start = insertion[0]
            i_tsd   = insertion[1]
            i_reads = insertion[6]
            l_reads = insertion[10]
            r_reads = insertion[11]
            t_supporting = insertion[7]
            l_supporting = insertion[8]
            r_supporting = insertion[9]
            t_orient     = insertion[13]
            
            if int(l_count) >= int(required_left_reads) or int(r_count) >= int(required_right_reads):
                #at lease need one end have junction reads
                coor       = int(i_start) + (len(i_tsd) - 1)
                coor_start = coor - (len(i_tsd) - 1)
                print >> READS, '%s\t%s:%s..%s\tJunction_reads\t%s' %(TE, usr_target, coor_start, coor, i_reads)
                print >> READS, '%s\t%s:%s..%s\tLeft_supporting_reads\t%s' %(TE, usr_target, coor_start, coor, l_reads)
                print >> READS, '%s\t%s:%s..%s\tRight_supporting_reads\t%s' %(TE, usr_target, coor_start, coor, r_reads)
                if int(l_count) > 0 and int(r_count) >0:
                    #both ends have junction reads
                    print >> NONREF, '%s\t%s\t%s\t%s\t%s..%s\t%s\tT:%s\tR:%s\tL:%s\tST:%s\tSR:%s\tSL:%s' %(repeat_family, i_tsd, exper, usr_target, coor_start, coor, t_orient, t_count, r_count, l_count, t_supporting, r_supporting, l_supporting)
                else:
                    #only have end with junction 
                    if int(r_supporting) >= 1 and int(l_supporting) >= 1:
                        #only one end with junction but both end with supporting reads
                        print >> NONREF, '%s\tsupporting_reads\t%s\t%s\t%s..%s\t%s\tT:%s\tR:%s\tL:%s\tST:%s\tSR:%s\tSL:%s' %(repeat_family, exper, usr_target, coor_start, coor, t_orient, t_count, r_count, l_count, t_supporting, r_supporting, l_supporting)
                    else:
                        #only one end with junction and only one end with supporting reads
                        if (int(r_supporting) >= 1 and int(l_count) >= 1) or (int(l_supporting) >= 1 and int(r_count) >= 1):
                            print >> NONREF, '%s\tsupporting_reads\t%s\t%s\t%s..%s\t%s\tT:%s\tR:%s\tL:%s\tST:%s\tSR:%s\tSL:%s' %(repeat_family, exper, usr_target, coor_start, coor, t_orient, t_count, r_count, l_count, t_supporting, r_supporting, l_supporting)
                        elif int(t_supporting) + int(t_count) == 1:
                            #singleton
                            print >> NONREF, '%s\tsingleton\t%s\t%s\t%s..%s\t%s\tT:%s\tR:%s\tL:%s\tST:%s\tSR:%s\tSL:%s' %(repeat_family, exper, usr_target, coor_start, coor, t_orient, t_count, r_count, l_count, t_supporting, r_supporting, l_supporting)
                        else:
                            #insufficient
                            print >> NONREF, '%s\tinsufficient_data\t%s\t%s\t%s..%s\t%s\tT:%s\tR:%s\tL:%s\tST:%s\tSR:%s\tSL:%s' %(repeat_family, exper, usr_target, coor_start, coor, t_orient, t_count, r_count, l_count, t_supporting, r_supporting, l_supporting)
            else:
                #no junction reads
                if int(r_supporting) >= 1 and int(l_supporting) >= 1:
                    #no junction reads, but both end with supporting reads
                    print >> NONREF, '%s\tsupporting_reads\t%s\t%s\t%s..%s\t%s\tT:%s\tR:%s\tL:%s\tST:%s\tSR:%s\tSL:%s' %(repeat_family, exper, usr_target, coor_start, coor, t_orient, t_count, r_count, l_count, t_supporting, r_supporting, l_supporting)
                elif int(t_supporting) == 1:
                    print >> NONREF, '%s\tsingleton\t%s\t%s\t%s..%s\t%s\tT:%s\tR:%s\tL:%s\tST:%s\tSR:%s\tSL:%s' %(repeat_family, exper, usr_target, coor_start, coor, t_orient, t_count, r_count, l_count, t_supporting, r_supporting, l_supporting)
                else:
                    print >> NONREF, '%s\tinsufficient_data\t%s\t%s\t%s..%s\t%s\tT:%s\tR:%s\tL:%s\tST:%s\tSR:%s\tSL:%s' %(repeat_family, exper, usr_target, coor_start, coor, t_orient, t_count, r_count, l_count, t_supporting, r_supporting, l_supporting)
                #pass


def TSD_from_read_depth(r, read_repeat, teReadClusters, teReadClusters_count, teReadClusters_depth, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found):
    #determine TSD from read depth at insertions site
    #count depth to find TSD in
    #if there are 5 reads (2 right, 3 left) they
    #should only be a depth of 5 at the TSD
    #teReadCluster = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : str()))))
    #teReadClusters[event]['read_inf'][name]['strand']= strand
    #teReadClusters_count[event]['read_count'] += 1
    #teReadClusters_depth[event]['read_inf']['depth'][i] += 1 
    for cluster in sorted(teReadClusters.keys(), key=int):
        read_total = teReadClusters_count[cluster]['read_count']
        TSD_len = 0 
        for chrs_pos in sorted(teReadClusters_depth[cluster]['read_inf']['depth'].keys(), key=int):
            depth = teReadClusters_depth[cluster]['read_inf']['depth'][chrs_pos]
            if float(depth) >= 0.8*float(read_total):
                TSD_len += 1
        if TSD_len > 0:
            TSD = '.'*TSD_len
            for name in teReadClusters[cluster]['read_inf'].keys():
                real_name = r.search(name).groups(0)[0] if r.search(name) else ''
                seq    = teReadClusters[cluster]['read_inf'][name]['seq']
                start  = teReadClusters[cluster]['read_inf'][name]['start']
                strand = teReadClusters[cluster]['read_inf'][name]['strand']
                #print '%s\t%s\t%s\t%s\t%s' %(cluster, start, name, TSD, strand)
                TSD_check(cluster, seq, start, real_name, read_repeat, name, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found)
            del teReadClusters[cluster]
            del teReadClusters_count[cluster]
            del teReadClusters_depth[cluster]
        else:
            pass
            #what if we can not find TSD? still could be insertions


def align_process(bin_ins, read_repeat, record, r, r_tsd, count, seq, start, end, name, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found, teReadClusters, teReadClusters_count, teReadClusters_depth, teSupportingReads):
    range_allowance = 500
    padded_start    = bin_ins[0] - range_allowance
    padded_end      = bin_ins[-1] + range_allowance 
    #insertions
    #print 'insertions: %s' %(name)
    if (int(start) >= padded_start and int(start) <= padded_end) or (int(end) >= padded_start and int(end) <= padded_end):
        bin_ins.extend([int(start), int(end)])
        bin_ins = sorted(bin_ins, key=int)
        if r.search(name):
            real_name = r.search(name).groups(0)[0]
            if not r_tsd.search(TSD):
                TSD_check(count, seq, start, real_name, read_repeat, name, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found)
            else:
                calculate_cluster_depth(count, seq, start, name, strand, teReadClusters, teReadClusters_count, teReadClusters_depth)
        elif not r.search(name) and not record.is_paired:
            #reads not matched to repeat and not mates of junctions
            #reads are mates of reads matched to middle of repeat
            #supporting reads
            teSupportingReads[count].append([name, seq, start, strand])
    else:
        #if start and end do not fall within last start and end
        #we now have a different insertion event
        count += 1
        if r.search(name):
            real_name = r.search(name).groups(0)[0]
            if not r_tsd.search(TSD):
                TSD_check(count, seq, start, real_name, read_repeat, name, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found)
            else:
                calculate_cluster_depth(count, seq, start, name, strand, teReadClusters, teReadClusters_count, teReadClusters_depth)
        elif not r.search(name) and not record.is_paired:
            #reads not matched to repeat and not mates of junctions
            #reads are mates of reads matched to middle of repeat
            #supporting reads
            teSupportingReads[count].append([name, seq, start, strand])
        #initial insertion site boundary
        bin_ins = [int(start), int(end)]
        #print '%s\t%s' %(count, bin_ins)
    return (bin_ins, count)

def existingTE(infile, existingTE_inf, existingTE_found):
    r = re.compile(r'(\S+)\t(\S+):(\d+)\.\.(\d+)')
    if infile != 'NONE':
        blat = 0
        with open(infile, 'r') as filehd:
            for line in filehd:
                line = line.rstrip() 
                if not blat:
                    if r.search(line):
                        te, chrs, start, end = r.search(line).groups(0)
                        start, end = sorted([start, end], key=int)
                        existingTE_inf[te]['start'][int(start)] = line
                        existingTE_inf[te]['end'][int(end)]     = line
                        existingTE_found[line]['start']= 0
                        existingTE_found[line]['end']  = 0

def existingTE_RM(infile, existingTE_inf):
    r_end = re.compile(r'\((\d+)\)')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\s+',line)
                #print line
                #print unit[5], unit[9], unit[12], unit[13], unit[14]
                if unit[9] == '+':
                    if int(unit[12]) == 1:
                        existingTE_inf[unit[10]]['start'][int(unit[6])] = 1
                        #print unit[10], 'start', unit[6]
                    if len(unit[14]) == 3:
                        unit[14] =re.sub(r'\(|\)', '', unit[14])
                        if int(unit[14]) == 0:
                            existingTE_inf[unit[10]]['end'][int(unit[7])] = 1
                            #print unit[10], 'end', unit[7]
                elif unit[9] == 'C':
                    if len(unit[12]) == 3:
                        unit[12] =re.sub(r'\(|\)', '', unit[12])
                        if int(unit[12]) == 0:
                            existingTE_inf[unit[10]]['start'][int(unit[6])] = 1
                            #print unit[10], 'start', unit[6]
                    if int(unit[14]) == 1:
                        existingTE_inf[unit[10]]['end'][int(unit[7])] = 1
                        #print unit[10], 'end', unit[7]
 

def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(seq)
    for i in range(len(bases)):
        bases[i] = complement[bases[i]] if complement.has_key(bases[i]) else bases[i]
    return ''.join(bases)

def reverse_complement(seq):
    return complement(seq[::-1])

def calculate_cluster_depth(event, seq, start, name, strand, teReadClusters, teReadClusters_count, teReadClusters_depth):
    teReadClusters_count[event]['read_count'] += 1
    teReadClusters[event]['read_inf'][name]['seq']   = seq
    teReadClusters[event]['read_inf'][name]['start'] = start
    teReadClusters[event]['read_inf'][name]['strand']= strand
    for i in range(int(start), int(start)+len(seq)):
        teReadClusters_depth[event]['read_inf']['depth'][i] += 1

def TSD_check(event, seq, start, real_name, read_repeat, name, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found):
    ##TSD already specified by usr, not unknown
    ##seq is entire trimmd read, not just the TSD portion of the read
    ##start is the first postition of the entire read match to ref
    repeat = read_repeat[real_name] # need to deal with any te, get infor from name of reads
    rev_com = reverse_complement(seq)
    result    = 0
    pos       = ''
    TE_orient = 0
    TSD_start = 0
    TSD_seq   = ''
    r5 = re.compile(r'start:[53]$')
    r3 = re.compile(r'end:[53]$')
    r5_tsd = re.compile(r'^(%s)' %(TSD))
    r3_tsd = re.compile(r'(%s)$' %(TSD))
    ##start means that the TE was removed from the start of the read
    ##5 means the trimmed end mapps to the 5prime end of the TE
    ##3 means the trimmed end mapps to the 3prime end of the TE
    if strand == '+':
        if r5.search(name) and (r5_tsd.search(seq) or r3_tsd.search(rev_com)):
            result    = 1
            TSD_seq   = r5_tsd.search(seq).groups(0)[0] if r5_tsd.search(seq) else 'UNK'
            pos       = 'right'
            TE_orient = '-' if name[-1] == '5' else '+'
            TSD_start = int(start)
        elif r3.search(name) and (r5_tsd.search(rev_com) or r3_tsd.search(seq)):
            result    = 1
            TSD_seq   = r3_tsd.search(seq).groups(0)[0] if r3_tsd.search(seq) else 'UNK'
            pos       = 'left'
            TE_orient = '+' if name[-1] == '5' else '-'
            TSD_start = int(start) + (len(seq)-len(TSD))
    elif strand == '-':
        if r5.search(name) and (r5_tsd.search(rev_com) or r3_tsd.search(seq)):
            result    = 1
            TSD_seq   = r3_tsd.search(seq).groups(0)[0] if r3_tsd.search(seq) else 'UNK'
            pos       = 'left'
            TE_orient = '+' if name[-1] == '5' else '-'
            TSD_start = int(start) + (len(seq)-len(TSD))
        elif r3.search(name) and (r5_tsd.search(seq) or r3_tsd.search(rev_com)):
            result    = 1
            TSD_seq   = r5_tsd.search(seq).groups(0)[0] if r5_tsd.search(seq) else 'UNK'
            pos       = 'right'
            TE_orient = '-' if name[-1] == '5' else '+'
            TSD_start = int(start)
    print '%s\t%s\t%s\t%s\t%s\t%s\t%s' %(event, name, TSD_seq, TSD_start, TE_orient, pos, repeat)
    if result and TE_orient:
        tir1_end, tir2_end = [0, 0]
        if pos == 'left':
            tir1_end = int(start) + len(seq)
            print 'tir1: %s' %(tir1_end)
        elif pos == 'right':
            tir2_end = int(start) - 1
            print 'tir2: %s' %(tir2_end)
        if tir1_end > 0 and existingTE_inf[repeat]['start'].has_key(tir1_end):
            te_id = existingTE_inf[repeat]['start'][tir1_end]
            #existingTE_found[te_id]['start'] += 1
            print 'tir1'
        elif tir2_end > 0 and existingTE_inf[repeat]['end'].has_key(tir2_end):
            te_id = existingTE_inf[repeat]['end'][tir2_end]
            #existingTE_found[te_id]['end'] += 1
            print 'tir2'
        else:
            print 'not match'
            ##non reference insertions
            teInsertions[event][TSD_start][TSD_seq]['count']   += 1   ## total junction reads
            teInsertions[event][TSD_start][TSD_seq][pos]       += 1   ## right/left junction reads
            teInsertions[event][TSD_start][TSD_seq][TE_orient] += 1   ## plus/reverse insertions
            #read_name = re.sub(r':start|:end', '', name)
            teInsertions_reads[event][TSD_start][TSD_seq]['read'].append(name)
            #print '1: %s\t 2: %s' %(read_name, teInsertions_reads[event][TSD_seq][TSD_start]['read'])
            #print 'C: %s\t%s\t%s\t%s\t%s' %(event, name, TSD_seq, TSD_start, TE_orient)

def convert_tag(tag):
    tags = {}
    for t in tag:
        tags[t[0]] = t[1]
    return tags

def find_insertion_cluster_bam(align_file, read_repeat, target, TSD, teInsertions, teInsertions_reads, teReadClusters, teReadClusters_count, teReadClusters_depth, existingTE_inf, existingTE_found, teSupportingReads):
    r = re.compile(r'(.*):(start|end):(5|3)')
    r_tsd = re.compile(r'UNK|UKN|unknown', re.IGNORECASE)
    r_cg  = re.compile(r'[SID]')
    bin_ins        = [0]
    count          = 0
    TSD_len        = len(TSD)

    ref  = 'None' if target == 'ALL' else target
    fsam = pysam.AlignmentFile(align_file, 'rb')
    for record in fsam.fetch(reference=ref, until_eof = True):
        if not record.is_unmapped:
            name   = record.query_name
            flag   = record.flag
            start  = int(record.reference_start) + 1
            MAPQ   = record.mapping_quality
            cigar  = record.cigarstring
            seq    = record.query_sequence
            tag    = record.tags if record.tags else []
            length = len(seq)
            end    = int(start) + int(length) - 1 #should not allowed for indel or softclip
            strand = ''
            # flag is 0 is read if read is unpaired and mapped to plus strand
            if int(flag) == 0:
                strand = '+'
            else:
                strand = '-' if record.is_reverse else '+'
            if r_cg.search(cigar):
                continue
            tags = convert_tag(tag)
            #print '%s\t%s\t%s' %(name, start, length)
            # filter low quality mapping reads: 
            # 1. paired-end reads at least have one reads unique mapped (MAPQ set to 0 for both reads if both are repeat, else should be > 0 at least one unique mapped)
            # 2. unpaired reads should unique mapped, no gap, mismatch <= 3 and no suboptimal alignment
            #print 'before: %s\t%s\t%s' %(name, count, bin_ins)
            #if record.is_proper_pair and (int(MAPQ) >= 29 or tags['XT'] == 'U'):
            if record.is_proper_pair and int(MAPQ) > 0:
                bin_ins, count = align_process(bin_ins, read_repeat, record, r, r_tsd, count, seq, start, end, name, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found, teReadClusters, teReadClusters_count, teReadClusters_depth, teSupportingReads)
            elif not record.is_paired:
                #if tags['XT'] == 'U' and int(tags['XO']) == 0 and int(tags['XM']) <= 3 and int(tags['X1']) == 0:
                if tags['XT'] == 'U' and int(tags['XO']) == 0 and (int(tags['XM']) <= 3 or int(tags['X1']) == 0):
                    bin_ins, count = align_process(bin_ins, read_repeat, record, r, r_tsd, count, seq, start, end, name, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found, teReadClusters, teReadClusters_count, teReadClusters_depth, teSupportingReads)
            #print 'after: %s\t%s\t%s' %(name, count, bin_ins)

    ###TSD not given we infer from read depth
    if r_tsd.search(TSD):
        TSD_from_read_depth(r, read_repeat, teReadClusters, teReadClusters_count, teReadClusters_depth, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found)        
        
def find_insertion_cluster_sam(sorted_align, TSD, teInsertions, teInsertions_reads, teReadClusters, teReadClusters_count, teReadClusters_depth, existingTE_inf, existingTE_found):
    align_seg = pysam.AlignedSegment()
    r = re.compile(r':(start|end):(5|3)')
    r_tsd = re.compile(r'UNK|UKN|unknown', re.IGNORECASE)
    bin_ins        = [0]
    count          = 0
    TSD_len        = len(TSD)
    ##go through the reads, cluster reads and determine junction reads and supporting reads
    for record in sorted_align.keys():
        name, flag, ref, start, MAPQ, cigar, MRNM, MPOS, TLEN, seq, qual, tag = ['']*12
        try:
            name, flag, ref, start, MAPQ, cigar, MRNM, MPOS, TLEN, seq, qual, tag = re.split(r'\t', record, 11)
        except:
            try:
                name, flag, ref, start, MAPQ, cigar, MRNM, MPOS, TLEN, seq, qual = re.split(r'\t', record, 10)
            except:
                print 'sam format error, can not split into 11 or 12 field'
                exit()
        #tags = re.split(r'\t', tags)
        if int(MAPQ) < 40:
            continue
        end = len(seq) + int(start) - 1 # should not allowed for indel or softclip
        align_seg.flag = int(flag)
        strand = ''
        if int(flag) == 0:
            strand = '+'
        else:
            strand = '-' if align_seg.is_reverse else '+'
        #print name, TLEN, tags[0]
        range_allowance = 0
        padded_start    = bin_ins[0] - range_allowance
        padded_end      = bin_ins[-1] - range_allowance 
        m = r.search(name)
        if m:
            #insertions
            #print 'insertions: %s' %(name)
            if (int(start) >= padded_start and int(start) <= padded_end) or (int(end) >= padded_start and int(end) <= padded_end):
                bin_ins.extend([int(start), int(end)])
                bin_ins = sorted(bin_ins, key=int)
                if not r_tsd.search(TSD):
                    TSD_check(count, seq, start, name, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found)
                else:
                    calculate_cluster_depth(count, seq, start, name, strand, teReadClusters, teReadClusters_count, teReadClusters_depth)
            else:
                #if start and end do not fall within last start and end
                #we now have a different insertion event
                count += 1
                if not r_tsd.search(TSD):
                    TSD_check(count, seq, start, name, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found)
                else:
                    calculate_cluster_depth(count, seq, start, name, strand, teReadClusters, teReadClusters_count, teReadClusters_depth)
                #initial insertion site boundary
                bin_ins = [int(start), int(end)]
        else:
            #supporting reads
            #print 'supportings: %s' %(name)
            pass
        #print name, flag, strand
    if r_tsd.search(TSD):
        TSD_from_read_depth(teReadClusters, teReadClusters_count, teReadClusters_depth, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found)

##default sam, not deal with paired end mapping and supporting reads 
def remove_redundant_sam(infile, target):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'@'):
                unit = re.split(r'\t',line)
                if target != 'ALL' and unit[2] != target:
                    continue
                else:
                    start = unit[3]
                    if not data.has_key(line):
                        data[line] = start
    data = OrderedDict(sorted(data.items(), key=lambda x: x[1]))
    return data

def createdir(dirname):
    if not os.path.exists(dirname):
        os.mkdir(dirname)

def parse_regex(infile):
    data = []
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\s+',line)
                unit[0] = re.sub(r'.(fq|fastq)','',unit[0])
                unit[1] = re.sub(r'.(fq|fastq)','',unit[1])
                unit[2] = re.sub(r'.(fq|fastq)','',unit[2])
                data = unit
    return data


def main():

    required_reads       = 1           ## rightReads + leftReads needs to be > to this value
    required_left_reads  = 1           ## needs to be >= to this value
    required_right_reads = 1           ## needs to be >= to this value
    align_file           = sys.argv[1] ## combined bowtie or bwa results, sam format only 
    usr_target           = sys.argv[2] ## chromosome to analyze: ALL or Chr1..N
    genome_path          = sys.argv[3] ## genome sequence
    TE                   = sys.argv[4] ## repeat to analyze: ALL or mPing/other te name 
    regex_file           = sys.argv[5] ## regex.txt
    exper                = sys.argv[6] ## prefix for output, title: HEG4
    flank_len            = sys.argv[7] ## length of seq flanking insertions to be returned: 100
    existing_TE          = sys.argv[8] ## existingTE.blatout
    mm_allow             = sys.argv[9] ## mismatches allowed: 0, 1, 2, 3
    bowtie2              = sys.argv[10] ## use bowtie2 or not: 1 or 0
    #relax_reference      = sys.argv[11]## relax mode for existing TE: 1 or 0
    #relax_align          = sys.argv[12]## relax mode for insertion: 1 or 0
    bowtie_sam           = 1           ## change to shift or remove in V2
    existingTE_inf       = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : str)))
    existingTE_found     = defaultdict(lambda : defaultdict(lambda : int))
    bwa                  = 0
    teInsertions         = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : int()))))
    teInsertions_reads   = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : list()))))
    teReadClusters       = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : str()))))
    teReadClusters_count = defaultdict(lambda : defaultdict(lambda : int()))
    teReadClusters_depth = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : int()))))
    teSupportingReads    = defaultdict(lambda : list())

    #read existing TE from file
    r_te = re.compile(r'repeatmasker|rm|\.out', re.IGNORECASE)
    if os.path.isfile(existing_TE) and os.path.getsize(existing_TE) > 0:
        if r_te.search(existing_TE):
            existingTE_RM(existing_TE, existingTE_inf)
        else:
            existingTE(existing_TE, existingTE_inf, existingTE_found)
    else:
        print 'Existing TE file does not exists or zero size'

    ##get the regelar expression patterns for mates and for the TE
    ##when passed on the command line as an argument, even in single
    ##quotes I lose special regex characters
    s         = re.compile(r'[\[.*+?]')
    mate_file = parse_regex(regex_file)
    TSD       = mate_file[3]
    TSDpattern= 1 if s.search(TSD) else 0

    
    ##remove redundant alignment, no need to if we used combined bam file from previous step
    #sorted_align = defaultdict(lambda : str)
    #if bwa == 1:
    #    sorted_align = remove_redundant_sam(align_file, usr_target)
    #elif bowtie_sam == 1:
    #    sorted_align = remove_redundant_sam(align_file, usr_target)
    #else:
    #    print 'We accept only sam format right now'
    #    usage() 
    #    exit(2)

    ##read -> repeat relation
    top_dir = re.split(r'/', os.path.dirname(os.path.abspath(align_file)))[:-1]
    result  = '%s/results' %('/'.join(top_dir))
    read_repeat_files = glob.glob('%s/te_containing_fq/*.read_repeat_name.txt' %('/'.join(top_dir)))
    read_repeat = read_repeat_name(read_repeat_files)

    ##cluster reads around insertions
    #find_insertion_cluster_sam(sorted_align, TSD, teInsertions, teInsertions_reads, teReadClusters, teReadClusters_count, teReadClusters_depth, existingTE_inf, existingTE_found)
    find_insertion_cluster_bam(align_file, read_repeat, usr_target, TSD, teInsertions, teInsertions_reads, teReadClusters, teReadClusters_count, teReadClusters_depth, existingTE_inf, existingTE_found, teSupportingReads)


    ##output insertions
    #top_dir = re.split(r'/', os.path.dirname(os.path.abspath(align_file)))[:-1]
    #result  = '%s/results' %('/'.join(top_dir))
    #read_repeat_files = glob.glob('%s/te_containing_fq/*.read_repeat_name.txt' %('/'.join(top_dir)))
    #read_repeat = read_repeat_name(read_repeat_files)
    write_output(result, read_repeat, usr_target, exper, TE, required_reads, required_left_reads, required_right_reads, teInsertions, teInsertions_reads, teSupportingReads)

 
if __name__ == '__main__':
    main()

