from src.classes_another import Bkp, Single_same_Multi_junc
from src.cluster_op import compare_two_bkp
from src.other_op import get_min_junc_support, get_min_multijunc_support

from Bio import pairwise2
from Bio.Seq import Seq

def get_seq(junc, bkp_index):
    if junc.junc_sub_type == 'INS':
        for aligns in junc.support_reads:
            for align in aligns:
                if align.chr == None:
                    return align.seq
    else:
        for aligns in junc.support_reads:
            for align in aligns:
                if align.chr == junc.include_bkps[bkp_index].chrom:
                    if junc.include_bkps[bkp_index].direct == 'L' and abs(junc.include_bkps[bkp_index].pos-align.ref_start)<1000:
                        return align.seq[:5000]
                    if junc.include_bkps[bkp_index].direct == 'S' and abs(junc.include_bkps[bkp_index].pos-align.ref_end)<1000:
                        return align.seq[len(align.seq)-5000:]
    return None
    

def compare_juncSe_for_multi_junc(juncSe1, juncSe2, seq_type):
    for i in range(len(juncSe1.include_juncs)):
        for l in range(2):
            for j in range(len(juncSe2.include_juncs)):
                for m in range(2):
                    if compare_two_bkp(juncSe1.include_juncs[i].include_bkps[l], juncSe2.include_juncs[j].include_bkps[m]):
                        if seq_type == 'ngs':
                            return 1, Bkp(juncSe1.include_juncs[i].include_bkps[l].chrom, int((juncSe1.include_juncs[i].include_bkps[l].pos+juncSe2.include_juncs[j].include_bkps[m].pos)/2), juncSe1.include_juncs[i].include_bkps[l].direct, juncSe1.include_juncs[i].include_bkps[l].type), juncSe1.include_juncs[i].include_bkps[1-l], juncSe2.include_juncs[j].include_bkps[1-m], []
                        seq1 = get_seq(juncSe1.include_juncs[i], 1-l)
                        seq2 = get_seq(juncSe2.include_juncs[j], 1-m)
                        if seq1 == None or seq2 == None:
                            return 0, None, None, None, []
                        seq1 = Seq(seq1)
                        seq2 = Seq(seq2)
                        alignments = pairwise2.align.globalxx(seq1, seq2)
                        srt_alignments = sorted(alignments, key=lambda x:x.score, reverse=True)
                        ratio = srt_alignments[0].score/min(len(seq1), len(seq2))
                        seqs = [srt_alignments[0].seqA, srt_alignments[0].seqB]
                        if ratio > 0.8:
                            return 1, Bkp(juncSe1.include_juncs[i].include_bkps[l].chrom, int((juncSe1.include_juncs[i].include_bkps[l].pos+juncSe2.include_juncs[j].include_bkps[m].pos)/2), juncSe1.include_juncs[i].include_bkps[l].direct, juncSe1.include_juncs[i].include_bkps[l].type), \
                                juncSe1.include_juncs[i].include_bkps[1-l], juncSe2.include_juncs[j].include_bkps[1-m], seqs
                        else:
                            return 0, Bkp(juncSe1.include_juncs[i].include_bkps[l].chrom, int((juncSe1.include_juncs[i].include_bkps[l].pos+juncSe2.include_juncs[j].include_bkps[m].pos)/2), juncSe1.include_juncs[i].include_bkps[l].direct, juncSe1.include_juncs[i].include_bkps[l].type), \
                                juncSe1.include_juncs[i].include_bkps[1-l], juncSe2.include_juncs[j].include_bkps[1-m], seqs

    return 0, None, None, None, []

def get_no_and_multi_juncSe(juncSe_list, seq_type):
    no_multi_juncSe = [] # for bridge search
    multi_juncSe = []
    multi_SV_juncSe = []

    for i in range(len(juncSe_list)):
        juncSe = juncSe_list[i]
        get_multi_flag = 0
        for target_juncSe in no_multi_juncSe:
            same_flag, same_bkp_info, juncSe1_info, juncSe2_info, seqs = compare_juncSe_for_multi_junc(target_juncSe, juncSe, seq_type)
            if same_flag == 1:
                multi_juncSe.append(Single_same_Multi_junc(same_bkp_info, juncSe1_info, juncSe2_info, [target_juncSe, juncSe], [target_juncSe.support_reads, juncSe.support_reads], seqs))
                get_multi_flag = 1
                break
            elif same_flag == 0 and same_bkp_info != None:
                multi_SV_juncSe.append(Single_same_Multi_junc(same_bkp_info, juncSe1_info, juncSe2_info, [target_juncSe, juncSe], [target_juncSe.support_reads, juncSe.support_reads], seqs))
                get_multi_flag = 1
                break
        if get_multi_flag == 0:
            no_multi_juncSe.append(juncSe)
        else:
            no_multi_juncSe.remove(target_juncSe)
    
    return no_multi_juncSe, multi_juncSe, multi_SV_juncSe

def get_no_and_multi_junc(juncSe_list, seq_type):
    no_multi_junc = []
    multi_junc = []

    all_juncs = []
    for juncSe in juncSe_list:
        for junc in juncSe.include_juncs:
            if len(junc.support_reads) < 2:
                continue
            junc.updated_support_reads = juncSe.get_support_num()
            all_juncs.append(junc)
    
    merge_juncs = []
    for junc in all_juncs:
        flag = 0
        for merge_junc in merge_juncs:
            if (compare_two_bkp(junc.include_bkps[0], merge_junc.include_bkps[0]) and compare_two_bkp(junc.include_bkps[1], merge_junc.include_bkps[1])) or \
                (compare_two_bkp(junc.include_bkps[0], merge_junc.include_bkps[1]) and compare_two_bkp(junc.include_bkps[1], merge_junc.include_bkps[0])):
                merge_junc.support_reads.extend(junc.support_reads)
                merge_junc.updated_support_reads = merge_junc.updated_support_reads + junc.updated_support_reads
        if flag == 0:
            merge_juncs.append(junc)
            
    del all_juncs
    print(len(merge_juncs))

    for junc in merge_juncs:
        if junc.include_bkps[0].chrom == 'chr1' and 152580000<junc.include_bkps[0].pos<152620000:
            print("nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn")

    for l in range(len(merge_juncs)):
        print(l)
        junc = merge_juncs[l]
        flag = 0
        for i in range(2):
            for k in no_multi_junc:
                target_junc = merge_juncs[k]
                for j in range(2):
                    if compare_two_bkp(junc.include_bkps[i], target_junc.include_bkps[j]):
                        # if seq_type == 'ngs':
                        #     return 1, Bkp(junc.include_bkps[0].chrom, int((junc.include_bkps[0].pos+target_junc.include_bkps[j].pos)/2), junc.include_bkps[i].direct, junc.include_bkps[i].type), junc.include_bkps[1-i], target_junc.include_bkps[1-j], []
                        seq1 = get_seq(junc, 1-i)
                        seq2 = get_seq(target_junc, 1-j)
                        # if seq1 == None or seq2 == None:
                        #     continue
                        seq1 = Seq(seq1)
                        seq2 = Seq(seq2)
                        try:
                            alignments = pairwise2.align.globalxx(seq1, seq2)
                        except Exception as e:
                            print(f"Error during alignment: {e}")
                            continue
                        # alignments = pairwise2.align.globalxx(seq1, seq2)
                        srt_alignments = sorted(alignments, key=lambda x:x.score, reverse=True)
                        ratio = srt_alignments[0].score/min(len(seq1), len(seq2))
                        seqs = [srt_alignments[0].seqA, srt_alignments[0].seqB]
                        flag = 1
                        if ratio > 0.8:
                            multi_junc.append(Single_same_Multi_junc(Bkp(junc.include_bkps[i].chrom, int((junc.include_bkps[i].pos+target_junc.include_bkps[j].pos)/2), junc.include_bkps[i].direct, junc.include_bkps[i].type), \
                                                                        [junc.include_bkps[1-i], target_junc.include_bkps[1-j]], junc.updated_support_reads+target_junc.updated_support_reads, seqs))
                            no_multi_junc.remove(k)
                            break
                if flag == 1:
                    break
            if flag == 1:
                break
            for multijunc in multi_junc:
                if compare_two_bkp(junc.include_bkps[i],multijunc.same_bkp_info):
                    new_seq = get_seq(junc, 1-i)
                    multijunc.add(junc.include_bkps[i], new_seq, junc.updated_support_reads) 
                    flag = 1
                    break
            if flag == 1:
                break
        if flag == 0:
            print(l)
            no_multi_junc.append(l)
        print("got")
    
    print("finished")
    print(len(no_multi_junc))
    
    filter_no_multi_junc = []
    filter_multi_junc = []
    for index in no_multi_junc:
        junc = merge_juncs[index]
        # min_support = get_min_junc_support(junc, seq_type)
        if junc.updated_support_reads >= 3:
            filter_no_multi_junc.append(junc)
    if seq_type == 'hifi':
        min_support = 5
    elif seq_type == 'ont':
        min_support = 10
    for juncs in multi_junc:
        # min_support = get_min_multijunc_support(juncs, seq_type)
        if juncs.support_reads_num >= min_support:
            filter_multi_junc.append(juncs)
    
    return filter_no_multi_junc, filter_multi_junc
            
