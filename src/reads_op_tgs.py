import pysam
import re
import os
import pickle

from src.classes_another import Align, Junction_Series
from src.detect_op_another import detect_intra_inter_trans, judge_noisy_aligns
from src.cluster_op import cluster_juncSe_list
from src.other_op import get_mismatch_ratio

def get_supplementary_aligns(align, bam, all_possible_chrs, min_mapq):

    try:
        sa_tag = align.get_tag("SA").split(";")
    except KeyError:
        return []

    supplementary_aligns = []
    # un_flag = 0
    # for element in sa_tag:
    #     fields = element.split(",")
        # if len(fields) != 6:
        #     continue
        # if fields[0] not in all_possible_chrs:
        #     un_flag = 1

    for element in sa_tag:

        fields = element.split(",")
        if len(fields) != 6:
            continue

        rname = fields[0]
        pos = int(fields[1])
        strand = fields[2]
        cigar = fields[3]
        mapq = int(fields[4])

        nm = int(fields[5])

        if mapq < min_mapq:
            continue

        if rname not in all_possible_chrs:
            continue

        # Generate an aligned segment from the information
        a = pysam.AlignedSegment()
        a.query_name = align.query_name
        a.query_sequence = align.query_sequence
        if strand == "+":
            a.flag = 2048
        else:
            a.flag = 2064
        a.reference_id = bam.get_tid(rname)

        a.reference_start = pos - 1

        # num = 0
        # opVals = re.findall(r'(\d+)([\w=])', cigar)
        # for i in range(len(opVals)):
        #     if opVals[i][1] != 'S' and opVals[i][1] != 'H':
        #         num += int(opVals[i][0])

        # if un_flag == 1 and num < 500:
        #     continue

        try:
            a.mapping_quality = mapq
        except OverflowError:
            a.mapping_quality = 0
        a.cigarstring = cigar
        a.next_reference_id = -1
        a.next_reference_start = -1
        a.template_length = 0
        a.query_qualities = align.query_qualities
        a.set_tags([("NM", nm, "i")])

        supplementary_aligns.append(a)

    return supplementary_aligns

def traverse_reads_tgs(bam_path, chrom, region_start, region_end, repeat_regions, options):
    """
    tranverse tgs reads to detect juncSe
    :return:
    """

    print("[Processing {}:{}-{}]".format(chrom, region_start, region_end))
    if not os.path.exists(options.out_path+'/reads_mismatch_ratio_data'):
        os.mkdir(options.out_path+'/reads_mismatch_ratio_data')
    reads_mismatch_ratio_file = open(options.out_path+'/reads_mismatch_ratio_data/'+'{}_{}_{}.txt'.format(chrom, region_start, region_end),'w')
    origin_juncSe_list = []

    bam_file = pysam.AlignmentFile(bam_path)
    aligns = bam_file.fetch(chrom, region_start, region_end)

    for primary in aligns:
        mapq = 0
        filter_flag = 'PASS'
        read_aligns = []

        # # STEP: filter, only proc primary aligns
        if primary.is_supplementary or primary.cigarstring is None or primary.is_unmapped or primary.is_secondary:
            continue
        # unmapped or secondary or low mapping quality, then pass this align
        mapq = primary.mapping_quality
        if primary.mapping_quality < options.min_mapq:
            filter_flag = 'FILTER'

        if options.mismatch_filter:
            align_mismatch_ration = get_mismatch_ratio(primary.cigarstring, primary.reference_start, primary.reference_end)
            # if align_mismatch_ration >= options.max_mismatch_ratio:
            #     filter_flag = 'FILTER'
        # align to a ref that not in genome reference
        ref_name = bam_file.getrname(primary.reference_id)
        if ref_name not in options.all_possible_chrs:
            print("[Warning]: '{0}' not in reference's .fa file, skip this read".format(ref_name))
            continue

        # # STEP: for primary align, get SA and adjust cords
        supplementaries = get_supplementary_aligns(primary, bam_file, options.all_possible_chrs, options.min_mapq)
        if len(supplementaries) > 3:
            continue
        if len(supplementaries) == 0:
            # continue
            if options.mismatch_filter:
                align_mismatch_ration_pr = get_mismatch_ratio(primary.cigarstring, primary.reference_start, primary.reference_end)
                # if align_mismatch_ration_pr >= options.max_mismatch_ratio:
                #     continue
            read_name = primary.query_name
            pm_seq = primary.query_sequence
            align_chr = bam_file.getrname(primary.reference_id)
            if align_chr not in options.all_possible_chrs:
                continue
            open_flag = 0
            opVals = re.findall(r'(\d+)([\w=])', primary.cigarstring)
            pre_sum = 0
            for op in opVals:
                if op[1] == 'I' and int(op[0])>10000:
                    align_obj1 = Align(read_name, align_chr, primary.reference_start, primary.reference_start+pre_sum, primary.query_alignment_start, primary.query_alignment_start+pre_sum, len(primary.query_sequence), 'T', pm_seq[primary.query_alignment_start: primary.query_alignment_start+pre_sum], "AP", filter_flag, mapq, align_mismatch_ration_pr)
                    align_obj2 = Align(read_name, align_chr, primary.reference_start+pre_sum, primary.reference_end, primary.query_alignment_start+pre_sum+int(op[0]), primary.query_alignment_end, len(primary.query_sequence), 'T', pm_seq[primary.query_alignment_start+pre_sum+int(op[0]): primary.query_alignment_end], "AP", filter_flag, mapq, align_mismatch_ration_pr)
                    read_aligns.append(align_obj1)
                    read_aligns.append(align_obj2)
                    open_flag = 1
                    break
                elif op[1] == 'D' and int(op[0])>10000:
                    align_obj1 = Align(read_name, align_chr, primary.reference_start, primary.reference_start+pre_sum, primary.query_alignment_start, primary.query_alignment_start+pre_sum, len(primary.query_sequence), 'T', pm_seq[primary.query_alignment_start: primary.query_alignment_start+pre_sum], "AP", filter_flag, mapq, align_mismatch_ration_pr)
                    align_obj2 = Align(read_name, align_chr, primary.reference_start+pre_sum+int(op[0]), primary.reference_end, primary.query_alignment_start+pre_sum, primary.query_alignment_end, len(primary.query_sequence), 'T', pm_seq[primary.query_alignment_start+pre_sum: primary.query_alignment_end], "AP", filter_flag, mapq, align_mismatch_ration_pr)
                    read_aligns.append(align_obj1)
                    read_aligns.append(align_obj2)
                    open_flag = 1
                    break
                pre_sum += int(op[0])
            if open_flag == 0:
                continue
        else: 
            pm_sa = [primary] + supplementaries

            primary_is_reverse = primary.is_reverse
            read_name = primary.query_name

            pm_seq = primary.query_sequence

            # # traverse all aligns
            for align in pm_sa:
                align_chr = bam_file.getrname(align.reference_id)
                if align_chr not in options.all_possible_chrs:
                    continue
                mapq = align.mapping_quality
                if options.mismatch_filter:
                    if align.is_supplementary:
                        align_mismatch_ration = align.get_tag("NM") / (align.reference_end - align.reference_start)
                    else:
                        align_mismatch_ration = get_mismatch_ratio(align.cigarstring, align.reference_start, align.reference_end)
                    # if align_mismatch_ration >= options.max_mismatch_ratio:
                    #     filter_flag = 'FILTER'
                # # focus on primary's forward, if other reads' forward is not same as primary, then adjust cords
                if align.is_reverse != primary_is_reverse:
                    align_q_start = align.query_length - align.query_alignment_end
                    align_q_end = align.query_length - align.query_alignment_start
                else:
                    align_q_start = align.query_alignment_start
                    align_q_end = align.query_alignment_end

                # reset seg's forward to the primary's forward
                # and those whose forward is the same as the primary, set reverse tag to False
                if align.is_reverse == primary_is_reverse:
                    align_direct = 'T'
                else:
                    align_direct = 'F'

                if align.reference_end - align.reference_start >= options.min_align_length:    
                    # # store the read's info in a dict
                    align_obj = Align(read_name, align_chr, align.reference_start, align.reference_end, align_q_start, align_q_end, len(align.query_sequence), align_direct, pm_seq[align_q_start: align_q_end], 'SA' if align.is_supplementary else "PM", filter_flag, mapq, align_mismatch_ration)
                    read_aligns.append(align_obj)
        
        # no align or one align means no bkp, then continue
        if len(read_aligns) <= 1:
            continue

        # sort by read and assign index
        srt_read_aligns = sorted(read_aligns, key=lambda aln: (aln.read_start, aln.read_end))

        if judge_noisy_aligns(srt_read_aligns):
            continue

        filtered_read_aligns = []
        for read_align in srt_read_aligns:
            # if read_align.chr != None:
            filtered_read_aligns.append(read_align)
        
        # if len(filtered_read_aligns) > 4:
        #     continue

        mapq_flag = 0
        for read in filtered_read_aligns:
            if read.mapq > 40:
                mapq_flag = 1
                break
        if mapq_flag == 0:
            continue
        for read in filtered_read_aligns:
            reads_mismatch_ratio_file.write("{};".format(read.mismatch_ratio))

        segments_direct = []
        for i in range(len(filtered_read_aligns)):
            filtered_read_aligns[i].set_index(i)
            segments_direct.append(filtered_read_aligns[i].direct)   
            
        # # STEP: detect from aligns
        included_juncs, end = detect_intra_inter_trans(filtered_read_aligns, options, repeat_regions, pm_seq, non_local_thresh=options.non_local_thresh)

        if len(included_juncs) == 0:
            continue

        all_juncs_support_reads = []
        for junc in included_juncs:
            all_juncs_support_reads.append(junc.support_reads)

        origin_juncSe_list.append(Junction_Series(included_juncs, all_juncs_support_reads, end))
    
    clu_juncSe_list = cluster_juncSe_list(origin_juncSe_list, 0)
    srt_origin_juncSe_list = sorted(clu_juncSe_list,key=lambda x:x.get_junc_num(),reverse=True)
    juncSe_list = cluster_juncSe_list(srt_origin_juncSe_list, 0)

    if options.seq_type == 'hifi':
        min_support = 5
    elif options.seq_type == 'ont':
        min_support = 10

    filtered_juncSe_list = []
    for juncSe in juncSe_list:
        if juncSe.get_support_num() >= min_support:
            filtered_juncSe_list.append(juncSe)

    # filtered_juncSe_list = []
    # for juncSe in juncSe_list:
    #     if juncSe.get_support_num() >= options.min_support:
    #         if len(juncSe.include_juncs) == 1:
    #             filtered_juncSe_list.append(juncSe)
    #             continue
    #         if len(juncSe.include_juncs) > 1:
    #             if len(juncSe.include_juncs[0].support_reads) < 2:
    #                 flag1 = 0
    #                 for i in range(1, len(juncSe.include_juncs)):
    #                     if len(juncSe.include_juncs[i].support_reads) >= 2:
    #                         juncSe.end[0] = juncSe.include_juncs[i-1].include_bkps[1]
    #                         juncSe.include_juncs = juncSe.include_juncs[i:]
    #                         filtered_juncSe_list.append(juncSe)
    #                         flag1 = 1
    #                         break
    #                 if flag1 == 0:
    #                     break
    #             if len(juncSe.include_juncs[0].support_reads) >= 2:
    #                 flag2 = 0
    #                 for i in range(1, len(juncSe.include_juncs)):
    #                     if len(juncSe.include_juncs[i].support_reads) < 2:
    #                         juncSe.end[1] = juncSe.include_juncs[i].include_bkps[0]
    #                         juncSe.include_juncs = juncSe.include_juncs[:i]
    #                         filtered_juncSe_list.append(juncSe)
    #                         flag2 = 1
    #                         break
    #                 if flag2 == 0:
    #                     filtered_juncSe_list.append(juncSe)
    
    if not os.path.exists(options.out_path+'/juncSe_data'):
        os.mkdir(options.out_path+'/juncSe_data')
    with open(options.out_path + '/juncSe_data'+'/juncSe_{}_{}_{}.pkl'.format(chrom, region_start, region_end), 'wb') as f:
            pickle.dump(filtered_juncSe_list, f)

    return filtered_juncSe_list