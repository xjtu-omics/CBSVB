import datetime
import os
import pysam

from src.classes import Align, Junction_Series
from src.cigar_op import cigar_to_list, cal_align_read_start_and_length_from_cigar
from src.detect_op import detect_intra_inter_trans, judge_noisy_aligns
from src.cluster_op import cluster_juncSe_list
from src.other_op import filter_juncs_trfregion, get_mismatch_ratio
from src.load_rmsk_trf_tree import load_interval_trees_from_json
import concurrent.futures
from intervaltree import Interval
import pickle
import logging

def traverse_reads_ngs(bam_path, chrom, region_start, region_end, repeat_regions, mode, options):

    """
    producer, used to traverse bam and collect read1 and read2 and their SAs to save 
    :param bam_path:
    :param trans_dis_threshold:
    :return:
    """
    start_time = datetime.datetime.now()

    print("[Processing {}:{}-{}]\n".format(chrom, region_start, region_end))

    # # STEP: collect read's read1 and read2
    # region_str = "{}-{}-{}".format(chrom, region_start, region_end)
    rmsk_trf_tree = None
    if os.path.exists('./data/rmsk_trf_tree/{}_{}_{}.json'.format(chrom, region_start, region_end)):
        rmsk_trf_tree = load_interval_trees_from_json('./data/rmsk_trf_tree/{}_{}_{}.json'.format(chrom, region_start, region_end))
    # if chrom == None and region_start == None and region_end == None:
    #     rmsk_intervals = load_interval_trees_from_json('./data/rmsk_trf_tree.json')
    #     print("got tree")

    read_pair_collection = {}
    bam_file = pysam.AlignmentFile(bam_path)

    if mode == "mated":
        aligns = bam_file.fetch(chrom, region_start, region_end)
        unmated_read_prefix = bam_path.split("/")[-1].replace(".bam", ".unmated_{}_{}_{}.bam".format(chrom, region_start, region_end))
        unmated_read_bam = os.path.join(options.out_path, unmated_read_prefix)
        unmated_read_bam_out = pysam.AlignmentFile(unmated_read_bam, "wb", header=bam_file.header)

    else:
        aligns = bam_file
        unmated_read_bam_out = None
        unmated_read_bam = None

    for pm_align in aligns:

        # # unmapped or secondary or low mapping quality, then pass this align
        if pm_align.cigarstring == None or pm_align.is_unmapped or pm_align.is_secondary:
            continue

        if pm_align.mapping_quality < options.min_mapq:
            continue

        # # skip supplementary align
        if pm_align.is_supplementary:
            continue
        if pm_align.mate_is_unmapped:
            continue

        if rmsk_trf_tree != None and pm_align.reference_name in rmsk_trf_tree:
            query_interval = Interval(pm_align.reference_start, pm_align.reference_end)
            # 使用 overlaps 方法检查重叠
            overlaps = rmsk_trf_tree[chrom].overlap(query_interval)
            if len(overlaps) > 50:
                continue
        
        # if chrom == None and region_start == None and region_end == None:
        #     overlapping_intervals = []
        #     # 如果染色体存在于树中
        #     if pm_align.reference_name in rmsk_intervals:
        #         tree = rmsk_intervals[pm_align.reference_name]
        #         # 查询重叠的区间
        #         overlapping_intervals = tree.overlap(pm_align.reference_start, pm_align.reference_end)
        #         if len(overlapping_intervals) > 0:
        #             overlap_start = max(list(overlapping_intervals)[0].begin, pm_align.reference_start)
        #             overlap_end = min(list(overlapping_intervals)[0].end, pm_align.reference_end)
        #             overlap_length = overlap_end - overlap_start
        #             if overlap_length > 50:
        #                 continue

        # if options.mismatch_filter:
        #     align_mismatch_ration = get_mismatch_ratio(pm_align.cigarstring, pm_align.reference_start, pm_align.reference_end)
        #     if align_mismatch_ration >= options.max_mismatch_ratio:
        #         continue

        read_name = pm_align.qname

        # # add to dict
        if read_name not in read_pair_collection:
            read_pair_collection[read_name] = [None, None]

        if pm_align.is_read1:
            read_pair_collection[read_name][0] = pm_align
        else:
            read_pair_collection[read_name][1] = pm_align

    partime_time = datetime.datetime.now()
    logging.info("collect cost {}s. read num {}".format((partime_time - start_time).seconds, len(read_pair_collection.keys())))

    # # STEP: traverse each read to get aligns and detect trans
    origin_juncSe_list = []

    not_mate_num = 0
    for read_name in read_pair_collection:
        read1_pm_align = read_pair_collection[read_name][0]
        read2_pm_align = read_pair_collection[read_name][1]
        
        if read1_pm_align is None:
            not_mate_num += 1
            if mode == "mated":
                unmated_read_bam_out.write(read2_pm_align)
            continue

        if read2_pm_align is None:
            not_mate_num += 1
            if mode == "mated":
                unmated_read_bam_out.write(read1_pm_align)
            continue
        
        read_aligns = []

        read1_chr = 0
        read2_chr = 0
        if read1_pm_align.reference_name != read2_pm_align.reference_name:
            if read1_pm_align.reference_name=='chrX':
                read1_chr = 23
            if read2_pm_align.reference_name=='chrX':
                read2_chr = 23
            if read1_chr == 0:
                read1_chr = int(read1_pm_align.reference_name[3:])
            if read2_chr == 0:
                read2_chr = int(read2_pm_align.reference_name[3:]) 
            if read1_chr > read2_chr:
                flag_2_1000 = 0

        # # dicide the order of the reads (flag_2_1000=1 means read2 in the back)
        flag_2_1000 = 1
        # put the read on the negative chain in the back
        if read1_pm_align.is_reverse != read2_pm_align.is_reverse and read1_pm_align.is_reverse == True:
            flag_2_1000 = 0
        if read1_pm_align.is_reverse == read2_pm_align.is_reverse and read1_pm_align.reference_name == read2_pm_align.reference_name and read1_pm_align.reference_start > read2_pm_align.reference_start:
            flag_2_1000 = 0
        if read1_pm_align.is_reverse == read2_pm_align.is_reverse and read1_pm_align.is_reverse == False and read1_chr > read2_chr:
            flag_2_1000 = 0
        if read1_pm_align.is_reverse == read2_pm_align.is_reverse and read1_pm_align.is_reverse == True and read1_chr < read2_chr:
            flag_2_1000 = 0

        for pm_align in [read1_pm_align, read2_pm_align]:
            pm_direct = 'F' if pm_align.is_reverse else 'T'
            # # get supplementarys
            try:
                sa_tags = pm_align.get_tag("SA").split(";")[0]
            except KeyError:
                sa_tags = ''

            # # deal supplementarys and push them to dict to save
            if len(sa_tags) != 0:
                sa_info_split = sa_tags.split(',')
                sa_chr = sa_info_split[0]
                sa_ref_start = int(sa_info_split[1])
                sa_direct = 'T' if sa_info_split[2] == '+' else 'F'
                sa_cigar_ops, sa_cigar_lens = cigar_to_list(sa_info_split[3])
                # calculate sa read start and length by cigar and direct
                sa_read_start, sa_length = cal_align_read_start_and_length_from_cigar(sa_cigar_ops, sa_cigar_lens, pm_direct, sa_direct)

                sa_ref_end = sa_ref_start + sa_length
                sa_read_end = sa_read_start + sa_length
                
                # push to dict to save
                if pm_align.is_read1:
                    if flag_2_1000 == 1:
                        read_aligns.append(Align(read_name, sa_chr, sa_ref_start, sa_ref_end, sa_read_start, sa_read_end, sa_length, sa_direct, pm_align.seq[sa_read_start: sa_read_end], 'read1_sa', read1_pm_align.isize))
                    elif flag_2_1000 == 0:
                        read_aligns.append(Align(read_name, sa_chr, sa_ref_start, sa_ref_end, 1000 + sa_read_start, 1000 + sa_read_end, sa_length, sa_direct, pm_align.seq[sa_read_start: sa_read_end], 'read1_sa', read1_pm_align.isize))
                else:
                    if flag_2_1000 == 1:
                        read_aligns.append(Align(read_name, sa_chr, sa_ref_start, sa_ref_end, 1000 + sa_read_start, 1000 + sa_read_end, sa_length, sa_direct, pm_align.seq[sa_read_start: sa_read_end], 'read2_sa', read1_pm_align.isize))
                    elif flag_2_1000 == 0:
                        read_aligns.append(Align(read_name, sa_chr, sa_ref_start, sa_ref_end, sa_read_start, sa_read_end, sa_length, sa_direct, pm_align.seq[sa_read_start: sa_read_end], 'read2_sa', read1_pm_align.isize))
           
            pm_chr = pm_align.reference_name
            pm_ref_start = pm_align.reference_start

            pm_direct = 'F' if pm_align.is_reverse else 'T'
            pm_cigar_ops, pm_cigar_lens = cigar_to_list(pm_align.cigarstring)

            # calculate align's read start by find the first 'M' in CIGAR
            pm_read_start, pm_length = cal_align_read_start_and_length_from_cigar(pm_cigar_ops, pm_cigar_lens, pm_direct, pm_direct)
            pm_seq = pm_align.query_sequence[pm_read_start: pm_read_start + pm_length]

            # calculate read and ref end
            pm_ref_end = pm_ref_start + pm_length
            pm_read_end = pm_read_start + pm_length
            # read1
            if pm_align.is_read1:
                if flag_2_1000 == 1:
                    read_aligns.append(Align(read_name, pm_chr, pm_ref_start, pm_ref_end, pm_read_start, pm_read_end, pm_length, pm_direct, pm_seq, 'read1', pm_align.isize))
                elif flag_2_1000 == 0:
                    read_aligns.append(Align(read_name, pm_chr, pm_ref_start, pm_ref_end, 1000 + pm_read_start, 1000 + pm_read_end, pm_length, pm_direct, pm_seq, 'read1', pm_align.isize))
            else:
                if flag_2_1000 == 1:
                    read_aligns.append(Align(read_name, pm_chr, pm_ref_start, pm_ref_end, 1000 + pm_read_start, 1000 + pm_read_end, pm_length, pm_direct, pm_seq, 'read2', pm_align.isize))
                elif flag_2_1000 == 0:
                    read_aligns.append(Align(read_name, pm_chr, pm_ref_start, pm_ref_end, pm_read_start, pm_read_end, pm_length, pm_direct, pm_seq, 'read2', pm_align.isize))

        # sort by read and assign index
        srt_read_aligns = sorted(read_aligns, key=lambda aln: (aln.read_start, aln.read_end))

        for i in range(len(srt_read_aligns)):
            srt_read_aligns[i].set_index(i)

        if judge_noisy_aligns(srt_read_aligns):
            continue

        # change the direct of the reads in the back
        for align in srt_read_aligns:
            if align.align_type == 'read1' or align.align_type == 'read2':
                read_type = align.align_type[4]
                break
        for align in srt_read_aligns:
            if align.align_type[4] != read_type:
                align.direct = 'T' if align.direct == 'F' else 'F'
            
        # # for inv read with sa, change the order
        if len(srt_read_aligns)>2:
            for i in range(len(srt_read_aligns)):
                if srt_read_aligns[i].direct == 'F':
                    if srt_read_aligns[i].align_type == 'read1':
                        if i<len(srt_read_aligns)-1 and srt_read_aligns[i+1].align_type == 'read1_sa':
                            srt_read_aligns[i], srt_read_aligns[i+1] = srt_read_aligns[i+1], srt_read_aligns[i]
                            break
                        elif i>0 and srt_read_aligns[i-1].align_type == 'read1_sa':
                            srt_read_aligns[i], srt_read_aligns[i-1] = srt_read_aligns[i-1], srt_read_aligns[i]
                            break
                    elif srt_read_aligns[i].align_type == 'read2':
                        if i<len(srt_read_aligns)-1 and srt_read_aligns[i+1].align_type == 'read2_sa':
                            srt_read_aligns[i], srt_read_aligns[i+1] = srt_read_aligns[i+1], srt_read_aligns[i]
                            break
                        elif i>0 and srt_read_aligns[i-1].align_type == 'read2_sa':
                            srt_read_aligns[i], srt_read_aligns[i-1] = srt_read_aligns[i-1], srt_read_aligns[i]
                            break
                    else:
                        continue
                    

        segments_direct = []
        for align in srt_read_aligns:
            segments_direct.append(align.direct)

        # detect trans and write to file  
        included_juncs, end = detect_intra_inter_trans(srt_read_aligns, options, repeat_regions, non_local_thresh=options.non_local_thresh)

        if len(included_juncs)>2 or len(included_juncs)==0:
            continue

        support_reads = []
        for junc in included_juncs:
            support_reads.append(junc.support_reads)
        origin_juncSe_list.append(Junction_Series(included_juncs, support_reads, end))
    
    partime_time2 = datetime.datetime.now()
    logging.info("detect cost {}s. detect num {}. unmated {}".format((partime_time2 - partime_time).seconds, len(origin_juncSe_list), not_mate_num))

    print("detect cost {}s. detect num {}. unmated {}".format((partime_time2 - partime_time).seconds, len(origin_juncSe_list), not_mate_num))

    if chrom != None:
        clusted_juncSe_list = cluster_juncSe_list(origin_juncSe_list, 0)
        srt_juncSe_list = sorted(clusted_juncSe_list, key=lambda x:x.get_junc_num(),reverse=True)
        juncSe_list = cluster_juncSe_list(srt_juncSe_list, 0)

        # filter junSe_list by support
        filtered_juncSe_list = []
        for juncSe in juncSe_list:
            if juncSe.get_support_num() >= options.min_support:
                filtered_juncSe_list.append(juncSe)
    else:
        filtered_juncSe_list = origin_juncSe_list

    # filtered_juncSe_list = []
    # for juncSe in juncSe_list:
    #     if juncSe.get_support_num() >= 3:
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
    
    # partime_time3 = datetime.datetime.now()    
    # logging.info("cluster cost {}s. after cluster {}, {}".format((partime_time3 - partime_time2).seconds, len(juncSe_list), len(filtered_juncSe_list)))
    
    # end_time = datetime.datetime.now()
    # logging.info("End cost {}s. {}".format((end_time - start_time).seconds, len(filtered_juncSe_list)))

    if mode == "mated":
        unmated_read_bam_out.close()
    logging.info("end")

    if not os.path.exists(options.out_path+'/juncSe_data'):
        os.mkdir(options.out_path+'/juncSe_data')
    if not os.path.exists(options.out_path+'/unmated_data'):
        os.mkdir(options.out_path+'/unmated_data')
    with open(options.out_path+'/juncSe_data'+'/juncSe_{}_{}_{}.pkl'.format(chrom, region_start, region_end), 'wb') as f:
        pickle.dump(filtered_juncSe_list, f)
    with open(options.out_path+'/unmated_data'+'/unmated_{}_{}_{}.pkl'.format(chrom, region_start, region_end), 'wb') as f:
        pickle.dump(unmated_read_bam, f)

    # return fitered_juncSe_list, single_same_list, unmated_read_bam
    return filtered_juncSe_list, unmated_read_bam
