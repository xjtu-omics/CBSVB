from collections import Counter
import os
import numpy as np

from src.classes_another import Junction_Series

from intervaltree import Interval, IntervalTree
import json

def load_interval_trees_from_json(filename):
    with open(filename, 'r') as json_file:
        tree_data = json.load(json_file)
    
    rmsk_trees = {}
    for chrom, intervals in tree_data.items():
        tree = IntervalTree()
        for start, end in intervals:
            tree.add(Interval(start, end))
        rmsk_trees[chrom] = tree
    
    return rmsk_trees

def filter_juncs_trfregion(include_juncs, end_info, rmsk_trf, data_type):
    if len(include_juncs) == 0:
        return False
    if data_type == 'tgs':
        judge_number = 1
    else:
        judge_number = 0
    repeat_region_num = 0
    for i in range(len(include_juncs)+1):
        if i==0:
            chrom = end_info[0].chrom
            pos1 = end_info[0].pos
            pos2 = include_juncs[0].include_bkps[0].pos
        elif i==len(include_juncs):
            chrom = end_info[1].chrom
            pos1 = include_juncs[len(include_juncs)-1].include_bkps[1].pos
            pos2 = end_info[1].pos
        elif i < len(include_juncs)-1:
            chrom = include_juncs[i].include_bkps[1].chrom
            pos1 = include_juncs[i].include_bkps[1].pos
            pos2 = include_juncs[i+1].include_bkps[0].pos
        else:
            continue
        start = min(pos1, pos2)
        end = max(pos1, pos2)
        start_flag = 0
        end_flag = 0
        for info in rmsk_trf[chrom]:
            if abs(info[0]-start)<30:
                start_flag = 1
            if abs(info[1]-end)<30:
                end_flag = 1
            if start_flag == 1 and end_flag == 1:
                repeat_region_num += 1
                break
            if info[1] - end > 300:
                break
    if repeat_region_num > judge_number:
        return False
    return True

def get_regions(options, mode):
    regions = dict()
    if mode == 'repeat_regions':
        file_path = options.repeat_region_path
    elif mode == 'rmsk_trf':
        file_path = options.rmsk_trf_path
    else:
        ValueError("Illegible mode")
    if file_path == None:
        return None
    with open(file_path) as region_file:
        for region in region_file:
            region_info = region.split('\n')[0].split('\t')
            if region_info[0] in options.all_possible_chrs:
                if region_info[0] not in regions.keys():
                    regions[region_info[0]] = []
                regions[region_info[0]].append((int(region_info[1]), int(region_info[2]), region_info[3]))
            else:
                ValueError("Illegible chromosome: {}".format(region_info[0]))
    region_file.close()
    return regions

def get_mismatch_ratio(cigarstring, reference_start, reference_end):
    op_counter = Counter(cigarstring)

    if 'X' not in op_counter:
        op_counter['X'] = 0
    
    if 'I' not in op_counter:
        op_counter['I'] = 0

    if 'D' not in op_counter:
        op_counter['D'] = 0
    
    mismatch_ratio = (op_counter['X']+op_counter['I']+op_counter['D']) / (reference_end - reference_start)
    
    return mismatch_ratio

def report_result(sample_name, juncSe_list, options):
    with open(os.path.join(options.out_path, "{}.txt".format(sample_name)), "w") as output_file:
        output_file.write('#id\tjunc_num\tsupport_reads\tdescription\n')
        for i in range(len(juncSe_list)):
            if i < len(juncSe_list)-1:
                output_file.write('Trans.{}\t{}\t{}\t{}\n'.format(i, juncSe_list[i].get_junc_num(), juncSe_list[i].get_support_num(), juncSe_list[i].to_string()))
            else:
                output_file.write('Trans.{}\t{}\t{}\t{}'.format(i, juncSe_list[i].get_junc_num(), juncSe_list[i].get_support_num(), juncSe_list[i].to_string()))
    output_file.close()

def count_pass(junc):
    flag = 0
    cnt_pass_reads = 0
    for reads in junc.support_reads:
        cnt_pass = 0
        for read in reads:
            if read.filter_flag == 'PASS':
                cnt_pass += 1
        if cnt_pass/len(reads) >= 0.5:
            cnt_pass_reads += 1
    if cnt_pass_reads/len(junc.support_reads) > 0.5:
        flag = 1
    return flag

def get_reliable_juncSe(srt_juncSe_list, threshold, seq_type):
    # if os.path.exists('./data/exp.rmsk.json'):
    #     rmsk_intervals = load_interval_trees_from_json('./data/exp.rmsk.json')
    #     print("got_rmsk_tree")
    # else:
    rmsk_intervals = None
    print("Warning: No rmsk_tree, run with no filter")
    # if os.path.exists('./data/region_tree_{}_chr1.json'.format(seq_type)):
    #     low_high_depth_intervals = load_interval_trees_from_json('./data/region_tree_{}_chr1.json'.format(seq_type))
    #     print("got_depth_tree")
    # else:
    low_high_depth_intervals = None
    print("Warning: No rmsk_tree, run with no filter")

    for juncSe in srt_juncSe_list:
        for junc in juncSe.include_juncs:
            for reads in junc.support_reads:
                for read in reads:
                    if read.mismatch_ratio > threshold:
                        read.filter_flag = 'FILTER'
        juncSe.update_support_reads()
    
    re_juncSe_list = []
    for juncSe in srt_juncSe_list:
        if len(juncSe.include_juncs) == 1:
            re_juncSe_list.append(juncSe)
            continue
        reliable_junc = []
        for i in range(len(juncSe.include_juncs)):
            if count_pass(juncSe.include_juncs[i]) and no_filter_region(juncSe.include_juncs[i], low_high_depth_intervals) and \
                no_filter_region(juncSe.include_juncs[i], rmsk_intervals):
                reliable_junc.append(i)
        if len(reliable_junc) == 1:
            if reliable_junc[0] == 0:
                end0 = juncSe.end[0]
            else:
                end0 = juncSe.include_juncs[reliable_junc[0]-1].include_bkps[1]
            if reliable_junc[0] == len(juncSe.include_juncs)-1:
                end1 = juncSe.end[1]
            else:
                end1 = juncSe.include_juncs[reliable_junc[0]+1].include_bkps[0]
            re_juncSe_list.append(Junction_Series(juncSe.include_juncs[reliable_junc[0]:reliable_junc[0]+1], juncSe.support_reads[reliable_junc[0]:reliable_junc[0]+1], [end0, end1]))
        else:
            last_index = 0
            for i in range(len(reliable_junc)-1):
                if reliable_junc[i+1]-reliable_junc[i]>1:
                    if last_index == 0:
                        end0 = juncSe.end[0]
                    else:
                        end0 = juncSe.include_juncs[last_index].include_bkps[1]
                    end1 = juncSe.include_juncs[reliable_junc[i]+1].include_bkps[0]
                    end = [end0, end1]
                    re_juncSe_list.append(Junction_Series(juncSe.include_juncs[last_index:reliable_junc[i]+1], juncSe.support_reads[last_index:reliable_junc[i]+1], end))
                    last_index = reliable_junc[i+1]-1
                    if i+1 == len(reliable_junc)-1:
                        end0 = juncSe.include_juncs[reliable_junc[i]].include_bkps[1]
                        end1 = juncSe.end[1]
                        end = [end0, end1]
                        re_juncSe_list.append(Junction_Series(juncSe.include_juncs[reliable_junc[i+1]:reliable_junc[i+1]+1], juncSe.support_reads[reliable_junc[i+1]:reliable_junc[i+1]+1], end))
                elif i==len(reliable_junc)-2 and reliable_junc[i+1]-reliable_junc[i] == 1:
                    if last_index == 0:
                        end0 = juncSe.end[0]
                    else:
                        end0 = juncSe.include_juncs[last_index].include_bkps[1]
                    end1 = juncSe.end[1]
                    end = [end0, end1]
                    re_juncSe_list.append(Junction_Series(juncSe.include_juncs[last_index:len(juncSe.include_juncs)], juncSe.support_reads[last_index:len(juncSe.include_juncs)], end))
    
    return re_juncSe_list

def get_mismatch_ratio_threshold(mismatch_ratio_data_path):
    all_mismatch_ratios = []
    files = os.listdir(mismatch_ratio_data_path)
    for file in files:
        f = open(mismatch_ratio_data_path+'/'+file, 'r')
        for line in f:
            mismatch_ratios = line.split(';')
            all_mismatch_ratios.extend([float(ratio) for ratio in mismatch_ratios if ratio != ""])

    # percentile_95 = np.percentile(all_mismatch_ratios, 95)
    # print("threshold:{}".format(percentile_95))
    # return percentile_95
    percentile_90 = np.percentile(all_mismatch_ratios, 90)
    print("threshold:{}".format(percentile_90))
    return percentile_90

def no_filter_region(junc, region_intervals):
    if region_intervals == None:
        return True
    for i in range(2):
        chrom = junc.include_bkps[i].chrom
        pos1 = junc.include_bkps[i].pos
        if junc.include_bkps[i].direct == 'S':
            pos2 = pos1 - 1000
        else:
            pos2 = pos1 + 1000
        if region_intervals != None:
            interval_start = min(pos1, pos2)
            interval_end = max(pos1, pos2)
            overlapping_intervals = []
            # 如果染色体存在于树中
            if chrom in region_intervals:
                tree = region_intervals[chrom]
                # 查询重叠的区间
                overlapping_intervals = tree.overlap(interval_start, interval_end)
                if len(overlapping_intervals) > 0:
                    overlap_start = max(list(overlapping_intervals)[0].begin, interval_start)
                    overlap_end = min(list(overlapping_intervals)[0].end, interval_end)
                    overlap_length = overlap_end - overlap_start
                    if overlap_length > 100:
                        return False
    return True

def get_min_junc_support(junc, type):
    depth1 = 0
    depth2 = 0
    chrom1 = junc.include_bkps[0].chrom
    pos1 = junc.include_bkps[0].pos
    chrom2 = junc.include_bkps[1].chrom
    pos2 = junc.include_bkps[1].pos
    files = os.listdir('./data/depths/')
    for file in files:
        if file.split('.')[0].split('_')[3] != type:
            continue
        if file.split('.')[0].split('_')[0] == chrom1:
            if int(file.split('.')[0].split('_')[1]) < pos1 < int(file.split('.')[0].split('_')[2]):
                with open('./data/depths/'+file, 'r') as json_file:
                    depth = json.load(json_file)
                depth1 = np.mean(np.array(depth[pos1-int(file.split('.')[0].split('_')[1])-500:pos1-int(file.split('.')[0].split('_')[1])+500]))
        if file.split('.')[0].split('_')[0] == chrom2:
            if int(file.split('.')[0].split('_')[1]) < pos2 < int(file.split('.')[0].split('_')[2]):
                with open('./data/depths/'+file, 'r') as json_file:
                    depth = json.load(json_file)
                depth2 = np.mean(np.array(depth[pos2-int(file.split('.')[0].split('_')[1])-500:pos2-int(file.split('.')[0].split('_')[1])+500]))
    min_support = (depth1 + depth2) *0.5 * 0.1 # depth的1/10作为最小support reads
    return min_support

def get_min_multijunc_support(juncs, type):
    depth = 0
    chrom = juncs.same_bkp_info.chrom
    pos = juncs.same_bkp_info.pos
    files = os.listdir('./data/depths/')
    for file in files:
        if file.split('.')[0].split('_')[3] != type:
            continue
        if file.split('.')[0].split('_')[0] == chrom:
            if int(file.split('.')[0].split('_')[1]) < pos < int(file.split('.')[0].split('_')[2]):
                with open('./data/depths/'+file, 'r') as json_file:
                    depth = json.load(json_file)
                depth = np.mean(np.array(depth[pos-int(file.split('.')[0].split('_')[1])-500:pos-int(file.split('.')[0].split('_')[1])+500]))
    min_support = depth * 0.1
    return min_support
    