
import os, sys, time, pysam, editdistance
from subprocess import Popen, PIPE
from spoa import poa
from multiprocessing import Pool

def input_organize(ref_fn, output_dir, input_fn):

    index_info = 'result,NNNNNNNNNN'

    ref_f = open(ref_fn).readlines()
    ref_name = ref_f[0].split()[0][1:]
    ref_seq = ''.join(ref_f[1:]).replace('\n', '').upper()
    half_pos = int(len(ref_seq)/2)
    
    umi_len = 0
    
    if ref_seq[:half_pos].find('(') != -1:
        front_pos = ref_seq.find('(')
        end_pos = ref_seq.find(')')
        forward_umi_pos = [front_pos, end_pos - 2]
        ref_seq = ref_seq[:front_pos] + ref_seq[front_pos + 1: end_pos] + ref_seq[end_pos + 1:]
        umi_len += end_pos - 2 - front_pos + 1
    else:
        forward_umi_pos = [0, -1]
    
    if ref_seq[half_pos:].find('(') != -1:
        front_pos = ref_seq.rfind('(')
        end_pos = ref_seq.rfind(')')
        back_umi_pos = [front_pos, end_pos - 2]
        ref_seq = ref_seq[:front_pos] + ref_seq[front_pos + 1: end_pos] + ref_seq[end_pos + 1:]
        umi_len += end_pos - 2 - front_pos + 1
    else:
        back_umi_pos = [0, -1]
    
    umi_pos = [forward_umi_pos, back_umi_pos]
    if umi_pos[0][1] == -1 and umi_pos[1][1] == -1:
        print('ERROR: There is any UMI sequence in reference, please check the reference and primer format.')
        sys.exit()
    if umi_pos[0][1] != -1:
        print(f'UMI position in the front part of sequence is from {umi_pos[0][0] + 1} to {umi_pos[0][1] + 1}.')
    else:
        print(f'There is no UMI at the front part of the sequence.')
    if umi_pos[1][1] != -1:
        print(f'UMI position in the back part of the sequence is from {umi_pos[1][0] + 1} to {umi_pos[1][1] + 1}.')
    else:
        print(f'There is no UMI at the back part of the sequence.')

    return umi_pos, umi_len, ref_name, ref_seq


# Extract UMI sequence

def extract_index_umi(sam_file, ref_file, output_dir, umi_front_pos, umi_back_pos, index_information='result,NNNNNNNNNN', umi_idx_r=1, idx_front_pos='0,-1', idx_back_pos=None):
    
    def cigar_len(cigar):
        l = 0
        for i in cigar:
            if i[0] in [0,2]:
                l += i[1]
        return l

    def rc(seq):
        return seq.translate(seq.maketrans('ATGCatgc', 'TACGtacg'))[::-1]

    def align_info(cigar, seq, align_st, align_len):
        seq_align = ['-' for i in range(align_st)]
        pos = 0
        for c in cigar:
            if c[0] == 0:
                for x in range(c[1]):
                    seq_align.append(seq[pos + x])
                if len(seq_align) > align_len:
                    return seq_align
                pos += c[1]
            elif c[0] == 1:
                seq_align[-1] += seq[pos: pos + c[1]]
                pos += c[1]
            elif c[0] == 2:
                for x in range(c[1]):
                    seq_align.append('-')
                if len(seq_align) > align_len:
                    return seq_align
            elif c[0] == 4:
                pos += c[1]

    def cut_align(seq, seq_qual, cigar): # remove soft clipping
        if cigar[0][0] == 4:
            seq = seq[cigar[0][1]:]
            seq_qual = seq_qual[cigar[0][1]:]
        if cigar[-1][0] == 4:
            seq = seq[:-1 * cigar[-1][1]]
            seq_qual = seq_qual[: -1 * cigar[-1][1]]
        return seq, seq_qual

    def demultiplex_idx(idx_seq, idx_dict): #find idx using editdistance
        scr = {}
        for i, idx_s in idx_dict.items():
            scr[i] = editdistance.eval(idx_s, idx_seq)
        val = len(idx_seq)
        multi_check = False
        fid = False
        for i, s in scr.items():
            if val > s:
                val = s
                fid = i
                multi_check = False
            elif val == s:
                multi_check = True
        if multi_check == True:
            fid = False
        return fid

    def extract_seq(align_seq, cut_st, cut_ed, r=15, cutoff=0.8):
        if align_seq[cut_st - r: cut_st].count('-') > r * cutoff:
            return ''
        elif align_seq[cut_ed: cut_ed + r].count('-') > r * cutoff:
            return ''
        return ''.join(align_seq[cut_st: cut_ed]).replace('-', '')

    ref_seq = ''.join(open(ref_file).readlines()[1:]).replace('\n','')
    ref_seq_len = len(ref_seq)
    
    idx_info = index_information.split(',')
    idx_pos = [[], []]
    idx_front = idx_front_pos
    for i in idx_front.split(','):
        idx_pos[0].append(int(i))
    idx_back = idx_back_pos
    if idx_back == None:
        idx_pos[1] = [ref_seq_len, ref_seq_len - 1]
    else:
        for i in idx_back.split(','):
            idx_pos[1].append(int(i))

    umi_pos = [[], []]
    umi_front = umi_front_pos
    for i in umi_front:
        umi_pos[0].append(int(i))
    umi_back = umi_back_pos
    if umi_back == None:
        umi_pos[1] = [ref_seq_len, ref_seq_len - 1]
    else:
        for i in umi_back:
            umi_pos[1].append(int(i))
    
    idx_dict = {}
    for i in range(int(len(idx_info)/2)):
        idx_dict[idx_info[2*i]] = idx_info[2*i+1]

    cnt_dict = {'all': 0, 'secondary': 0, 'unmap': 0,'filted_short': 0, 'filted_idx_len': 0, 'pass': 0, 'split_pass': 0} # "All" means the whole number of Aligned reads

    with pysam.AlignmentFile(sam_file) as f:

        if umi_pos[0][1] < idx_pos[0][1]:
            front_align_cut = idx_pos[0][1]
        else:
            front_align_cut = umi_pos[0][1]

        if umi_pos[1][0] < idx_pos[1][0]:
            back_align_cut = umi_pos[1][0]
        else:
            back_align_cut = idx_pos[1][0]
        
        if idx_pos[0][0] == 0 and umi_pos[0][0] == 0:
            front_indi_pos = 0
        elif idx_pos[0][0] != 0 and umi_pos[0][0] != 0:
            front_indi_pos = min([idx_pos[0][0], umi_pos[0][0]])
        elif umi_pos[0][0] == 0:
            front_indi_pos = idx_pos[0][0]
        else:
            front_indi_pos = umi_pos[0][0]

        if idx_pos[1][1] == ref_seq_len and umi_pos[1][1] == ref_seq_len:
            back_indi_pos = ref_seq_len
        elif idx_pos[1][1] != 0 and umi_pos[1][1] != 0:
            back_indi_pos = min([idx_pos[1][1], umi_pos[1][1]])
        elif umi_pos[1][1] == ref_seq_len:
            back_indi_pos = idx_pos[1][1]
        else:
            back_indi_pos = umi_pos[1][1]
        
        demulti_fw_dict = {}
        demulti_umi_dict = {}

        for idx_i in idx_dict.keys():
            demulti_fw_dict[idx_i] = open(f'{output_dir}/demultiplexing/{idx_i}.fastq', 'w')
            demulti_umi_dict[idx_i] = open(f'{output_dir}/demultiplexing/umi_{idx_i}.fasta', 'w')
        
        idx_len = [i[1] - i[0] + 1 for i in idx_pos]
        umi_len = [i[1] - i[0] + 1 for i in umi_pos]

        for line in f:
            
            lines = [line]
            if line.has_tag('SA') == True:
                for i in range(len(line.get_tag('SA').split(';')) - 1):
                    lines.append(next(f))
            
            seq = line.query_sequence
            qual = line.query_qualities

            partial_lines = []

            line_id = 0
            for line in lines:
                
                cnt_dict['all'] += 1

                if line.is_secondary:
                    cnt_dict['secondary'] += 1
                    continue

                if line.is_unmapped:
                    cnt_dict['unmap'] += 1
                    continue
        
                pos_refst = line.reference_start
                pos_refed = line.reference_end
            
                if pos_refst > front_indi_pos and pos_refed < back_indi_pos:
                    cnt_dict['filted_short'] += 1
                    continue

                elif pos_refst > front_indi_pos or pos_refed < back_indi_pos:
                    partial_lines.append(line)
                    continue
                
                cigar = line.cigar
                
                front_align = align_info(cigar, line.seq, pos_refst, front_align_cut + 25)
                back_align = align_info(cigar[::-1], rc(line.seq), ref_seq_len - pos_refed, ref_seq_len - back_align_cut - 1 + 25)

                idx_seq = ['', '']
                umi_seq = ['', '']


                idx_seq[0] = extract_seq(front_align, idx_pos[0][0], idx_pos[0][1] + 1)
                umi_seq[0] = extract_seq(front_align, umi_pos[0][0], umi_pos[0][1] + 1)
                
                idx_seq[1] = extract_seq(back_align, ref_seq_len - idx_pos[1][1] - 1, ref_seq_len - idx_pos[1][0])
                umi_seq[1] = extract_seq(back_align, ref_seq_len - umi_pos[1][1] - 1, ref_seq_len - umi_pos[1][0])
            
                check_len = False

                for i in range(2):
                    if len(idx_seq[i]) - idx_len[i] > umi_idx_r*2:
                        check_len = True
                        break
                    if len(umi_seq[i]) - umi_len[i] > umi_idx_r*2:
                        check_len = True
                        break
            
            
                if check_len == True:
                    cnt_dict['filted_idx_len'] += 1
                    continue
                
                cut_seq, cut_seq_qual = cut_align(line.seq, line.query_qualities, cigar)
                mean_qual = sum(cut_seq_qual)/len(cut_seq_qual)
                qual_string = "".join(chr(q + 33) for q in cut_seq_qual)

                if line.is_forward:
                    strand = '+'
                else:
                    strand = '-'
            
                if idx_dict != {'result': 'NNNNNNNNNN'}:
                    fid = demultiplex_idx(''.join(idx_seq), idx_dict)
                    if fid == False:
                        cnt_dict['filted_idx_len'] += 1
                        continue

                else:
                    fid = 'result'
            
                demulti_fw_dict[fid].write(f'@{line.query_name}_{line_id}\n{cut_seq}\n+\n{qual_string}\n')
                demulti_umi_dict[fid].write(f'>{line.query_name}_{line_id};{strand};{cut_seq};{mean_qual}\n{"".join(umi_seq)}\n')
                line_id += 1
                cnt_dict['pass'] += 1   
            
            # for not fully aligned reads
            
            if len(partial_lines) == 1:
                cnt_dict['filted_short'] += 1
                continue

            partial_line_info = []

            for line in partial_lines:

                pos_refst = line.reference_start
                pos_refed = line.reference_end
                pos_in_ref = False
                partial_seq = line.query_sequence
                align_strand_in_refseq = False
                align_strand_to_seq = False
                
                cigar = line.cigar
                align_len = cigar_len(cigar)

                if cigar[0][0] == 4:
                    pos_in_seq = cigar[0][1]
                    pos_in_seq_ed = pos_in_seq + len(partial_seq)
                    align_strand_to_seq = 1
                
                else:
                    pos_in_seq = seq.find(partial_seq)
                    align_strand_to_seq = 1
                    if pos_in_seq == -1:
                        pos_in_seq = seq.find(rc(partial_seq))
                        align_strand_to_seq = -1
        
                    pos_in_seq_ed = pos_in_seq + len(partial_seq)
    
                if pos_in_seq == -1:
                    cnt_dict['filted_short'] += 1
                    continue


                if line.is_reverse == False:
                    align_strand_in_refseq = 1
                elif line.is_reverse == True:
                    align_strand_in_refseq = -1
                else:
                    cnt_dict['filted_short'] += 1
                    continue

                if pos_refst < front_indi_pos and align_len > front_align_cut + 25:
                    pos_in_ref = 1
                    front_align = align_info(line.cigar, line.seq, pos_refst, front_align_cut + 25)
                    adaptive_idx = extract_seq(front_align, idx_pos[0][0], idx_pos[0][1] + 1)
                    adaptive_umi = extract_seq(front_align, umi_pos[0][0], umi_pos[0][1] + 1)
                elif pos_refed > back_indi_pos and align_len > ref_seq_len - back_align_cut + 25:
                    pos_in_ref = -1
                    back_align = align_info(line.cigar[::-1], rc(line.seq), ref_seq_len - pos_refed, ref_seq_len - back_align_cut - 1 + 25)
                    adaptive_idx = extract_seq(back_align, ref_seq_len - idx_pos[1][1] - 1, ref_seq_len - idx_pos[1][0])
                    adaptive_umi = extract_seq(back_align, ref_seq_len - umi_pos[1][1] - 1, ref_seq_len - umi_pos[1][0])
                
                else:
                    cnt_dict['filted_short'] += 1
                    continue
                
                partial_line_info.append([pos_in_ref, pos_in_seq, pos_in_seq_ed, adaptive_idx, adaptive_umi, align_strand_in_refseq, align_strand_to_seq])
            
            partial_line_info = sorted(partial_line_info, key=lambda x: x[1])

            for i in range(len(partial_line_info) - 1):
                info = partial_line_info[i: i+2]
                if info[0][5] != info[1][5]:
                    continue
                if info[0][5] * info[0][0] == -1:
                    continue
                if info[0][6] != info[1][6]:
                    continue
                if info[0][6] == 1:
                    cut_seq_qual = qual[info[0][1]:info[1][2]]
                    cut_seq = seq[info[0][1]:info[1][2]]
                    idx_seq = [info[0][3], info[1][3]]
                    umi_seq = [info[0][4], info[1][4]]
                if info[0][6] == -1:
                    cut_seq_qual = qual[info[0][1]:info[1][2]][::-1]
                    cut_seq = rc(seq[info[0][1]:info[1][2]])
                    idx_seq = [info[1][3], info[0][3]]
                    umi_seq = [info[1][4], info[0][4]]
                print(line.query_name)
                print(info)
                qual_string = "".join(chr(q + 33) for q in cut_seq_qual)
                mean_qual = sum(cut_seq_qual)/len(cut_seq_qual)
                if idx_dict != {'result': 'NNNNNNNNNN'}:
                    fid = demultiplex_idx(''.join(idx_seq), idx_dict)
                    if fid == False:
                        cnt_dict['filted_idx_len'] += 1
                        continue
                else:
                    fid = 'result'
                demulti_fw_dict[fid].write(f'@{line.query_name}_{line_id}_Splited\n{cut_seq}\n+\n{qual_string}\n')
                demulti_umi_dict[fid].write(f'>{line.query_name}_{line_id};+;{cut_seq};{mean_qual};Splited\n{"".join(umi_seq)}\n')
                line_id += 1
                cnt_dict['split_pass'] += 1
            
    print(cnt_dict)

def clustering_umi(input_file, umi_len, output_clust, output_consensus, threads, clust_cutoff, umi_r=1):

    print('Running 1st vsearch ...\r', end='')
    os.system(' '.join(['vsearch',
            '--clusterout_id',
            '--clusters', f'{output_clust}/1st_clusters/',
            '--consout', f'{output_clust}/1st_consensus.fasta',
            '--minseqlength', str(umi_len - umi_r),
            '--maxseqlength', str(umi_len + umi_r),
            '--qmask', 'none',
            '--threads', str(threads), 
            '--cluster_fast', input_file,
            '--clusterout_sort',
            '--gapopen', '0E/5I',
            '--gapext', '0E/2I',
            '--mismatch', '-8',
            '--match', '6',
            '--iddef', '0',
            '--minwordmatches', '0',
            '-id', '0.9']))
    print('Running 1st vsearch ... completed')

    print('Running 2nd vsearch ...\r', end='')   
    os.system(' '.join(['vsearch',
            '--clusterout_id',
            '--clusters', f'{output_clust}/2nd_clusters/',
            '--minseqlength', str(umi_len - umi_r),
            '--maxseqlength', str(umi_len + umi_r),
            '--qmask', 'none',
            '--threads', str(threads), 
            '--cluster_fast', f'{output_clust}/1st_consensus.fasta',
            '--clusterout_sort',
            '--gapopen', '0E/5I',
            '--gapext', '0E/2I',
            '--mismatch', '-8',
            '--match', '6',
            '--iddef', '0',
            '--minwordmatches', '0',
            '-id', '0.95']))
    print('Running 2nd vsearch ... completed')

    print('Reorganinzing clustering files ...\r', end='')
    cluster_cnt = 0
    filted_cnt = 0
    for fn in os.listdir(f'{output_clust}/2nd_clusters'):
        if fn.find('.') != -1:
            continue
        fw_list = []
        with open(f'{output_clust}/2nd_clusters/{fn}') as f2:
            for line in f2:
                v1_cluster = line.strip().replace('Splited;', '').split(';')[5].replace('clusterid=','')
                with open(f'{output_clust}/1st_clusters/{v1_cluster}') as f1:
                    for line in f1:
                        line = line.strip().split(';')
                        sid = line[0][1:]
                        seq = line[2]
                        fw_list.append('>' + sid)
                        fw_list.append(seq)
                        next(f1)
                next(f2)
        if len(fw_list) <= clust_cutoff*2:
            filted_cnt += 1
            continue
        cluster_cnt += 1
        with open(f'{output_clust}/medaka_input/{fn}.fasta', 'w') as fw:
            fw.write('\n'.join(fw_list))
    print('Reorganinzing clustering files ... completed')
    print(f'The number of clustering : {cluster_cnt} (filted cluster : {filted_cnt}')

    print('Generating consensus sequence ...\r', end='') 
    medaka_input_list = []
    for i in os.listdir(f'{output_clust}/medaka_input/'):
        if i[0] != '.' and i[-6:] == '.fasta':
            medaka_input_list.append(f'{output_clust}/medaka_input/{i}')
    fw = open(output_consensus + '/consensus.fasta', 'w')
    t = time.time()
    for fn_n, fn in enumerate(medaka_input_list):
        print(f'Generating consensus sequence [{fn_n} / {len(medaka_input_list)}] {round(time.time() - t,2)} s\r', end='') 
        t = time.time()
        reads = []
        with open(fn) as f:
            for line in f:
                reads.append(next(f).strip())
        min_cutoff = int(len(reads)*0.1)
        if min_cutoff < 3:
            min_cutoff = 3
        consensus, msa = poa(reads[:30], min_coverage=min_cutoff, genmsa=False)
        fw.write(f'>{fn[:fn.find(".")]}\n{consensus}\n')
    fw.close()
    print('Generating consensus sequence ... completed')


def extract_mutation(info): 
	cigar, is_reverse, ref_pos, align_end, seq, read_name, ref_seq, window_st, window_ed, large_ins_cut, large_del_cut, SA_tag = info
	if SA_tag != False:
		SA_tag = SA_tag.split(',')
		SA_read_strand = SA_tag[2]
		if [is_reverse, SA_read_strand] in [[False, '-'], [True, '+']]:
			return False
		SA_read_st = int(SA_tag[1]) - 1
		SA_read_cigar = SA_tag[3]
		SA_read_ed = SA_read_st
		p = '' 
		for i in SA_read_cigar:
			if i in '0123456789':
				p += i
			else:
				p = int(p)
				if i in 'MD':
					SA_read_ed += p
				p = ''
		if ref_pos < SA_read_st:
			LD_st = align_end
			LD_ed = SA_read_st - 1
		else:
			LD_st = SA_read_ed
			LD_ed = ref_pos
		if LD_ed < window_st or LD_st > window_ed:
			return False
		else:
			return [f'{LD_st}_{LD_ed}:Large_Del_{LD_ed - LD_st + 1}']
	seq_pos = 0
	mut_info = [read_name]
	for i in cigar:
		if i[0] == 0:
			if ref_pos <= window_st <= ref_pos + i[1] or ref_pos <= window_ed <= ref_pos + i[1]:
				for x in range(i[1]):
					if window_st <= ref_pos + x <= window_ed:
						if ref_seq[ref_pos + x] != seq[seq_pos + x]:
							mut_info.append(f'{ref_pos}:Sub_{ref_seq[ref_pos+x]}>{seq[seq_pos+x]}')							  
			ref_pos += i[1]
			seq_pos += i[1]
		elif i[0] == 1:
			if window_st <= ref_pos <= window_ed:
				if i[1] > large_ins_cut:
					mut_info.append(f'{ref_pos+1}_{ref_pos+2}:Large_Ins_{i[1]}')
				else:
					mut_info.append(f'{ref_pos+1}_{ref_pos+2}:Ins_{i[1]}')
			seq_pos += i[1]
		elif i[0] == 2:
			if ref_pos + i[1] < window_st or ref_pos > window_ed:
				pass
			else:
				if i[1] > large_del_cut:
					mut_info.append(f'{ref_pos+1}_{ref_pos+i[1]}:Large_Del_{i[1]}')
				else:
					mut_info.append(f'{ref_pos+1}_{ref_pos+i[1]}:Del_{i[1]}')
			ref_pos += i[1]
		elif i[0] == 4:
			seq_pos += i[1]
	return mut_info


def mutation_analysis(ref_seq, target_seq, cleavage_pos, window_r, input_file, output_dir, threads=1, large_ins_cut=20, large_del_cut=100):

    def rc(s): return s.translate(s.maketrans('ATGC','TACG'))[::-1]

    # cleavage position is defined to next position to cleavage site
    if ref_seq.find(target_seq) != -1:
        st_pos = ref_seq.find(target_seq)
        cv_pos = st_pos + cleavage_pos
    elif ref_seq.find(rc(target_seq)) != -1:
        st_pos = ref_seq.find(rc(target_seq))
        cv_pos = st_pos + (len(target_seq) - cleavage_pos)
    else:
        print('ERROR: Can not find target seuqence in reference sequence')
        sys.exit()

    file_name = input_file[:input_file.find('.')]
    cnt_dict = {'all': 0, 'ins': 0, 'del': 0, 'large_ins': 0, 'large_del': 0, 'wt': 0, 'sub':0, 'complex':0}
    window_st = cv_pos - window_r
    window_ed = cv_pos + window_r
    pool = Pool(1)
    pool_input = []
    with pysam.AlignmentFile(input_file) as f, open(f'{output_dir}/read_classification.txt', 'w') as fw:
        fw.write('Read_id\tClassification\tMutation_info\n')
        for line_n, line in enumerate(f):
            if line.is_unmapped or line.is_secondary or line.is_supplementary:
                continue
            if line_n % 1000 == 0:
                pool_output = pool.map(extract_mutation, pool_input)
                for i in pool_output:
                    if i == False:
                        continue
                    fw.write(i[0])
                    cnt_dict['all'] += 1
                    if len(i) == 1:
                        cnt_dict['wt'] += 1
                        fw.write(f'\tWT\tNone\n')
                    elif len(i) == 2:
                        mut = i[1][i[1].find(':')+1:i[1].rfind('_')]
                        cnt_dict[mut.lower()] += 1
                        fw.write(f'\t{mut}\t{i[1]}\n')
                    elif len(i) > 2:
                        cnt_dict['complex'] += 1
                        fw.write(f'\tComplex\t{",".join(i[1:])}\n')
                pool_input = []
            SA_tag = False
            if line.has_tag('SA'):
                SA_tag = line.get_tag('SA').split(';')
                if len(SA_tag) != 2:
                    SA_tag = False
                elif len(SA_tag) > 2:
                    continue
                else:
                    SA_tag = SA_tag[0]
            pool_input.append([line.cigar, line.is_reverse, line.reference_start, line.reference_end, line.query_sequence, line.query_name, ref_seq, window_st, window_ed, large_ins_cut, large_del_cut, SA_tag])
        pool_output = pool.map(extract_mutation, pool_input)
        for i in pool_output:
            if i == False:
                continue
            fw.write(i[0])
            cnt_dict['all'] += 1
            if len(i) == 1:
                cnt_dict['wt'] += 1
                fw.write(f'\tWT\tNone\n')
            elif len(i) == 2:
                mut = i[1][i[1].find(':')+1:i[1].rfind('_')]
                cnt_dict[mut.lower()] += 1
                fw.write(f'\t{mut}\t{i[1]}\n')
            elif len(i) > 2:
                cnt_dict['complex'] += 1
                fw.write(f'\tComplex\t{",".join(i[1:])}\n')
    
    print(f"All reads : {cnt_dict['all']}")
    for i in ['wt', 'sub', 'ins', 'large_ins', 'del', 'large_del', 'complex']:
        print(f"{i} : {cnt_dict[i]}  {round(cnt_dict[i]*100/cnt_dict['all'], 2)} %")

    return cv_pos


