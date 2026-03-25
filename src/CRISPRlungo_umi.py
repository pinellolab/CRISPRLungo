
import os, sys, time, pysam, editdistance
from subprocess import Popen, PIPE
from spoa import poa
from multiprocessing import Pool
import csv

def rc(s): return s.translate(s.maketrans('ATGC','TACG'))[::-1]

def check_in_quality(read):
	
	if read.mapping_quality > 30:
		return True
	if read.has_tag('SA'):
		if read.get_tag('NM') < read.query_alignment_length * 0.05:
			return True
	return False

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

	cnt_dict = {'all': 0, 'unmapped': 0, 'filted_align_short': 0, 'filted_umi_len': 0, 'used': 0, 'used_with_supple': 0} 

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
			fid = idx_i
		
		idx_len = [i[1] - i[0] + 1 for i in idx_pos]
		umi_len = [i[1] - i[0] + 1 for i in umi_pos]

		for line in f:

			cnt_dict['all'] += 1

			if line.is_unmapped:	
				cnt_dict['unmapped'] += 1
				continue

			if line.is_secondary:
				continue
			
			lines = [line]

			if line.has_tag('SA') == True:
				SA_n = len(line.get_tag('SA').split(';'))
				while len(lines) < SA_n:
					next_read = next(f)
					if next_read.is_supplementary:
						lines.append(next_read)

			#for i in range(len(line.get_tag('SA').split(';')) - 1):
			#		lines.append(next(f))

			seq = lines[0].query_sequence
			qual = lines[0].query_qualities

			partial_lines = []

			line_id = 0
			counted_check = False

			for line in lines:

				pos_refst = line.reference_start
				pos_refed = line.reference_end

				if pos_refst > front_indi_pos and pos_refed < back_indi_pos:
					#cnt_dict['filted_align_short'] += 1
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
					cnt_dict['filted_umi_len'] += 1
					continue
				
				cut_seq, cut_seq_qual = cut_align(line.seq, line.query_qualities, cigar)
				mean_qual = sum(cut_seq_qual)/len(cut_seq_qual)
				qual_string = "".join(chr(q + 33) for q in cut_seq_qual)

				if line.is_forward:
					strand = '+'
				else:
					strand = '-'

				demulti_fw_dict[fid].write(f'@{line.query_name}_{line_id}\n{cut_seq}\n+\n{qual_string}\n')
				demulti_umi_dict[fid].write(f'>{line.query_name}_{line_id};{strand};{cut_seq};{mean_qual}\n{"".join(umi_seq)}\n')
				line_id += 1
				cnt_dict['used'] += 1	
				counted_check = True
			
			# for not fully aligned reads

			if len(partial_lines) == 1 and counted_check == False:
				cnt_dict['filted_align_short'] += 1
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

				if line.is_supplementary == False:
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
					#cnt_dict['filted_align_short'] += 1
					continue


				if line.is_reverse == False:
					align_strand_in_refseq = 1
				elif line.is_reverse == True:
					align_strand_in_refseq = -1
				else:
					#cnt_dict['filted_align_short'] += 1
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
					#cnt_dict['filted_align_short'] += 1
					continue
				
				partial_line_info.append([pos_in_ref, pos_in_seq, pos_in_seq_ed, adaptive_idx, adaptive_umi, align_strand_in_refseq, align_strand_to_seq])
			
			"""if 'UMI1267' in lines[0].query_name:
				print('---')
				for i in partial_line_info:
					print(i)
				input()"""

			partial_line_info = sorted(partial_line_info, key=lambda x: x[1])

			counted_check_2 = False

			for i in range(len(partial_line_info) - 1):
				info = partial_line_info[i: i+2]
				if info[0][5] != info[1][5]:
					continue
				#if info[0][5] * info[0][0] == -1:
				#	continue
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

				qual_string = "".join(chr(q + 33) for q in cut_seq_qual)
				mean_qual = sum(cut_seq_qual)/len(cut_seq_qual)

				demulti_fw_dict[fid].write(f'@{line.query_name}_{line_id}_Splited\n{cut_seq}\n+\n{qual_string}\n')
				demulti_umi_dict[fid].write(f'>{line.query_name}_{line_id};+;{cut_seq};{mean_qual};Splited\n{"".join(umi_seq)}\n')
				line_id += 1
				cnt_dict['used_with_supple'] += 1
				counted_check = True
			if counted_check_2 == False and counted_check == False:
				cnt_dict['filted_align_short'] += 1

			"""if '0b4fc67e61bc' in lines[0].query_name:
				print(counted_check)
				print(counted_check_2)
				input()"""


	with open(output_dir + '/results/preprocess_count.txt', 'w') as fw:
		for i in cnt_dict.keys():
			fw.write(i + '\t')
		fw.write('\n')
		for i in cnt_dict.values():
			fw.write(f'{i}\t')
			
	print(cnt_dict)


def run_spoa(info):
	consensus, msa =  poa(info[0], min_coverage=info[1], genmsa=info[2])
	return consensus, info[3]

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

	with open(f'{output_clust}/1st_consensus.fasta') as f, open(f'{output_clust}/confirmed_1st_consensus.fasta', 'w') as fw:
		for sid in f:
			umi = next(f).strip()
			if abs(len(umi) - umi_len) > umi_r:
				for i in sid.split(';'):
					if i.split('=')[0] == 'clusterid':
						cluster_file_n = i.split('=')[1].strip()
						f_clust = open(f'{output_clust}/1st_clusters/{cluster_file_n}').readlines()
						cluster_d = {}
						for x in range(int(len(f_clust)/2)):
							u = f_clust[2*x+1].strip()
							if u not in cluster_d:
								cluster_d[u] = [1, abs(len(u) - umi_len)]
							else:
								cluster_d[u][0] += 1
						umi = min(cluster_d, key=lambda k: (cluster_d[k][1], -cluster_d[k][0]))
			fw.write(f'>{sid.strip()}\n{umi}\n')
					

	print('Running 2nd vsearch ...\r', end='')	 
	os.system(' '.join(['vsearch',
			'--clusterout_id',
			'--clusters', f'{output_clust}/2nd_clusters/',
			'--minseqlength', str(umi_len - umi_r),
			'--maxseqlength', str(umi_len + umi_r),
			'--qmask', 'none',
			'--threads', str(threads), 
			'--cluster_fast', f'{output_clust}/confirmed_1st_consensus.fasta',
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
	with Pool(threads) as pool:
		pool_input = []
		print(f'Generating consensus sequence [0 / {len(medaka_input_list)}] {round(time.time() - t,2)} s\r', end='') 
		for fn_n, fn in enumerate(medaka_input_list):

			t = time.time()
			reads = []
			with open(fn) as f:
				for line in f:
					reads.append(next(f).strip())
			min_cutoff = int(len(reads)*0.1)
			if min_cutoff < 3:
				min_cutoff = 3
			pool_input.append([reads[:30], min_cutoff, False, fn])
			if fn_n % 100 == 0:
				for spoa_res in pool.map(run_spoa, pool_input):
					fn = spoa_res[1]
					fw.write(f'>{fn[fn.rfind("/")+1:fn.rfind(".")]}\n{spoa_res[0]}\n')
					print(f'Generating consensus sequence [{fn_n} / {len(medaka_input_list)}]						\r', end='') 
				pool_input = []
		for spoa_res in pool.map(run_spoa, pool_input):
			fn = spoa_res[1]
			fw.write(f'>{fn[fn.rfind("/")+1:fn.rfind(".")]}\n{spoa_res[0]}\n')
			
	fw.close()
	print('Generating consensus sequence ... completed')
















