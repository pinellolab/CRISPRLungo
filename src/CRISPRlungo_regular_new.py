#!/usr/bin/env python
# coding: utf-8

import pysam
from scipy.stats import fisher_exact
from scipy.stats import chi2_contingency
import math
import csv
from Bio import SeqIO
import sys
import time
from multiprocessing import Pool, Pipe

import bisect
import numpy as np
import pandas as pd
from collections import Counter

def cal_test (input_val):
	table, key, umi_clustered = input_val
	control_count = table[0][0]
	treated_count = table[1][0]
	# Perform Fisher's Exact Test if control count is less than  or equal to 5 or chi2 if control_count >5
	if control_count == 0 and treated_count != 0:
		if key[2] < 5 and table[1][0]/sum(table[1])  < 0.001 and umi_clustered == False:
			p_value = 2
		else:
			p_value = 0
	elif control_count <= 5:
		_, p_value = fisher_exact(table, alternative='two-sided')
	if control_count > 5:
		chi2, p_value, dof, expected = chi2_contingency(table)
	return p_value, key, table[1][0] / sum(table[1])

def check_in_quality(read):
	
	if read.mapping_quality > 30:
		return True
	if read.has_tag('SA'):
		if read.get_tag('NM') < read.query_alignment_length * 0.05:
			return True
	return False

def check_in_window(mutation, window_filter, cv_pos, cv_pos_2, window, check_window_between_targets):

	if window_filter == False:
		return True
	
	window_st = cv_pos - window
	window_ed = cv_pos + window
	
	mut_st = mutation[1]
	if mutation[0] == 'insertion':
		mut_ed = mutation[4]
	elif mutation[0] == 'substitution':
		mut_ed = mut_st + mutation[2] - 1
	else:
		mut_ed = mut_st + mutation[2]

	check_position_in_window = False
	
	if mut_ed < window_st or mut_st > window_ed:
		pass
	else:
		check_position_in_window = True
	
	if cv_pos_2 != False:

		if check_window_between_targets:
			window_st2 = cv_pos
		else:
			window_st2 = cv_pos_2 - window

		window_ed2 = cv_pos_2 + window
		if mut_ed < window_st2 or mut_st > window_ed2:
			pass
		else:
			check_position_in_window = True

	return check_position_in_window

def reverse_complementary(s):
		return s.translate(s.maketrans('ATGCatgc','TACGtacg'))[::-1]

def cigar_len(cigar):
	length = 0
	for i in cigar:
		if i[0] in [0, 2]:
			length += i[1]
	return length

def analyze_SA_reads(partial_reads, query_seq, reference_sequence, reference_length, query_name, range_align_end=100, v = False):

	# Initialize an empty list to store information about partial reads
	
	partial_read_info = []
	
	for read in partial_reads:

		pos_refst = read.reference_start
		pos_refed = read.reference_end
		
		pos_in_ref = False
		partial_seq = read.query_sequence
		query_seq_len = len(query_seq)
		align_strand_in_refseq = 2 # 2 = not aligned, -1 = reverse, 1 = forward, 3 = full intact read
		align_strand_to_seq = 2 # 2 = not aligned, -1 = reverse, 1 = forward

		cigar = read.cigar
		align_len = cigar_len(cigar)

		pos_in_seq = -1
		
		align_strand_info = read.query_name.split('alignStrandInfo_')[1].split('_')
		align_strand_to_seq = 1
		for i in align_strand_info[1:]:
			align_strand_to_seq *= int(i)
	
		if query_seq.find(partial_seq) != -1:
			if cigar[0][0] in [4, 5]:
				pos_in_seq = cigar[0][1]
			else:
				pos_in_seq = 0
			if cigar[-1][0] in [4, 5]:
				pos_in_seq_ed = query_seq_len - cigar[-1][1] - 1
			else:
				pos_in_seq_ed = query_seq_len - 1
			#align_strand_to_seq = 1
		elif query_seq.find(reverse_complementary(partial_seq)) != -1:
			if cigar[0][0] in [4, 5]:
				pos_in_seq_ed = query_seq_len - cigar[0][1] - 1
			else:
				pos_in_seq_ed = query_seq_len - 1
			if cigar[-1][0] in [4, 5]:
				pos_in_seq = cigar[-1][1]
			else:
				pos_in_seq = 0
			#align_strand_to_seq = -1

		# If sequence position is not found, continue to the next read


		if pos_in_seq == -1:
			continue

		if pos_refst < range_align_end or pos_refed > reference_length -range_align_end:
			pass
		else:
			continue
		
		# Check if the read spans almost the entire reference sequence

		if pos_refst < range_align_end and pos_refed > reference_length - range_align_end:
			align_strand_in_refseq = 3	# Mark as an intact read
		elif not read.is_reverse:  # If the read is aligned in the forward direction
			align_strand_in_refseq = 1
		elif read.is_reverse:  # If the read is aligned in the reverse direction
			align_strand_in_refseq = -1
		else:
			continue

		partial_read_info.append([pos_refst, pos_refed, pos_in_seq, pos_in_seq_ed, align_strand_in_refseq, align_strand_to_seq, cigar, partial_seq])
	
	partial_read_info = sorted(partial_read_info, key=lambda x: x[2])
	partial_read_info.append((2, 2, 2, 2, 2, 2, '')) # To check last partial read is fully aligned

	mutations_in_reads = []


	# Loop through the list of partial reads' information

	for i in range(len(partial_read_info) - 1):
		mutations_in_read = []
		info = partial_read_info[i: i+2]  # Take two consecutive reads for comparison

		if info[1][4] == 3:
			continue

		# Check if the read is an intact read (spans most of the reference)
		
		if info[0][4] == 3:
			ref_pos = info[0][0]  # Starting position in the reference
			seq_pos = 0
			using_query_seq = info[0][7]  # Use the query sequence directly

			for operation, length in info[0][6]:
				if operation == 0:	
					sub_tmp_list = []
					for x in range(length):
						if using_query_seq[seq_pos + x] != reference_sequence[ref_pos + x]:
							mutations_in_read.append(('substitution', ref_pos+x, 1, reference_sequence[ref_pos + x], using_query_seq[seq_pos + x]))
							"""
							if len(sub_tmp_list) > 0 and sub_tmp_list[-1][1] + sub_tmp_list[-1][2] == ref_pos + x:
								sub_tmp_list[-1][2] += 1
								sub_tmp_list[-1][3] += reference_sequence[ref_pos + x]
								sub_tmp_list[-1][4] += query_seq[seq_pos + x]
							else:
								sub_tmp_list.append(['substitution', ref_pos + x, 1, reference_sequence[ref_pos + x], query_seq[seq_pos + x]])
					for sub in sub_tmp_list:
						mutations_in_read.append(tuple(sub))
							"""
				if operation == 1:	
					mutations_in_read.append(('insertion', ref_pos, length, using_query_seq[seq_pos: seq_pos + length], ref_pos + 1,seq_pos, seq_pos + length))
				elif operation == 2:  
					mutations_in_read.append(('deletion', ref_pos, length))
				
				# Adjust the reference and sequence positions after operations
				if operation in [0, 2, 3]: 
					ref_pos += length
				if operation in [0, 1,4 ]:
					seq_pos += length
			
			# Append the detected mutations in the current read
			mutations_in_reads.append(mutations_in_read)
			continue
		
		# Process reads aligned to the forward strand
		#if info[0][5] == 1 and info[1][5] == 1:

		if info[0][5] == -1:
			used_query_seq = reverse_complementary(query_seq)
		else:
			used_query_seq = query_seq


		if info[0][5] == info[1][5] and info[0][4] == info[1][4] and info[0][1] < info[1][0]:
			
			if info[0][1] < info[1][0] - 100:  # Check for a deletion between the two reads
				mutations_in_read.append(('deletion', info[0][1], info[1][0] - info[0][1]))
			if abs(info[0][3] - info[1][2]) - 1 > 20:  # Check for an insertion between the reads
				ins_end_pos = info[1][0]
				if ins_end_pos <= info[0][1]:
					ins_end_pos = info[0][1] + 1
				if query_seq.find(info[0][7]) != -1:
					ins_seq = query_seq[info[0][3] + 1: info[0][3] + 1 + abs(info[0][3] - info[1][2]) - 1]
				else:
					ins_seq = reverse_complementary(query_seq[info[0][3] + 1: info[0][3] + 1 + abs(info[0][3] - info[1][2]) - 1])
				mutations_in_read.append(('insertion', info[0][1], abs(info[0][3] - info[1][2]) - 1,  ins_seq, ins_end_pos, info[0][3], info[0][3] + abs(info[0][3] - info[1][2])))

			# Iterate through both reads and check for mutations
			for sub_info in info:
				ref_pos = sub_info[0]  
				seq_pos = 0
				used_query_seq = sub_info[7]
				for operation, length in sub_info[6]:
					if operation == 0:	
						sub_tmp_list = []
						for x in range(length):
							"""
							if query_seq[seq_pos + x] != reference_sequence[ref_pos + x]:
								if len(sub_tmp_list) > 0 and sub_tmp_list[-1][1] + sub_tmp_list[-1][2] == ref_pos + x:
									sub_tmp_list[-1][2] += 1
									sub_tmp_list[-1][3] += reference_sequence[ref_pos + x]
									sub_tmp_list[-1][4] += query_seq[seq_pos + x]
								else:
									sub_tmp_list.append(['substitution', ref_pos + x, 1, reference_sequence[ref_pos + x], query_seq[seq_pos + x]])
						for sub in sub_tmp_list:
							mutations_in_read.append(tuple(sub))
							"""
							if used_query_seq[seq_pos + x] != reference_sequence[ref_pos + x]:
								mutations_in_read.append(('substitution', ref_pos+x, 1, reference_sequence[ref_pos + x], used_query_seq[seq_pos + x]))	
					elif operation == 1:  
						mutations_in_read.append(('insertion', ref_pos, length, used_query_seq[seq_pos: seq_pos + length], ref_pos + 1, seq_pos, seq_pos + length))
					elif operation == 2: 
						mutations_in_read.append(('deletion', ref_pos, length))
					
					if operation in [0, 2, 3]:
						ref_pos += length
					if operation in [0, 1, 4]:
						seq_pos += length
			
			mutations_in_reads.append(mutations_in_read)

		# Process reads aligned to the reverse strand
		#elif info[0][5] == -1 and info[1][5] == -1:
		elif info[0][5] ==	info[1][5] and info[0][4] == info[1][4]  and info[1][1] < info[0][0]:
			if info[1][1] < info[0][0] - 100:  # Check for a deletion between the two reads
				mutations_in_read.append(('deletion', info[1][1], info[0][0] - info[1][1] - 1))
			if abs(info[1][2] - info[0][3]) - 1 > 20:  # Check for an insertion between the reads
				ins_end_pos = info[0][0]
				if ins_end_pos <= info[1][1]:
					ins_end_pos = info[1][1] + 1
				if query_seq.find(info[0][7]) != -1:
					ins_seq = query_seq[info[0][3] + 1: info[0][3] + 1 + abs(info[1][2] - info[0][3]) - 1]
				else:
					ins_seq = reverse_complementary(query_seq[info[0][3] + 1: info[0][3] + 1 + abs(info[1][2] - info[0][3]) - 1])
				mutations_in_read.append(('insertion', info[1][1], abs(info[1][2] - info[0][3]) - 1, ins_seq, ins_end_pos, info[0][3], info[0][3] + (info[1][2] - info[0][3] - 1)))

			# Reverse complement query sequence for reverse strand processing
			query_seq_len = len(query_seq)

			# Process mutations for both reads
			for sub_info in info:
				ref_pos = sub_info[0]  
				seq_pos = 0
				used_query_seq = sub_info[7]
				
				for operation, length in sub_info[6]:
					if operation == 0:	
						sub_tmp_list = []
						for x in range(length):
							if used_query_seq[seq_pos + x] != reference_sequence[ref_pos + x]:
								mutations_in_read.append(('substitution', ref_pos+x, 1, reference_sequence[ref_pos + x], used_query_seq[seq_pos + x]))
								"""
								if len(sub_tmp_list) > 0 and sub_tmp_list[-1][1] + sub_tmp_list[-1][2] == ref_pos + x:
									sub_tmp_list[-1][2] += 1
									sub_tmp_list[-1][3] += reference_sequence[ref_pos + x]
									sub_tmp_list[-1][4] += query_seq[query_pos + x]
								else:
									sub_tmp_list.append(['substitution', ref_pos + x, 1, reference_sequence[ref_pos + x], query_seq[query_pos + x]])
						for sub in sub_tmp_list:
							mutations_in_read.append(tuple(sub))
								"""
					elif operation == 1: 
						mutations_in_read.append(('insertion', ref_pos, length, used_query_seq[seq_pos: seq_pos + x], ref_pos + 1, seq_pos, seq_pos+x))
					elif operation == 2: 
						mutations_in_read.append(('deletion', ref_pos, length))
					
					if operation in [0, 2, 3]:
						ref_pos += length
					if operation in [0, 1, 4]:
						seq_pos += length
			
			mutations_in_reads.append(mutations_in_read)

	return mutations_in_reads


def analysis_function(control, edited, refernce, output_dir, cv_pos, cv_pos_2, window, check_window_between_targets, induced_mutations, range_align_end=100, use_all_mutations=False, length_min = 10, allowance_value = 0.05, pooling = True, largeins_cutlen=50, largedel_cutlen=50, Filter1 = True, window_filter = True, threads=1, umi_clustered = False):

	# Path to SAM files
	samfile_path1 = edited
	samfile_path2 = control

	fasta_file = refernce

	def find_sub(read,fasta_file):
		reference_sequence = ""
		for record in SeqIO.parse(fasta_file, "fasta"):
			reference_sequence = str(record.seq).upper()
			break

		# Open the SAM file
		substitution_list = []
	   
		read_seq = read.query_sequence
		
		# Iterate over the aligned portion of the read
		for read_idx, ref_idx in enumerate(read.get_reference_positions(full_length=True)):
			if ref_idx is None:
				continue
			read_base = read_seq[(read_idx)]
			if read_base != reference_sequence[ref_idx]:
				substitution_list.append(("substitution",ref_idx,1,reference_sequence[ref_idx],read_base))

		# Output the counts for each position
		return(substitution_list)

	def combine_deletions(tuples):
		# Function to check if two deletions should be combined
		if len(tuples)<=1:return(tuples)

		def should_combine(pos1, len1, pos2,len2):
			return (abs(pos2 - (pos1 + len1)) <= 0.05 * (len1+len2)/2 and len1>=5 and len2>=5)

		how_many_pop =0
		combined_tuples = []
		for i in range(0,len(tuples)):
			tuple1 = tuples[i-how_many_pop]
			if list(tuple1)[2]<=5: 
				tuples.pop(i-how_many_pop)
				how_many_pop+=1
				combined_tuples.append(tuple1)
		i = 0
		while i < len(tuples)-1:
			if tuples[i][0] == "insertion" or tuples[i][0] == "deletion" and tuples[i+1][0] == tuples[i][0]:
				pos1, len1 = tuples[i][1], tuples[i][2]
				j = i + 1
				while j < (len(tuples)-1) and tuples[j][0] == tuples[i][0] and should_combine(pos1, len1, tuples[j][1], tuples[j][2]):
					len1 += tuples[j][2] + (tuples[j][1] - (pos1 + len1))
					j += 1
				combined_tuples.append((tuples[j][0], pos1, len1))
				i = j
			else:
				combined_tuples.append(tuples[i])
				i += 1
		return combined_tuples
	
	def quant_unique_indels(samfile_path):
		# Open the SAM file
		fw = open('short.txt', 'w')
		samfile = pysam.AlignmentFile(samfile_path, 'r')
		reference_sequence = ""
		for record in SeqIO.parse(fasta_file, "fasta"):
			reference_sequence = str(record.seq).upper()
		reference_length = len(reference_sequence)

		# Initialize a dictionary to hold mutation information
		# Key: (type_of_mutation, position_of_mutation, length_of_mutation)
		# Value: count of occurrences
		mutations = {}

		#initilize dict where read number is key and value is list of mutations
		reads_cnt = {'all_reads': 0, 'unmapped': 0, 'low_quality': 0, 'short':0, 'used': 0}
		counter = 0
		# Iterate over reads in the SAM file
		for read in samfile.fetch():
			mutations_in_read = []
			ref_pos = read.reference_start
			# Skip unmapped reads

			if read.is_secondary or read.is_supplementary:
				continue

			reads_cnt['all_reads'] += 1

			if read.is_unmapped:
				reads_cnt['unmapped'] += 1
				continue

			if not check_in_quality(read):
				reads_cnt['low_quality'] += 1
				continue
			
			# Current position in the reference sequence
			"""if read.reference_start > 100 and read.reference_end < reference_length - 100:
				reads_cnt['short'] += 1
				continue"""

			
			if read.has_tag('SA'):
				
				SA_reads = [read]
				SA_n = len(read.get_tag('SA').split(';'))
				ori_query_seq = read.query_sequence.upper()
				
				while len(SA_reads) < SA_n:
					next_read = next(samfile)
					if next_read.is_supplementary:
						SA_reads.append(next_read)
				
				SA_reads_muts =  analyze_SA_reads(SA_reads, ori_query_seq, reference_sequence, reference_length, read.query_name, range_align_end)

				if SA_reads_muts == []:
					fw.write(read.query_name + '\n')
					reads_cnt['short'] += 1
					continue

				for SA_mutations in SA_reads_muts:
					reads_cnt['used'] += 1
					for mut in SA_mutations:
						if check_in_window(mut, window_filter, cv_pos, cv_pos_2, window, check_window_between_targets) == True:
							if mut[0] == 'insertion':
								mut = mut[:4]
							mutations[mut] = mutations.get(mut, 0) + 1
				
			else:
				if read.reference_start>=range_align_end or read.reference_end<(len(reference_sequence)-range_align_end):
					fw.write(read.query_name + '\n')
					reads_cnt['short'] += 1
					continue
				if not check_in_quality(read):
					reads_cnt['low_quality'] += 1
					continue
				query_seq = read.query_sequence
				query_pos = 0
				reads_cnt['used'] += 1
				# Iterate over CIGAR operations
				for operation, length in read.cigar:
					# Check for insertion (I) or deletion (D)
					if operation == 0:	
						sub_tmp_list = []
						for x in range(length):
							if query_seq[query_pos + x] != reference_sequence[ref_pos + x]:
								mutations_in_read.append(('substitution', ref_pos+x, 1, reference_sequence[ref_pos + x], query_seq[query_pos + x]))
								"""
								if len(sub_tmp_list) > 0 and sub_tmp_list[-1][1] + sub_tmp_list[-1][2] == ref_pos + x:
									sub_tmp_list[-1][2] += 1
									sub_tmp_list[-1][3] += reference_sequence[ref_pos + x]
									sub_tmp_list[-1][4] += query_seq[query_pos + x]
								else:
									sub_tmp_list.append(['substitution', ref_pos + x, 1, reference_sequence[ref_pos + x], query_seq[query_pos + x]])
						for sub in sub_tmp_list:
							mutations_in_read.append(tuple(sub))
								"""
					elif operation == 1:	# Insertion
						#key = ('insertion', ref_pos, length, query_seq[query_pos: query_pos + length], ref_pos + 1)
						key = ('insertion', ref_pos, length, query_seq[query_pos: query_pos+length], ref_pos + 1, query_pos, query_pos+length)
						mutations_in_read.append(key)
					elif operation == 2 :  # Deletion
						key = ('deletion', ref_pos, length)
						mutations_in_read.append(key)
					# Update ref_pos based on operation
					if operation in [0, 2, 3]:	# Match/Mismatch, Deletion, N (Skipped region)
						ref_pos += length
					if operation in [0, 1, 4]:
						query_pos += length

					counter +=1

				#combined_mutations_in_read = combine_deletions(mutations_in_read)
				#if len(combined_mutations_in_read) != len(mutations_in_read):
				   # combined_mutations_in_read = combine_deletions(combined_mutations_in_read)
			
				#substitution counter

				#subs_in_read = find_sub(read,fasta_file)
				#for sub in subs_in_read: mutations_in_read.append(sub)
				for mut in mutations_in_read:
					if check_in_window(mut, window_filter, cv_pos, cv_pos_2, window, check_window_between_targets) == True:
						if mut[0] == 'insertion':
							mut = mut[:4]
						mutations[mut] = mutations.get(mut, 0) + 1
				
				

			
		# Close the SAM file
		samfile.close()
		print('Done')

		return mutations, reads_cnt


	def custom_round(num,margin,length):
		if math.floor(length*margin)  == 0: return(num)
		closest_value = round(num/math.floor(length*margin))*math.floor(length*margin)
		return int(closest_value)

	def pool(mut_dict, allowance = allowance_value):
		outputdict = {}
		for key in mut_dict.keys():

			if math.floor(list(key)[2]*allowance)  == 0: 
				outputdict[key] = mut_dict[key]
			else: 
				newlength = custom_round(list(key)[2],allowance,list(key)[2])
				newpos = custom_round(list(key)[1],allowance,list(key)[2])
				type = list(key)[0]
				mutation = (type,newpos,newlength)
				outputdict[mutation] = outputdict.get(mutation, 0) + mut_dict[key]
		return(outputdict)

	def significant_mutations(control_dict, edited_dict, control_reads, edited_reads, use_all_mutations = None, length_min=10):
		"""
		Performs a Fisher's Exact Test on each mutation to compare its frequency between
		a control and an edited dataset. Returns the keys of the significant mutations.

		Parameters:
		- control_dict: Dictionary of mutations for the control dataset.
		- edited_dict: Dictionary of mutations for the edited dataset.
		- control_reads: Total number of reads in the control dataset.
		- edited_reads: Total number of reads in the edited dataset.

		Returns:
		- A list of keys for the mutations that are significantly different.
		"""

		def _lengths_from_control(control_dict):
			del_lens, ins_lens, inv_lens = [], [], []
			for mut, mut_n in control_dict.items():
				L = int(mut[2])
				t = mut[0]
				if t in ('deletion'):
					del_lens.append(L)
				elif t == 'insertion':
					ins_lens.append(L)
				#elif t == 'inversion':
				#	inv_lens.append(L)
			return sorted(del_lens), sorted(ins_lens)

		# Ensure counters and sizes exist
		ctrl_counts = control_dict
		edit_counts = edited_dict
		total_ctrl = control_reads
		total_edit = edited_reads

		# Recompute distributions to guarantee the deletion list
		ctrl_del_sorted, ctrl_ins_sorted = _lengths_from_control(control_dict)

		def probability_based_on_length(mutation, ins_sorted, del_sorted):#, inv_sorted):
			"""
			Tail-percentile by length within control of the same class.
			- insertions vs ins_sorted
			- (small) deletions vs del_sorted
			- inversions vs inv_sorted
			Substitutions get 1.0.
			Uses (length-1) with bisect_right to count strictly-longer events.
			"""
			mutation = list(mutation)
			length = mutation[2]
			mtype = mutation[0]
			
			if mtype.startswith('substitution'):
				return 1.0
			if mtype == 'insertion':
				arr = ins_sorted
			elif mtype in ('deletion', 'long_deletion'):
				arr = del_sorted
			#elif mtype == 'inversion':
			#	arr = inv_sorted
			else:
				return 0.0

			if not arr:
				return 0.0

			idx = bisect.bisect_right(arr, length - 1)  # strictly longer
			return 1 - idx / len(arr)

		def probability_based_on_count(mutation, ctrl_counts, edit_counts, total_ctrl, total_edit):
			"""
			2x2 enrichment test. Only proceeds when the edited fraction > control fraction.
			Fisher if any expected <= 5; otherwise Chi-square.
			"""
			a = ctrl_counts.get(mutation, 0)
			c = edit_counts.get(mutation, 0)
			b = total_ctrl - a
			d = total_edit - c

			# Only continue if mutation is enriched in edited
			if a/total_ctrl >= c/total_edit:
				return 1.0

			row1, row2 = a + b, c + d
			col1, col2 = a + c, b + d
			N = row1 + row2
			exp = [row1*col1/N, row1*col2/N, row2*col1/N, row2*col2/N]

			table = [[a, b], [c, d]]
			if any(e <= 5 for e in exp):
				_, p = fisher_exact(table)
			else:
				_, p, _, _ = chi2_contingency(table)
			return p

		# 3) Score each unique edited mutation with lightweight % progress
		records = []
		total = len(edit_counts)
		interval = max(1, total // 100)
		t0 = time.time()

		for idx, (mut, cnt_edit) in enumerate(edit_counts.items(), start=1):
			p_len = probability_based_on_length(mut, ctrl_ins_sorted, ctrl_del_sorted)
			cnt_ctrl = ctrl_counts.get(mut, 0)
			p_cnt = probability_based_on_count(mut, ctrl_counts, edit_counts, total_ctrl, total_edit)

			records.append({
				"mutation": mut,
				"prob_length": p_len,
				"count_edited": cnt_edit,
				"count_control": cnt_ctrl,
				"prob_count": p_cnt
			})

			if idx % interval == 0 or idx == total:
				pct = idx / total * 100
				print(f"\rScored mutations: {idx}/{total} ({pct:.1f}%)", end='', flush=True)

		print(f"\nScored {total} mutations in {time.time() - t0:.1f}s")

		# (2) build output records and assign double_min
		out_records = []
		len_df = 0
		for rec in records:
			mut = rec["mutation"]
			cnt_ctrl = rec["count_control"]


			if cnt_ctrl == 0 and mut[2] > length_min:
				double_min = -1   
			else:
				double_min = rec["prob_count"]
				len_df += 1

			rec["double_min"] = double_min
			rec["bh_threshold"] = None
			rec["significant_bh"] = False
	
			out_records.append(rec)

		#len_df = len(out_records)


		valid = [r for r in out_records if r["double_min"] is not None]
		valid_sorted = sorted(valid, key=lambda x: x["double_min"])
		alpha = 0.05 

		n = 0
		for i in sorted(out_records, key= lambda x: x['double_min']):
			if i['double_min'] == -1:
				continue
			n += 1
			i["bh_threshold"] = (n / len_df) * alpha
			i["significant_bh"] = (i['double_min'] <= i["bh_threshold"])

		if n != 0:
			p_cutoff = 0
			for i in out_records:
				if i['double_min'] != -1 and i['significant_bh'] == True:
					if i['double_min'] > p_cutoff:
						p_cutoff = i['double_min']
						
			for i in out_records:
				if (i["double_min"] != -1 and i["double_min"] <= p_cutoff) or (i["double_min"] == -1) or (use_all_mutations == True):
					i['significant'] = True
				else:
					i['significant'] = False
		else:
			p_cutoff = None
			for i in out_records:
				if (i["double_min"] == -1) or (use_all_mutations == True):
					i['significant'] = True
				else:
					i['significant'] = False

		print(f"p_cutoff: {p_cutoff}")

		significant_keys = [i['mutation'] for i in out_records if i['significant']]

		#df.to_csv(output_dir + "/mutation_pattern_p_values.txt", sep="\t", index=False)
		#significant_keys = signif_mutations['mutation']
		#pvalues = signif_mutations['double_min']

		with open(output_dir + "/mutation_pattern_p_values.txt", "w") as fw:	
			header = ["mutation", "prob_length", "count_edited",
                      "count_control", "prob_count", "double_min",
                      "bh_threshold", "significant_bh", 'significant']
			fw.write("\t".join(header) + "\n")
			for i in out_records:
				row = [str(i.get(h, "")) for h in header]
				fw.write("\t".join(row) + "\n")


		print(f"Wrote CSV in {time.time() - t0:.1f}s")
		print(f"Total runtime: {time.time() - start_total:.1f}s" if 'start_total' in locals() else "Done.")

		return significant_keys


	# Separate insertion and deletion counts for histogram plotting
	 

	def creat_dict_analysis(samfile_path, significant_keys, induced_mutations):
		# Open the SAM file
		reference_sequence = ""
		for record in SeqIO.parse(fasta_file, "fasta"):
			reference_sequence = str(record.seq).upper()
		reference_length = len(reference_sequence)
		samfile = pysam.AlignmentFile(samfile_path, 'r')

		if Filter1:
			print("Proceeding With Statistical Tests... ")
		else:
			print("Reporting All Mutations... ")

		# Initialize a dictionary to hold mutation information
		# Key: (type_of_mutation, position_of_mutation, length_of_mutation)
		# Value: count of occurrences
		mutations = {}
		#initilize dict where read number is key and value is list of mutations
		dict_of_reads = {}
		total_reads = -1
		set_sig_mut = set(significant_keys)
		# Iterate over reads in the SAM file
		reads_that_passed = []
		count_printed = 0
		a = 0
		b = 0
		for read in samfile.fetch():
			mutations_in_read = []
			# Skip unmapped reads
			if read.is_unmapped or read.is_secondary or read.is_supplementary:
				continue

			if not check_in_quality(read):
				continue

			#create more varible way to quantify end
			v = False

			if read.has_tag('SA'):
				SA_reads = [read]
				SA_n = len(read.get_tag('SA').split(';'))
				ori_query_seq = read.query_sequence.upper()

				while len(SA_reads) < SA_n:
					next_read = next(samfile)
					if next_read.is_supplementary:
						SA_reads.append(next_read)

				SA_reads_muts = analyze_SA_reads(SA_reads, ori_query_seq, reference_sequence, reference_length, read.query_name, range_align_end, v = True)

				if len(SA_reads_muts) == 1:
					split_reads_check = True
				else:
					split_reads_check = False

				for split_n, mutations_in_read in enumerate(SA_reads_muts):

					total_reads += 1
					dict_of_reads[total_reads] = [[], []]
					if split_reads_check:
						reads_that_passed.append(read.query_name)
					else:
						reads_that_passed.append(read.query_name + '_split_' + str(split_n))
						
					if Filter1:
						for mut in mutations_in_read:
							allowance = allowance_value
							mut = list(mut)
							length1 = mut[2]
							pos = mut[1]
							mut_type = mut[0]
							if pooling:
								if (math.floor(length1*allowance)  == 0) :
									#pos = mut[1]
									#mut_type = mut[0]
									if mut[0] == 'insertion':
										mutation = tuple(mut[:4])
									else:
										mutation = tuple(mut)
								else:
									length1 = custom_round(length1,allowance,length1)
									pos = custom_round(list(mut)[1],allowance,length1)
									mut_type = mut[0]
									mutation = (mut_type,pos,length1)

							if check_in_window(mut, window_filter, cv_pos, cv_pos_2, window, check_window_between_targets) == True:
								if mutation in set_sig_mut or tuple(mut) in induced_mutations: 
									dict_of_reads[total_reads][0].append(mut)
								else:
									dict_of_reads[total_reads][1].append(mut)
					else:
						for mut in mutations_in_read:
							if check_in_window(mut, window_filter, cv_pos, cv_pos_2, window, check_window_between_targets) == True:
								dict_of_reads[total_reads][0].append(mut)

			else:

				if read.reference_start>range_align_end or read.reference_end<(len(reference_sequence)-range_align_end):
					continue
				if not check_in_quality(read):
					continue

				reads_that_passed.append(read.query_name)

				total_reads+=1
				#if total_reads%1000 == 0: print(read.reference_start)
				dict_of_reads[total_reads] = [[], []]
				# Current position in the reference sequence
				ref_pos = read.reference_start

				query_seq = read.query_sequence
				query_pos = 0
				
				# Iterate over CIGAR operations
				for operation, length in read.cigar:
					# Check for insertion (I) or deletion (D)
					if operation == 0:	
							sub_tmp_list = []
							for x in range(length):
								if query_seq[query_pos + x] != reference_sequence[ref_pos + x]:
									mutations_in_read.append(('substitution', ref_pos+x, 1, reference_sequence[ref_pos + x], query_seq[query_pos + x]))
									"""
									if len(sub_tmp_list) > 0 and sub_tmp_list[-1][1] + sub_tmp_list[-1][2] == ref_pos + x:
										sub_tmp_list[-1][2] += 1
										sub_tmp_list[-1][3] += reference_sequence[ref_pos + x]
										sub_tmp_list[-1][4] += query_seq[query_pos + x]
									else:
										sub_tmp_list.append(['substitution', ref_pos + x, 1, reference_sequence[ref_pos + x], query_seq[query_pos + x]])
							for sub in sub_tmp_list:
								mutations_in_read.append(tuple(sub))
									"""
					if operation == 1:	# Insertion
						key = ('insertion', ref_pos, length, query_seq[query_pos: query_pos + length], ref_pos + 1, query_pos, query_pos+length)
					elif operation == 2:  # Deletion
						key = ('deletion', ref_pos, length)

					if operation in [0, 2, 3]:	# Match/Mismatch, Deletion, N (Skipped region)
						ref_pos += length
					if operation in [0, 1, 4]:
						query_pos += length
	
					if operation not in [1,2]: continue
	
					mutations_in_read.append(key)
					
				#subs_in_read = find_sub(read,fasta_file)
				#for sub in subs_in_read: 
				#	mutations_in_read.append(sub)
				
				#combined_mutations_in_read = combine_deletions(mutations_in_read)
				#if len(combined_mutations_in_read) != mutations_in_read: 
					#combined_mutations_in_read = combine_deletions(combined_mutations_in_read)
				

				if Filter1:
					for mut in mutations_in_read:
						allowance = allowance_value
						mut = list(mut)
						length1 = mut[2]
						pos = mut[1]
						mut_type = mut[0]
						if pooling:
							if (math.floor(length1*allowance)  == 0 or math.floor(list(mut)[1]*allowance)==0) :
								#pos = mut[1]
								#mut_type = mut[0]
								if mut[0] == 'insertion':
									mutation = tuple(mut[:4])
								else:
									mutation = tuple(mut)
							else:
								length1 = custom_round(length1,allowance,list(mut)[2])
								pos = custom_round(mut[1],allowance,list(mut)[2])
								mut_type = mut[0]
								mutation = (mut_type,pos,length1)
						if check_in_window(mut, window_filter, cv_pos, cv_pos_2, window, check_window_between_targets) == True and (mutation in set_sig_mut or tuple(mut) in induced_mutations):
							dict_of_reads[total_reads][0].append(mut)
						else:
							dict_of_reads[total_reads][1].append(mut)

				else: 
					for mut in mutations_in_read:
						if check_in_window(mut, window_filter, cv_pos, cv_pos_2, window, check_window_between_targets) == True:
							dict_of_reads[total_reads][0].append(mut)

				# Update ref_pos based on operation

		for read_n in dict_of_reads.keys():
			tmp_mut = dict_of_reads[read_n]
			dict_of_reads[read_n] = [sorted(tmp_mut[0], key= lambda x: x[1]), sorted(tmp_mut[1], key= lambda x: x[1])]


		# Close the SAM file
		samfile.close()
		return(dict_of_reads,reads_that_passed)

	
	print("reading Control SAM files...\r", end='')
	control1, control_reads_cnt = quant_unique_indels(samfile_path2)
	print(f"Control align : {control_reads_cnt}")
	print("reading Treated SAM files... Done")
	edited1, edited_reads_cnt = quant_unique_indels(samfile_path1)
	print(f"Treated align : {edited_reads_cnt}")
	print("reading CRISPR-treated SAM file...\r", end='')
	
	
	print("reading CRISPR-treated SAM file... Done")

	if pooling: 
		print('Binning mutation information in Control file ...\r', end='')
		control1 = pool(control1)
		print('Binning mutation information in Treated file ... \r', end='')
		edited1 = pool(edited1)
		print('Binning mutation information ... Done!			  ')

	
	print(window)
	significant_keys = significant_mutations(control1, edited1, control_reads_cnt['used'], edited_reads_cnt['used'], use_all_mutations = use_all_mutations, length_min = length_min)
	print("Calculating mutation's significant... Done!			   ")


	print("Calculating mutation in each reads.. \r", end ='')
	edited_dict_reads, reads_that_passed = creat_dict_analysis(samfile_path1, significant_keys, induced_mutations)

	"""for x,y in edited_dict_reads.items():
		if reads_that_passed[x].find('807e023b')!= -1:
			print(y)
			input()"""
	
	#edited_dict_reads2, _ = creat_dict_analysis(samfile_path2,significant_keys, induced_mutations)
	edited_dict_reads2 = {}
	print("Calculating mutation in each reads... Done!		 ")

	"""sub_count_edited = 0
	sub_count_control = 0
	for value in edited_dict_reads.values():
		for mutation in value:
			if mutation[0][0:3] == "sub": sub_count_edited+=1
	#for value in edited_dict_reads2.values():
	#	for mutation in value:
	#		if mutation[0][0:3] == "sub": sub_count_control+=1"""

	return edited_dict_reads, edited_dict_reads2, reads_that_passed, control_reads_cnt, edited_reads_cnt



def classify_mut_mild(mutations, induced_mutations, largeins_cutlen, largedel_cutlen, partial_induce_cutoff=0.8): #((mut_type, position, length))
	if mutations == []:
		return 'WT', 'None', '', 'X'
	muts = ''
	mut_info = ''
	
	whole_induced_muts = len(induced_mutations)
	match_induced_mutations_cnt = 0
	non_induced_mutations_cnt = 0
	induced_mut_type = ''

	for mutation in mutations:

		mut_type = mutation[0]
		pos = mutation[1]
		length = mutation[2]
		insert_info = ''

		if tuple(mutation) in induced_mutations:
			match_induced_mutations_cnt += 1
		else:
			non_induced_mutations_cnt += 1
		if mut_type == 'deletion':
			if length <= largedel_cutlen:
				muts += 'Del,'
				mut_info += f"{pos}_{pos+length-1}:Del_{length},"
			else:
				muts += 'LargeDel,'
				mut_info +=  f"{pos}_{pos+length-1}:LargeDel_{length},"
		elif mut_type == 'insertion':
			length = mutation[2]
			if length <= largeins_cutlen:
				muts += 'Ins,'
				mut_info += f"{pos}_{mutation[4]}:Ins_{length}_{mutation[3]},"
			else:
				muts += 'LargeIns,'
				mut_info += f"{pos}_{mutation[4]}:LargeIns_{length}_{mutation[3]},"
			if len(mutation) == 8:
				for insert in mutation[7]:
					insert_info += f"{insert[0]}_{insert[1]}to{insert[2]}inSeq_{insert[3]}to{insert[4]}inRef_{insert[5]},"
		else: 
			muts += "Sub,"
			mut_info += f"{pos}_{pos+length-1}:Sub_{length}_{mutation[3]}>{mutation[4]},"

	mut_info = mut_info[:-1]
	
	if induced_mutations != []:
		if match_induced_mutations_cnt/whole_induced_muts > partial_induce_cutoff:
			if match_induced_mutations_cnt == whole_induced_muts:
				if non_induced_mutations_cnt == 0:
					induced_mut_type = 'Precise'
				elif non_induced_mutations_cnt > 0:
					induced_mut_type = 'Precise_with_mutations'
			else:
				if non_induced_mutations_cnt == 0:
					induced_mut_type = 'Partial'
				elif non_induced_mutations_cnt > 0:
					induced_mut_type = 'Partial_with_mutations'
		else:
			induced_mut_type = 'X'


	if 'inversion' in insert_info:
		return 'Inv', mut_info, insert_info, induced_mut_type
	elif 'Large' in muts:
		if 'LargeDel' in muts and 'LargeIns' in muts:
			return 'ComplexWithLargeMut', mut_info, insert_info, induced_mut_type
		elif 'LargeDel' in muts:
			return 'LargeDel', mut_info, insert_info, induced_mut_type
		elif 'LargeIns' in muts:
			return 'LargeIns', mut_info, insert_info, induced_mut_type
	elif 'Del' in muts and 'Ins' in muts:
		return 'Complex', mut_info, insert_info, induced_mut_type
	elif 'Del' in muts:
		return 'Del', mut_info, insert_info, induced_mut_type
	elif 'Ins' in muts:
		return 'Ins', mut_info, insert_info, induced_mut_type
	else:
		if 'Sub' in muts:
			return 'Sub', mut_info, insert_info, induced_mut_type
		else:
			return 'Complex', mut_info, insert_info, induced_mut_type

"""def classify_mutations(mutations):
	if not mutations:
		return 'WT', 'None'
	
	def qualifier(type,length):
		if type == 'deletion':
			if length < 100:
				return 'Del'
			else:
				return 'LargeDel'
		elif type == 'insertion':
			if length < 20:
				return 'Ins'
			else:
				return 'LargeIns'
		else: return(type)

	#if len(mutations) > 1:

	return classify_mut_mild(mutations), ','.join(f"{pos}_{pos+length-1}:{qualifier(type,length)}_{length}" for type, pos, length in mutations)
	
	type, pos, length = mutations[0]

	if type == 'deletion':
		if length < 100:
			return 'Del', f"{pos}_{pos+length-1}:Del_{length}"
		else:
			return 'LargeDel', f"{pos}_{pos+flength-1}:LargeDel_{length}"
	elif type == 'insertion':
		if length < 20:
			return 'Ins', f"{pos}_{pos+length-1}:Ins_{length}"
		else:
			return 'LargeIns', f"{pos}_{pos+length-1}:LargeIns_{length}"
	elif type == 'inversion':
		return 'Inv', f"{pos}_{pos+length-1}:Inv_{length}"
	else: return "Sub", f"{pos}_{pos+length-1}:Sub_{length}"
	"""
	

def process_mutations(mutations_dict, output_file, ids, induced_mutations, partial_induce_cutoff, largeins_cutlen, largedel_cutlen):
	with open(output_file, 'w', newline='') as file:
		writer = csv.writer(file, delimiter='\t')
		writer.writerow(['Read_id', 'Classification', 'Mutation_info', 'Integration_info', 'Induce_type', 'whole_mutation'])
		
		for read_id, mutations in mutations_dict.items():

			read_id = ids[read_id]
			classification, mutation_info, insert_info, induce_type = classify_mut_mild(mutations[0], induced_mutations, largeins_cutlen, largedel_cutlen, partial_induce_cutoff=partial_induce_cutoff)
			filtered_mutation_info = []
			if mutations[1] != []:
				tmp_classification, filtered_mutation_info, tmp_insert_info, tmp_induce_type = classify_mut_mild(mutations[1], induced_mutations, largeins_cutlen, largedel_cutlen, partial_induce_cutoff=partial_induce_cutoff)
			if	induce_type == ':':
				induce_type = '-'
			if induced_mutations != []:
				writer.writerow([read_id, classification, mutation_info, insert_info, induce_type, filtered_mutation_info])
			else:
				writer.writerow([read_id, classification, mutation_info, insert_info, '-', filtered_mutation_info])

def write_cnt_file(control_reads_cnt, edited_reads_cnt, output_file):
	with open(output_file, 'w') as fw:
		for f_type in ['Control', 'Treated']:
			for i in control_reads_cnt.keys():
				fw.write(f'{f_type}_{i}\t')
			fw.write('\t')
		fw.write('\n')
		for i in control_reads_cnt.values():
			fw.write(f'{i}\t')
		fw.write('\t')
		for i in edited_reads_cnt.values():
			fw.write(f'{i}\t')



def get_induced_mutation(sam_file_path, fasta_file, cv_pos, cv_pos_2, window, check_window_between_targets, range_align_end, largedel_cutlen, largeins_cutlen, full_window=False):

	reference_sequence = ""
	for record in SeqIO.parse(fasta_file, "fasta"):
		reference_sequence = str(record.seq).upper()
	reference_length = len(reference_sequence)
	samfile = pysam.AlignmentFile(sam_file_path, 'r')

	mutation_information = []

	for read in samfile.fetch():
		mutations_in_read = []
		# Skip unmapped reads
		if read.is_unmapped or read.is_secondary or read.is_supplementary:
			continue
		if not check_in_quality(read):
			continue
		#create more varible way to quantify end

		if read.has_tag('SA'):
			SA_reads = [read]
			SA_n = len(read.get_tag('SA').split(';'))
			ori_query_seq = read.query_sequence.upper()

			while len(SA_reads) < SA_n:
				next_read = next(samfile)
				if next_read.is_supplementary:
					SA_reads.append(next_read)

			if len(SA_reads) != SA_n:
				continue
			SA_reads_muts = analyze_SA_reads(SA_reads, ori_query_seq, reference_sequence, reference_length, read.query_name, range_align_end)

			if len(SA_reads_muts) == 1:
				split_reads_check = True
			else:
				split_reads_check = False

			for split_n, mutations_in_read in enumerate(SA_reads_muts):

				for mut in mutations_in_read:

					mut = list(mut)
					length1 = mut[2]
					pos = mut[1]
					mut_type = mut[0]
				
					window_check = check_in_window(mut, True, cv_pos, cv_pos_2, window, check_window_between_targets)

					if window_check == True: 
						if mut not in mutation_information:
							mutation_information.append(mut)
					else:
						print("ERROR: Induced mutation is outside the window range. Please increase the window option (--window).")
				
		else:

			ref_pos = read.reference_start

			query_seq = read.query_sequence
			query_pos = 0
			
			# Iterate over CIGAR operations
			for operation, length in read.cigar:
				# Check for insertion (I) or deletion (D)
				if operation == 0:	
						sub_tmp_list = []
						for x in range(length):
							"""
							if query_seq[query_pos + x] != reference_sequence[ref_pos + x]:
								if len(sub_tmp_list) > 0 and sub_tmp_list[-1][1] + sub_tmp_list[-1][2] == ref_pos + x:
									sub_tmp_list[-1][2] += 1
									sub_tmp_list[-1][3] += reference_sequence[ref_pos + x]
									sub_tmp_list[-1][4] += query_seq[query_pos + x]
								else:
									sub_tmp_list.append(['substitution', ref_pos + x, 1, reference_sequence[ref_pos + x], query_seq[query_pos + x]])
						for sub in sub_tmp_list:
							mutations_in_read.append(tuple(sub))
							"""
							if query_seq[query_pos + x] != reference_sequence[ref_pos + x]:
								mutations_in_read.append(('substitution', ref_pos+x, 1, reference_sequence[ref_pos + x], query_seq[query_pos + x]))	
				if operation == 1:	# Insertion
					key = ('insertion', ref_pos, length, query_seq[query_pos: query_pos + length], ref_pos + 1, query_pos, query_pos+length)
				elif operation == 2:  # Deletion
					key = ('deletion', ref_pos, length)

				if operation in [0, 2, 3]:	# Match/Mismatch, Deletion, N (Skipped region)
					ref_pos += length
				if operation in [0, 1, 4]:
					query_pos += length

				if operation not in [1,2]: continue
				
				mutations_in_read.append(key)
			
			for mut in mutations_in_read:

				mut = list(mut)
				length1 = mut[2]
				pos = mut[1]
				mut_type = mut[0]
	
				mutation = (mut_type,pos,length1)

				if check_in_window(mut, True, cv_pos, cv_pos_2, window, check_window_between_targets) == True:
					mutation_information.append(tuple(mut))
				else:
					print("ERROR: Induced mutation is outside the window range. Please increase the window option (--window).")
	
	induced_mutation_str = ''
	
	for mutation in mutation_information:
		mut_type = mutation[0]
		pos = mutation[1]
		length = mutation[2]
		insert_info = ''

		if mut_type == 'deletion':
			if length <= largedel_cutlen:
				induced_mutation_str += f"{pos}_{pos+length-1}:Del_{length},"
			else:
				induced_mutation_str +=  f"{pos}_{pos+length-1}:LargeDel_{length},"
		elif mut_type == 'insertion':
			length = mutation[2]
			if length <= largeins_cutlen:
				induced_mutation_str += f"{pos}_{mutation[4]}:Ins_{length}_{mutation[3]},"
			else:
				induced_mutation_str += f"{pos}_{mutation[4]}:LargeIns_{length}_{mutation[3]},"
			if len(mutation) == 8:
				for insert in mutation[7]:
					insert_info += f"{insert[0]}_{insert[1]}to{insert[2]}inSeq_{insert[3]}to{insert[4]}inRef_{insert[5]},"
		else: 
			induced_mutation_str += f"{pos}_{pos+length-1}:Sub_{length}_{mutation[3]}>{mutation[4]},"


	return mutation_information, induced_mutation_str[:-1]



