#!/usr/bin/env python
# coding: utf-8


import argparse, sys
import time, os, pysam
from subprocess import Popen, PIPE

import CRISPRlungo_regular as regular_py
import CRISPRlungo_visualization as visual
import CRISPRlungo_umi
import CRISPRlungo_insert_analysis 


def main():

	start_time = time.time()

	if len(sys.argv) == 1:
		print('usage: CRISPRlungo [regular | umi]')
		sys.exit()

	parser = argparse.ArgumentParser(description='Analysis CRISPR mutaion pattern for Long-read sequencing')
	subparsers = parser.add_subparsers(dest='command')

	regular_parser = subparsers.add_parser('regular', help='Analyze CRISPR mutation for Long-read sequencing data without UMI')
	regular_parser.add_argument('ref', type=str,  help='Reference FASTA file')
	regular_parser.add_argument('control', type=str,  help='Control FASTQ file of Long-read sequencing data')
	regular_parser.add_argument('treated', type=str,  help='Treated FASTQ file of Long-read sequencing data')
	regular_parser.add_argument('output_dir', type=str, help='Output directory name')
	regular_parser.add_argument('target', type=str, help='Target sequence')
	regular_parser.add_argument('--window', type=int, default=10, help='The range for mutation analysis around cleavage site')
	regular_parser.add_argument('-t', '--threads', type=int, default=8, help='Output file name')
	regular_parser.add_argument('--cleavage_pos', type=int, default=16, help='Cleavage position in target sequence [default=16]')
	regular_parser.add_argument('--no_statistical_filter', action='store_false', default=True, help='Disable the statistical filter')
	regular_parser.add_argument('--no_cleavage_site_filter', action='store_false', default=True, help='Disable the cleavage site filter')
	regular_parser.add_argument('--just_visualization', action='store_true', help='If you just want to analyze mutation and the consensus generation is already done using CRISPRlungo, use this option [default: False]')
	regular_parser.add_argument('--additional_target', type=str, default=None, help='Addtional target sequence')
	regular_parser.add_argument('--induced_sequence_path', type=str, default=None, help='FASTA file for induced sequence')
	regular_parser.add_argument('--align_sa_len_threshold', type=int, default=200, help='FASTA file for induced sequence')
	regular_parser.add_argument('--p_value_threshold', type=float, default=0.002, help='Statistical threshold')
	
	regular_parser.add_argument('--whole_window_between_targets', action='store_true', help='Include the region for two targets into window range')
	regular_parser.add_argument('--min_read_cnt', type=int, default=0, help='After counting based on mutation pattern, reads with counts less than the value are removed')
	regular_parser.add_argument('--min_read_freq', type=float, default=0, help='After counting based on mutation pattern, reads with frequency less than the value are removed')
	regular_parser.add_argument('--allele_plot_window', type=float, default=20, help='Window for allele plot, [default: window + 10]')


	


	umi_parser = subparsers.add_parser('umi', help='Analyze CRISPR mutation for Long-read sequencing data with UMI')
	umi_parser.add_argument('ref', type=str,  help='Reference FASTA file')
	umi_parser.add_argument('input_file', type=str,  help='FASTQ file of Long-read sequencing data')
	umi_parser.add_argument('output_dir', type=str,  help='Directory for output files')
	umi_parser.add_argument('target', type=str,  help='Target sequence without PAM')
	umi_parser.add_argument('-t', '--threads', type=int, default=1, help='Output file name')
	umi_parser.add_argument('-c', '--clust_cutoff', type=int, default=5, help='The minimum of UMI cluster size')
	umi_parser.add_argument('--cleavage_pos', type=int, default=16, help='Cleavage position in target sequence')
	umi_parser.add_argument('--window', type=int, default=10, help='The range for mutation analysis around cleavage site')
	umi_parser.add_argument('--mix_tag', type=bool, default=False, help='In mutation classification, when there are multiple mutations, "False" prioritizes indels, labeling them as ins or del, while "True" labels them as Complex. [default: False]')
	umi_parser.add_argument('--just_visualization', type=bool, default=False, help='If you just want to analyze mutation and the consensus generation is already done using CRISPRlungo, use this option [default: False]')
	umi_parser.add_argument('--additional_target', type=str, default=None, help='Addtional target sequence')
	umi_parser.add_argument('--induced_sequence_path', type=str, default=None, help='FASTA file for induced sequence')
	
	umi_parser.add_argument('--whole_window_between_targets', action='store_true', help='Include the region for two targets into window range')
	umi_parser.add_argument('--min_read_cnt', type=int, default=0, help='After counting based on mutation pattern, reads with counts less than the value are removed')
	umi_parser.add_argument('--min_read_freq', type=float, default=0, help='After counting based on mutation pattern, reads with frequency less than the value are removed')
	umi_parser.add_argument('--allele_plot_window', type=float, default=20, help='Window for allele plot, [default: window + 10]')
	umi_parser.add_argument('--align_sa_len_threshold', type=int, default=200, help='FASTA file for induced sequence')


	args = parser.parse_args()
	read_res_cnt = {'fastq': 0, 'filtered': 0, 'aligned': 0}


	def run_minimap2(ref_file, input_file, output_file, longjoin_bandwidth, chaining_bandwidth, threads, len_cutoff=200):
		align_st_time = time.time()
		print(f'minimap2 aligning {input_file} ...\r', end='')
		p = Popen(f'minimap2 -ax map-ont -t {threads} -p 0.5 -r {chaining_bandwidth},{longjoin_bandwidth} {ref_file} {input_file} -o {output_file}', shell=True, stderr=PIPE, stdout=PIPE).communicate()
		if p[1].decode('utf-8').find('ERROR') != -1:
			print(p[1].decode('utf-8'))
			sys.exit()
		else:
			print(f'minimap2 aligning {input_file} ... Done {round(time.time() - align_st_time, 2)} s')			
	
	def soft_clipped(align_file_1, output_file, cut_len_threshold, fasta_check=False):
		print('with', cut_len_threshold)
		check_SA = False
		read_d = []
		with pysam.AlignmentFile(align_file_1) as f1, open(output_file, 'w') as fw:
			for line in f1:
				if line.is_unmapped or line.is_secondary or line.is_supplementary:
					continue
				clipped_seq = [[],[]]
				cigar = line.cigar
				read_d.append(line)
				if 'alignStrandInfo' not in line.query_name:
					line.query_name += f'_alignStrandInfo'
				
				if line.is_reverse:
					strand = -1
				else:
					strand = 1
				line.query_name += f'_{strand}'
				if fasta_check == False:
					if cigar[0][0] == 4 and cigar[0][1] > cut_len_threshold:
						fw.write(f'@{line.query_name}\n{line.query_sequence[:cigar[0][1]]}\n+\n{line.qual[:cigar[0][1]]}\n')
						check_SA = True
					if cigar[-1][0] == 4 and cigar[-1][1] > cut_len_threshold:
						fw.write(f'@{line.query_name}\n{line.query_sequence[-cigar[-1][1]:]}\n+\n{line.qual[-cigar[-1][1]:]}\n')
						check_SA = True
				else:
					if cigar[0][0] == 4 and cigar[0][1] > cut_len_threshold:
						fw.write(f'>{line.query_name}\n{line.query_sequence[:cigar[0][1]]}\n')
						check_SA = True
					if cigar[-1][0] == 4 and cigar[-1][1] > cut_len_threshold:
						fw.write(f'>{line.query_name}\n{line.query_sequence[-cigar[-1][1]:]}\n')
						check_SA = True

				
		return check_SA, read_d

	def run_triple_minimap2(ref_file, input_file, output_file, longjoin_bandwidth, chaining_bandwidth, threads, len_cutoff = 200, fasta_check=False):
		output_path = output_file[:output_file.rfind('/')]
		check_SA = True
		n = 0
		read_d = {}
		file_format = 'fastq'
		if fasta_check:
			file_format = 'fasta'
		while check_SA == True:
			n += 1
			if n == 1:
				fastq_file = input_file
			else:
				fastq_file = output_path + f'/{n}_soft.{file_format}'
			print(output_path)
			print(output_file)
			run_minimap2(ref_file, fastq_file, output_path + f'/{n}_align.sam' , longjoin_bandwidth, chaining_bandwidth, threads)
			check_SA, read_out = soft_clipped(output_path + f'/{n}_align.sam', output_path + f'/{n+1}_soft.{file_format}', len_cutoff, fasta_check)
			for read in read_out:
				query_name = read.query_name
				strand_info = query_name.split('_alignStrandInfo_')[1].split('_')
				strand = 1
				for x in strand_info:
					strand *= int(x)
				if strand == -1:
					read.is_reverse = True
				else:
					read.is_reverse = False
				query_name = query_name[:query_name.find('_alignStrandInfo_')]
				if n != 1:
					read.is_supplementary = True
				if query_name not in read_d:
					read_d[query_name] = [read]
				else:
					ori_seq = read_d[query_name][0].query_sequence
					ori_strand = -1 if read_d[query_name][0].is_reverse else 1
					if ori_strand * strand == -1:
						part_st = ori_seq.find(reverse_complementary(read.query_alignment_sequence))
						part_ed = ori_seq.find(reverse_complementary(read.query_alignment_sequence)) + len(read.query_alignment_sequence)
					else:
						part_st = ori_seq.find(read.query_alignment_sequence)
						part_ed = ori_seq.find(read.query_alignment_sequence) + len(read.query_alignment_sequence)
					
					if not fasta_check: q = read.query_qualities

					if read.cigar[0][0] == 4:
						read.query_sequence = read.query_alignment_sequence
						if not fasta_check: q = q[read.cigar[0][1]:]
						read.cigar = [(5, part_st)] + read.cigar[1:]
						
					else:
						read.cigar = [(5, part_st)] + read.cigar

					if read.cigar[-1][0] == 4:
						read.query_sequence = read.query_alignment_sequence
						if not fasta_check: q = q[:-read.cigar[-1][1]]
						read.cigar = read.cigar[:-1] + [(5, len(ori_seq) - part_ed)]
					else:
						read.cigar += [(5, len(ori_seq) - part_ed)]
					if not fasta_check: read.query_qualities = q

					read_d[query_name].append(read)
		fw = open(output_file, 'w') 
		with open(output_path + '/1_align.sam') as f:
			for line in f:
				if line[0] == '@':
					fw.write(line)
				else:
					break
		for x, reads in read_d.items():
			read_n = 0
			if len(reads) == 1:
				tags = reads[0].get_tags()
				filtered_tags = [t for t in tags if t[0] != "SA"]
				reads[0].set_tags(filtered_tags)
			else:
				reads[0].set_tag("SA", ('N;'*(len(reads)))[:-1], value_type='Z')
			for read in reads:
				read_n += 1
				fw.write(read.to_string() +f'\t\n')
		fw.close()
	
	
	def create_dir(dir_name):
		if not os.path.exists(dir_name):
			os.makedirs(dir_name)

	def reverse_complementary(s):
		return s.translate(s.maketrans('ATGCatgc','TACGtacg'))[::-1]
	
	def rev_comp(s):
		return s.translate(s.maketrans('ATGC','TACG'))[::-1]	

	output_dir = args.output_dir
	create_dir(output_dir)
	create_dir(f'{output_dir}/results')

	if args.additional_target == None and args.whole_window_between_targets == True:
		print('ERROR: whole_window_between_targets is declared, but additional_target is not declared.')
		sys.exit()

	if args.allele_plot_window == 0:
		allele_plot_window = args.window + 10
	else:
		allele_plot_window = args.allele_plot_window
	
	ref_file = args.ref
	ref_file_list = open(args.ref).readlines()
	ref_seq = ''.join(ref_file_list[1:]).replace('\n','').upper()
	ref_name = ref_file_list[0].strip()[1:]
	ref_len = len(ref_seq)

	original_target = 1
	
	if args.target != None:
		target = args.target.upper()
		if ref_seq.find(target) != -1:
			cv_pos = ref_seq.find(target) + args.cleavage_pos
			strand = 1
		elif ref_seq.find(reverse_complementary(target)) != -1:
			cv_pos = ref_seq.find(reverse_complementary(target)) + len(target) - args.cleavage_pos - 2
			strand = -1
		else:
			print('ERROR: Can not find target sequence in reference !!')
			sys.exit()

	if args.additional_target != None:
		target_2 = args.additional_target.upper()
		if ref_seq.find(target_2) != -1:
			cv_pos_2 = ref_seq.find(target_2) + args.cleavage_pos
			strand_2 = 1
		elif ref_seq.find(reverse_complementary(target_2)) != -1:
			cv_pos_2 = ref_seq.find(reverse_complementary(target_2)) + len(target) - args.cleavage_pos - 2
			strand_2 = -1
		else:
			print('ERROR: Can not find additional targt sequence in reference !!')
			sys.exit()
	else:
		target_2 = False
		cv_pos_2 = False
		strand_2 = False

	if cv_pos_2 != False and cv_pos > cv_pos_2:
		tmp_pos = [cv_pos, target, strand]
		cv_pos, target, strand = [cv_pos_2, target_2, strand_2]
		cv_pos_2, target_2, strand_2 = tmp_pos
		original_target = 2

	longjoin_bandwidth = int(ref_len * 0.3)
	if longjoin_bandwidth < 500:
		chaining_bandwidth = longjoin_bandwidth
	else:
		chaining_bandwidth = 500

	if args.command == 'regular':
		
		print('Start to analysis using REGULAR mode ...')
			
		threads = args.threads
		
		# Run Miniseq2 to align received files
		
		

		if args.just_visualization == False:

			ref_indexed_file = os.path.splitext(ref_file)[0] + '.mmi'
			Popen(f'minimap2 -d {ref_indexed_file} {ref_file}', shell=True, stderr=PIPE, stdout=PIPE).communicate()

			if args.induced_sequence_path:
				run_minimap2(ref_indexed_file, args.induced_sequence_path,  output_dir + '/induced_mutation_reference.sam', 20000, 500, 1, len_cutoff=args.align_sa_len_threshold)
				induced_mutations = regular_py.get_induced_mutation(output_dir + '/induced_mutation_reference.sam', ref_file, cv_pos, cv_pos_2, args.window, args.whole_window_between_targets)
			else:
				induced_mutations = []
			

			run_minimap2(ref_indexed_file, args.control,  output_dir + '/Control_alignment.sam', longjoin_bandwidth, chaining_bandwidth, threads, len_cutoff=args.align_sa_len_threshold)
			run_minimap2(ref_indexed_file, args.treated,  output_dir + '/Treated_alignment.sam', longjoin_bandwidth, chaining_bandwidth, threads, len_cutoff=args.align_sa_len_threshold)

			# Run Statistical anlaysis.py
			print(args.no_statistical_filter)
			edited_dictionary, controll_dictionary, List_of_valid_IDs, control_reads_cnt, edited_reads_cnt = regular_py.analysis_function(output_dir + '/Control_alignment.sam', 
				output_dir + '/Treated_alignment.sam', 
				ref_file, 
				cv_pos, 
				cv_pos_2, 
				args.window, 
				args.whole_window_between_targets, 
				induced_mutations, 
				threads=args.threads, 
				Filter1=args.no_statistical_filter, 
				window_filter=args.no_cleavage_site_filter, 
				p_limit_value=args.p_value_threshold, allowance_value=0)
			
			c = 0
			for read_id, mutations in edited_dictionary.items():
				read_id = List_of_valid_IDs[read_id]
				classification, mutation_info, insert_info, induce_type = regular_py.classify_mut_mild(mutations, induced_mutations)
				if read_id.find('Inv_3000') != -1:
					c += 1

			print('\nbefore count', c)
			print('\n\n')

			edited_dictionary = CRISPRlungo_insert_analysis.confirm_insertion_seq(edited_dictionary, ref_seq, ref_name, './possible_insertion.fasta', output_dir, threads)

			c = 0
			for read_id, mutations in edited_dictionary.items():
				read_id = List_of_valid_IDs[read_id]
				classification, mutation_info, insert_info, induce_type = regular_py.classify_mut_mild(mutations, induced_mutations)
				if read_id.find('Inv_3000') != -1:
					c += 1

			print('\nafter count', c)
			print('\n\n')

			regular_py.process_mutations(edited_dictionary, f'{output_dir}/results/read_classification.txt', List_of_valid_IDs, induced_mutations)
			regular_py.write_cnt_file(control_reads_cnt, edited_reads_cnt, f'{output_dir}/results/preprocess_count.txt')


	elif args.command == 'umi':

		print('Start to analyze using UMI mode ...')
		
		ref_file = args.ref
		input_file = args.input_file
		threads = args.threads

		# Get UMI position information
		umi_pos, umi_len, ref_name, ref_seq = CRISPRlungo_umi.input_organize(ref_file, output_dir,input_file)
		ref_len = len(ref_seq)
		
		index_info = 'result,NNNNNNNNNN' #For index demultiplexing, not yet supported.
		
		index_info_list = index_info.split(',')
		index_names = []
		for i in range(int(len(index_info_list)/2)):
			i = index_info_list[2*i]
			index_names.append(i)

		if args.just_visualization == False:
	
			ref_dir = f'{output_dir}/ref_seq'
			create_dir(ref_dir)
		
			with open(ref_dir + '/ref.fasta', 'w') as f:
				f.write(f'>{ref_name}\n{ref_seq}')
		
			create_dir(output_dir + '/align')

			# Align recieved FASTQ file
			longjoin_bandwidth = int(ref_len * 0.3)
			if longjoin_bandwidth < 500:
				chaining_bandwidth = longjoin_bandwidth
			else:
				chaining_bandwidth = 500
				
			p = Popen(f'minimap2 -d {ref_dir}/ref.mmi {ref_dir}/ref.fasta', shell=True, stderr=PIPE, stdout=PIPE).communicate()
			if args.induced_sequence_path:
				run_minimap2(f'{ref_dir}/ref.mmi', args.induced_sequence_path,  output_dir + '/induced_mutation_reference.sam', 20000, 500, 1, len_cutoff=args.align_sa_len_threshold)
				induced_mutations = regular_py.get_induced_mutation(output_dir + '/induced_mutation_reference.sam', ref_file, cv_pos, cv_pos_2, args.window, args.whole_window_between_targets)
			else:
				induced_mutations = []
				
			run_minimap2(ref_dir + "/ref.mmi", args.input_file,  output_dir + '/align/input_fastq_align.sam', 20000, 500, threads, len_cutoff=args.align_sa_len_threshold)

			# Get UMI sequence
			create_dir(f'{output_dir}/demultiplexing')

			CRISPRlungo_umi.extract_index_umi(output_dir + '/align/input_fastq_align.sam', ref_file, output_dir, umi_pos[0], umi_pos[1])

			# Clustering UMI sequence

			create_dir(f'{output_dir}/clustering')
			create_dir(f'{output_dir}/consensus')

			for i in index_names:
			
				create_dir(f'{output_dir}/clustering/{i}')
				create_dir(f'{output_dir}/clustering/{i}/1st_clusters')
				create_dir(f'{output_dir}/clustering/{i}/2nd_clusters')
				create_dir(f'{output_dir}/clustering/{i}/medaka_input')
				create_dir(f'{output_dir}/consensus/{i}')
			
			for fn in index_names:
				CRISPRlungo_umi.clustering_umi( f'{output_dir}/demultiplexing/umi_{fn}.fasta', umi_len, f'{output_dir}/clustering/{fn}', f'{output_dir}/consensus/{fn}', threads, args.clust_cutoff)
				run_minimap2(ref_dir + "/ref.mmi", f'{output_dir}/consensus/{fn}/consensus.fasta', f'{output_dir}/align/consensus_result.sam', 20000, 500, threads, fasta_check=True, len_cutoff=args.align_sa_len_threshold)

		# Mutation Analysis
		print('Mutation analysis...')

		for fn in index_names:

			#create_dir(f'{output_dir}/results/{fn}')
			CRISPRlungo_umi.mutation_analysis(ref_seq, ref_name, cv_pos, strand, cv_pos_2, strand_2, args.window, f'{output_dir}/align/consensus_{fn}.sam', f'{output_dir}/results', args.whole_window_between_targets, induced_mutations, threads=threads, mix_tag=args.mix_tag)

	#visualization

	full_html = """
	<!DOCTYPE html>
	<html>
	<head>
		<title>Combined Graphs</title>
		<style>
			body {margin: 0;padding: 0;display: flex;justify-content: center;align-items: center;background-color: #f4f4f4;}
    		#graph-container {min-width: 500px; max-width: 1000px;width: 80vw;background-color: #ffffff;border: 1px solid #ddd;box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);overflow: hidden; }
			img {max-width: 100%;max-height: 100%;object-fit: contain;}
		</style>
	</head>
	<body>
		<div id='graph-container'>"""

	print('Drawing graphs ... \r', end='')

	tsv_file = output_dir + '/results/read_classification.txt'
	graph_output_dir = output_dir + '/results/'

	plots = {}

	if args.command == 'regular':
		read_per_position = visual.visualization_preprocess_regular(output_dir + '/Treated_alignment.sam', ref_file)
		plot_html = visual.align_count_plot(f'{output_dir}/results/preprocess_count.txt', f'{output_dir}/results')
		plots['treated_align'] = plot_html[0]
		plots['control_align'] = plot_html[1]
		visual.regular_accuracy_plot(ref_seq, read_per_position, graph_output_dir)
	else:
		read_per_position = visual.visualization_preprocess_regular(output_dir + '/align/consensus_result.sam', ref_file)
		plot_html = visual.align_count_plot(f'{output_dir}/results/preprocess_count.txt', f'{output_dir}/results')
		plots['treated_align'] = plot_html[0]
		plots['control_align'] = ''
	
	read_cnt_file = f'{output_dir}/results/mutation_patter_count.txt'
	visual.write_read_count(tsv_file,  f'{output_dir}/results/preprocess_count.txt', read_cnt_file, f'{output_dir}/results/mutation_summary_count.txt', args.min_read_cnt, args.min_read_freq, args.induced_sequence_path)
	if args.target != None:
		print('Drawing graphs ... Allele plot                         \r', end='')
		visual.allele_plot(ref_seq, cv_pos, cv_pos_2, strand, strand_2, read_cnt_file, graph_output_dir, args.cleavage_pos, target, target_2, original_target, args.min_read_cnt, args.min_read_freq, allele_plot_window)
	print('Drawing graphs ... pie plot                         \r', end='')
	plots['mutation_pie'], plots['pattern_pie'], plots['allele_pie'] = visual.mutation_pie_chart(read_cnt_file, graph_output_dir)
	print('Drawing graphs ... indel plot                         \r', end='')
	plots['indel_per_pos'] = visual.indel_per_position(read_cnt_file, ref_seq, graph_output_dir)
	print('Drawing graphs ... insertion plot                         \r', end='')
	plots['insertion_len'] = visual.Insertion_length(read_cnt_file, graph_output_dir)
	print('Drawing graphs ... deletion plot                         \r', end='')
	plots['deletion_len'] = visual.Deletion_length(read_cnt_file, graph_output_dir)
	plots['deletion_count_len'] = visual.Deletion_count_length(read_cnt_file, graph_output_dir)	
	print('Drawing graphs ... large deletion tornado plot                         \r', end='')
	visual.LD_tornado(read_cnt_file, cv_pos, ref_len, strand, graph_output_dir)
	print('Drawing graphs ... base proportion plot                         \r', end='')
	plots['base_proportion'] = visual.base_proportion(read_per_position, graph_output_dir, ref_seq)
	
	visual.write_html(plots, args.command, args.target, output_dir)

if __name__=='__main__':
	main()

