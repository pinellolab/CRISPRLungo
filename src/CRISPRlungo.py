#!/usr/bin/env python

import argparse, sys
import time, os, pysam
from subprocess import Popen, PIPE

import CRISPRlungo_regular_new as regular_py
import CRISPRlungo_visualization as visual
import CRISPRlungo_umi
import CRISPRlungo_insert_analysis 
from CRISPRlungo_minimap import *



def main():

	start_time = time.time()
	
	parser = argparse.ArgumentParser(description='Analysis CRISPR mutaion pattern for Long-read sequencing')

	parser.add_argument('ref', type=str,  help='Reference FASTA file, if your file has UMI, the reference includes UMI sequence and marks UMI using ( and ).')
	parser.add_argument('treated', type=str,  help='Treated FASTQ file of Long-read sequencing data')
	parser.add_argument('output_dir', type=str, help='Output directory name')
	parser.add_argument('target', type=str, help='Target sequence')

	parser.add_argument('--umi', action='store_true', help='If UMI sequencing is used, use this option.')

	parser.add_argument('--window', type=int, default=5, help='The range for mutation analysis around cleavage site')
	parser.add_argument('--whole_window_between_targets', action='store_true', help='Include the region for two targets into window range')
	parser.add_argument('--cleavage_pos', type=int, default=16, help='Cleavage position in target sequence [default=16]')
	parser.add_argument('--additional_target', type=str, default=None, help='Addtional target sequence')
	parser.add_argument('--control', type=str, default=None, help='When a control file is input, background filtering is performed using the control file')
	parser.add_argument('--induced_sequence_path', type=str, default=None, help='When a desired sequence file is input, additional classification of desired mutations is performed')
	parser.add_argument('--integration_file', type=str, default=False, help='FASTA file for possible sequences that can be integrated')
	parser.add_argument('--merge_substitution', action='store_true', help='A continuous substitution is considered as one mutation.')

	parser.add_argument('--largeins_cutlen', default=20, type=int, help='The minimum length for large deletions')
	parser.add_argument('--largedel_cutlen', default=100, type=int, help='The minimum length for large insertion')

	parser.add_argument('--min_read_cnt', type=int, default=0, help='After counting based on mutation pattern, reads with counts less than the value are removed')
	parser.add_argument('--min_read_freq', type=float, default=0, help='After counting based on mutation pattern, reads with frequency less than the value are removed')	
	parser.add_argument('--mix_tag', type=bool, default=False, help='In mutation classification, when there are multiple mutations, "False" prioritizes indels, labeling them as ins or del, while "True" labels them as Complex. [default: False]')
	parser.add_argument('--min_mut_freq_no_control_refmut', type=float, default = 0.5, help='In mutation classification, when there are multiple mutations, "False" prioritizes indels, labeling them as ins or del, while "True" labels them as Complex. [default: False]')
	parser.add_argument('--induced_paritial_similiarity', type=float, default=0.8, help='The If the mutation pattern similarity is higher than this value, but not completely identical, it is considered partially induced.')
	parser.add_argument('--range_both_end_region', type=int, default=100, help='If the reads were not aligned this range from both end, the read is considered as short fragment.')

	parser.add_argument('--align_sa_len_threshold', type=int, default=100, help='FASTA file for induced sequence')
	parser.add_argument('--p_value_threshold', type=float, default=0.002, help='Statistical threshold')
	parser.add_argument('--mut_freq_threshold', type=float, default=0, help='muation frequency threshold, if you want more harsh filteration, you can use this')
	parser.add_argument('-c', '--clust_cutoff', type=int, default=5, help='The minimum of UMI cluster size')

	parser.add_argument('--just_visualization', action='store_true', help='If you just want to analyze mutation and the consensus generation is already done using CRISPRlungo, use this option [default: False]')
	parser.add_argument('--allele_plot_window', type=int, default=20, help='Window for allele plot, [default: window + 10]')
	parser.add_argument('--allele_plot_lines', type=int, default=20, help='Window for allele plot, [default: window + 10]')
	parser.add_argument('--show_all_between_allele', action='store_true', help='Draw all sequences between two targets in an allele plot')

	parser.add_argument('-t', '--threads', type=int, default=8, help='Output file name')

	args = parser.parse_args()
	read_res_cnt = {'fastq': 0, 'filtered': 0, 'aligned': 0}
	
	current_dir = os.path.dirname(os.path.abspath(__file__))

	def create_dir(dir_name):
		if not os.path.exists(dir_name):
			os.makedirs(dir_name)


	output_dir = args.output_dir
	threads = args.threads
	treated_file_path = args.treated
	if args.control:
		control_file_path = args.control
	


	create_dir(output_dir)
	create_dir(f'{output_dir}/results')
	create_dir(f'{output_dir}/custom_results')

	if args.additional_target == None and args.whole_window_between_targets == True:
		print('ERROR: whole_window_between_targets is declared, but additional_target is not declared.')
		sys.exit()

	if args.allele_plot_window == 0:
		allele_plot_window = args.window + 10
	else:
		allele_plot_window = args.allele_plot_window

	# Get reference information 

	ref_file = args.ref

	ref_dir = f'{output_dir}/ref_seq'
	create_dir(ref_dir)

	range_align_end = args.range_both_end_region

	
	if args.umi:

		umi_pos, umi_len, ref_name, ref_seq = CRISPRlungo_umi.input_organize(ref_file, output_dir, treated_file_path)
		ref_len= len(ref_seq)
		
		ref_seq_w_umi = ref_seq
		ref_seq = ref_seq[umi_pos[0][1]+1: umi_pos[1][0]]

		fw = open(output_dir + '/ref_seq/ref_wo_umi.fasta', 'w')
		fw.write(f'>{ref_name}\n{ref_seq}\n')
		fw.close()

		ref_file_w_umi = ref_dir + '/ref.fasta'
		p = Popen(f'minimap2 -d {output_dir}/ref_seq/ref_wo_umi.mmi {output_dir}/ref_seq/ref_wo_umi.fasta', shell=True, stderr=PIPE, stdout=PIPE).communicate()

		fw = open(output_dir + '/ref_seq/ref.fasta', 'w')
		fw.write(f'>{ref_name}\n{ref_seq_w_umi}\n')
		fw.close()
		p = Popen(f'minimap2 -d {output_dir}/ref_seq/ref.mmi {output_dir}/ref_seq/ref.fasta', shell=True, stderr=PIPE, stdout=PIPE).communicate()

	else:

		ref_file_list = open(args.ref).readlines()
		ref_seq = ''.join(ref_file_list[1:]).replace('\n','').upper().replace('(', '')
		ref_name = ref_file_list[0].strip()[1:]
		ref_len = len(ref_seq)

		fw = open(ref_dir + '/ref_wo_umi.fasta', 'w')
		fw.write(f'>{ref_name}\n{ref_seq}\n')
		fw.close()

		p = Popen(f'minimap2 -d {ref_dir}/ref_wo_umi.mmi {ref_dir}/ref_wo_umi.fasta', shell=True, stderr=PIPE, stdout=PIPE).communicate()
		ref_index_path = f'{ref_dir}/ref_wo_umi.mmi'
		

	original_target = 1
	
	final_mutation_analysis_file = args.treated

	# Get target sequence information

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


	# get induced mutation patterns

	if args.induced_sequence_path:
		run_triple_minimap2(f'{ref_dir}/ref_wo_umi.mmi', args.induced_sequence_path,  output_dir + '/induced_mutation_reference.sam', longjoin_bandwidth, chaining_bandwidth, 1, len_cutoff=args.align_sa_len_threshold)
		induced_mutations, induced_mutation_str = regular_py.get_induced_mutation(output_dir + '/induced_mutation_reference.sam', ref_file, cv_pos, cv_pos_2, args.window, args.whole_window_between_targets, range_align_end, args.largedel_cutlen, args.largeins_cutlen)
	else:
		induced_mutations = []
		induced_mutation_str = False

	create_dir(output_dir + '/align')

	# UMI Clustering
	

	if args.umi:

		print('Generating treated consensus file...')

		index_info = 'result,NNNNNNNNNN' #For index demultiplexing, not yet supported.

		index_info_list = index_info.split(',')
		index_names = []
		for i in range(int(len(index_info_list)/2)):
			i = index_info_list[2*i]
			index_names.append(i)

		if args.just_visualization == False:

			# Align recieved FASTQ file
			run_triple_minimap2(ref_dir + "/ref.mmi", treated_file_path,  output_dir + '/align/input_fastq_align_to_umiref.sam', longjoin_bandwidth, chaining_bandwidth, threads, len_cutoff=args.align_sa_len_threshold)

			# Get UMI sequence
			create_dir(f'{output_dir}/demultiplexing')

			CRISPRlungo_umi.extract_index_umi(output_dir + '/align/input_fastq_align_to_umiref.sam', ref_file_w_umi, output_dir, umi_pos[0], umi_pos[1], index_information='treated,NNNNNNNNNN')

			# Clustering UMI sequence

			create_dir(f'{output_dir}/clustering')
			create_dir(f'{output_dir}/consensus')

			create_dir(f'{output_dir}/clustering/treated')
			create_dir(f'{output_dir}/clustering/treated/1st_clusters')
			create_dir(f'{output_dir}/clustering/treated/2nd_clusters')
			create_dir(f'{output_dir}/clustering/treated/medaka_input')
			create_dir(f'{output_dir}/consensus/treated')

			CRISPRlungo_umi.clustering_umi( f'{output_dir}/demultiplexing/umi_treated.fasta', umi_len, f'{output_dir}/clustering/treated', f'{output_dir}/consensus/treated', threads, args.clust_cutoff)

			treated_file_path = f'{output_dir}/consensus/treated/consensus.fasta'

			if args.control:
				
				print('Generating control consensus file...')

				run_triple_minimap2(ref_dir + "/ref.mmi", args.control,  output_dir + '/align/control_fastq_align_to_umiref.sam', longjoin_bandwidth, chaining_bandwidth, threads, len_cutoff=args.align_sa_len_threshold)

				# Get UMI sequence

				CRISPRlungo_umi.extract_index_umi(output_dir + '/align/control_fastq_align_to_umiref.sam', ref_file_w_umi, output_dir, umi_pos[0], umi_pos[1], index_information='control,NNNNNNNNNN')

				# Clustering UMI sequence
				
				create_dir(f'{output_dir}/clustering/control')
				create_dir(f'{output_dir}/clustering/control/1st_clusters')
				create_dir(f'{output_dir}/clustering/control/2nd_clusters')
				create_dir(f'{output_dir}/clustering/control/medaka_input')
				create_dir(f'{output_dir}/consensus/control')

				CRISPRlungo_umi.clustering_umi( f'{output_dir}/demultiplexing/umi_control.fasta', umi_len, f'{output_dir}/clustering/control', f'{output_dir}/consensus/control', threads, args.clust_cutoff)

				control_file_path = f'{output_dir}/consensus/control/consensus.fasta'
			
			ref_index_path = f'{ref_dir}/ref_wo_umi.mmi'
			final_mutation_analysis_file = f'{output_dir}/consensus/treated/consensus.fasta'
			
	create_dir(f'{output_dir}/results/')
	create_dir(f'{output_dir}/css/')

	if args.control:

		print('Start filtering module ...')

		if args.just_visualization == False:
			
			run_triple_minimap2(f'{ref_dir}/ref_wo_umi.mmi', control_file_path,  output_dir + '/align/Control_alignment.sam', longjoin_bandwidth, chaining_bandwidth, threads, len_cutoff=args.align_sa_len_threshold, fasta_check=True)
			run_triple_minimap2(f'{ref_dir}/ref_wo_umi.mmi', treated_file_path,  output_dir + '/align/Treated_alignment.sam', longjoin_bandwidth, chaining_bandwidth, threads, len_cutoff=args.align_sa_len_threshold, fasta_check=True)

			# Run Statistical anlaysis.py

			edited_dictionary, controll_dictionary, List_of_valid_IDs, control_reads_cnt, edited_reads_cnt = regular_py.analysis_function(
				output_dir + '/align/Control_alignment.sam', 
				output_dir + '/align/Treated_alignment.sam', 
				f'{ref_dir}/ref_wo_umi.fasta', 
				output_dir + '/results/',
				cv_pos, 
				cv_pos_2, 
				args.window, 
				args.whole_window_between_targets, 
				induced_mutations, 
				range_align_end = range_align_end,
				threads=threads, 
				largeins_cutlen = args.largeins_cutlen,
				largedel_cutlen = args.largedel_cutlen,
				p_limit_value=args.p_value_threshold, 
				mut_freq_value = args.mut_freq_threshold,
				allowance_value=0,
				umi_clustered = args.umi)

			edited_dictionary = CRISPRlungo_insert_analysis.confirm_insertion_seq(edited_dictionary, ref_seq, ref_name, args.integration_file, output_dir, threads)

			regular_py.process_mutations(edited_dictionary, 
				f'{output_dir}/results/read_classification.txt', 
				List_of_valid_IDs, 
				induced_mutations, 
				args.induced_paritial_similiarity, args.largeins_cutlen, args.largedel_cutlen,)
			regular_py.write_cnt_file(control_reads_cnt, edited_reads_cnt, f'{output_dir}/results/preprocess_count.txt')
	
	else:
		ref_index = ref_dir + '/ref_wo_umi.mmi'
		run_triple_minimap2(ref_index, final_mutation_analysis_file, f'{output_dir}/align/Treated_alignment.sam', longjoin_bandwidth, chaining_bandwidth, threads, len_cutoff=args.align_sa_len_threshold, fasta_check=True)
		if not args.umi and not args.control:
			write_cnt = True
		else:
			write_cnt = False
		CRISPRlungo_umi.mutation_analysis(ref_seq, ref_name, cv_pos, strand, cv_pos_2, strand_2, args.window, 
			f'{output_dir}/align/Treated_alignment.sam', 
			f'{output_dir}/results', 
			args.whole_window_between_targets, induced_mutations, current_dir,
			threads=threads, mix_tag=args.mix_tag, 
			write_cnt = write_cnt, partial_induce_cutoff=args.induced_paritial_similiarity, 
			range_align_end=range_align_end, 
			largeins_cutlen = args.largeins_cutlen,
			largedel_cutlen = args.largedel_cutlen)


	#Visualization

	if args.just_visualization:
		edited_reads_cnt = {}
		f = open(f'{output_dir}/results/preprocess_count.txt').readlines()
		menu = f[0].strip().split()
		val = f[1].strip().split()
		for i in zip(menu, val):
			i = list(i)
			if 'Treated' in i[0]:
				i[0] = '_'.join(i[0].split('_')[1:])
				edited_reads_cnt[i[0]] = int(i[1])

	print('Drawing graphs ... \r', end='')

	tsv_file = output_dir + '/results/read_classification.txt'
	graph_output_dir = output_dir + '/results/'

	plots = {}

	read_per_position = visual.visualization_preprocess_regular(output_dir + '/align/Treated_alignment.sam', ref_file)

	if args.control:
		visual.regular_accuracy_plot(ref_seq, read_per_position, graph_output_dir)
	
	fw_input = open(f'{output_dir}/results/input_summary.txt', 'w')
	fw_input.write(f'Target_1 :{target}\n')
	fw_input.write(f'Target_2 :{target_2}\n')
	fw_input.write(f'Ref_seq :{ref_seq}\n')
	fw_input.write(f'CleavagePos_1 :{cv_pos}\n')
	fw_input.write(f'CleavageStrand_1 :{strand}\n')
	fw_input.write(f'CleavagePos_2 :{cv_pos_2}\n')
	fw_input.write(f'CleavageStrand_2 :{strand_2}\n')
	fw_input.write(f'Window : {args.window}\n')
	fw_input.write(f'Window_between_cleavage : {args.whole_window_between_targets}\n')
	fw_input.write(f'induced_mut :{induced_mutation_str}\n')
	fw_input.write(f'cut_pos_in_target :{args.cleavage_pos}\n')
	fw_input.write(f'original_target :{original_target}\n')
	fw_input.write(f'minimum_read_count : {args.min_read_cnt}\n')
	fw_input.write(f'minimum_read_frequency : {args.min_read_freq}\n')
	fw_input.write(f'induced_sequence_path : {args.induced_sequence_path}\n')
	fw_input.write(f'largeins_cutlen : {args.largeins_cutlen}\n')
	fw_input.write(f'largedel_cutlen : {args.largedel_cutlen}\n')
	fw_input.close()

	read_cnt_file = f'{output_dir}/results/mutation_patter_count.txt'
	mut_cnt, precise_cnt = visual.write_read_count(tsv_file,  f'{output_dir}/results/preprocess_count.txt', read_cnt_file, f'{output_dir}/results/mutation_summary_count.txt', args.min_read_cnt, args.min_read_freq, args.induced_sequence_path)
	plot_html = visual.align_count_plot(f'{output_dir}/results/preprocess_count.txt', f'{output_dir}/results/mutation_summary_count.txt', f'{output_dir}/results')
	plots['treated_align'] = plot_html[1]
	plots['control_align'] = ''
	if args.control:
		plots['control_align'] = plot_html[0]
	#read_cnt_file = f'{output_dir}/results/mutation_patter_count.txt'
	if args.target != None:
		print('Drawing graphs ... Allele plot                         \r', end='')
		visual.allele_plot(ref_seq, cv_pos, cv_pos_2, strand, strand_2, read_cnt_file, graph_output_dir, args.cleavage_pos, target, target_2, original_target, args.min_read_cnt, args.min_read_freq, allele_plot_window, args.allele_plot_lines, induced_mutation_str, args.show_all_between_allele)
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
	plots['base_proportion'] = visual.base_proportion(read_per_position, graph_output_dir, ref_seq, cv_pos, cv_pos_2, allele_plot_window, args.show_all_between_allele)
	
	visual.write_html(plots, args.control, args.target, output_dir, mut_cnt, precise_cnt, edited_reads_cnt)


if __name__=='__main__':
	main()
