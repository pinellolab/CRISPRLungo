import time, pysam, sys, os
from subprocess import Popen, PIPE

import CRISPRlungo_log as log


def reverse_complementary(s):
	return s.translate(s.maketrans('ATGCatgc','TACGtacg'))[::-1]

def run_minimap2(ref_file, input_file, output_file, longjoin_bandwidth, chaining_bandwidth, threads, minimap2_opt, task=None):
		# `task` lets a caller (run_triple_minimap2) fold every pass into one
		# line of output instead of printing one line per pass.
		own_task = task is None
		if own_task:
			task = log.task(f'minimap2 : {os.path.basename(input_file)}')

		cmd = f'minimap2 -ax map-ont -t {threads} -z 100 -p 0.5 -r {chaining_bandwidth},{longjoin_bandwidth} {minimap2_opt} {ref_file} {input_file} -o {output_file}'
		proc = Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
		_, stderr = proc.communicate()
		stderr_text = stderr.decode('utf-8', errors='replace')

		# Check the exit status, not just the text of stderr. Scanning for the
		# word "ERROR" alone let silent failures (minimap2 missing from PATH, a
		# killed process, a write that never happened) through, and the run then
		# died much later with a confusing "no such file" from pysam.
		if proc.returncode != 0:
			reason = f'minimap2 exited with code {proc.returncode}'
		elif 'ERROR' in stderr_text:
			reason = 'minimap2 reported an error'
		elif not os.path.exists(output_file):
			reason = 'minimap2 produced no alignment file'
		else:
			reason = None

		if reason:
			task.fail('minimap2 failed')
			log.error(f'{reason} while aligning {input_file}.\n  Command: {cmd}',
					  stderr_text or 'minimap2 produced no error message. '
					  'Check that minimap2 is installed in the active environment (which minimap2) '
					  'and that there is free disk space in the output directory.')

		if own_task:
			task.done()

# Reads that minimap2 could not place are dropped by soft_clipped() below, so
# they never reach the alignment file the read counters are computed from. That
# made "all reads" the number of *aligned* reads and left the "unmapped" column
# permanently at zero. The aligner is the only place that still sees them, so it
# records the two numbers in a small file beside the SAM.
ALIGN_STATS_SUFFIX = '.align_stats'


def align_stats_path(sam_path):
	return sam_path + ALIGN_STATS_SUFFIX


def write_align_stats(sam_path, stats):
	try:
		with open(align_stats_path(sam_path), 'w') as fw:
			for key in ('input_reads', 'unmapped'):
				fw.write(f'{key}\t{stats.get(key, 0)}\n')
	except OSError as e:
		# Counting is not worth failing a run over.
		log.warn(f'Could not write alignment statistics beside {sam_path}: {e}')


def read_align_stats(sam_path):
	"""Input read count / unmapped count recorded when sam_path was written.

	Returns None when the file is absent (an alignment produced by an older
	run, or by a code path that does not go through run_triple_minimap2), so
	callers can fall back to counting the alignment file alone.
	"""
	path = align_stats_path(sam_path)
	if not os.path.exists(path):
		return None
	stats = {}
	try:
		with open(path) as f:
			for line in f:
				key, _, value = line.strip().partition('\t')
				if key:
					stats[key] = int(value)
	except (OSError, ValueError):
		return None
	return stats or None


def soft_clipped(align_file_1, output_file, cut_len_threshold, fasta_check=False):

	check_SA = False
	read_d = []
	stats = {'input_reads': 0, 'unmapped': 0}
	with pysam.AlignmentFile(align_file_1) as f1, open(output_file, 'w') as fw:
		for line in f1:
			if line.is_secondary or line.is_supplementary:
				continue
			# minimap2 emits exactly one primary-or-unmapped record per input
			# read, so this is the input read count.
			stats['input_reads'] += 1
			if line.is_unmapped:
				stats['unmapped'] += 1
				continue
			if line.mapping_quality < 30 and align_file_1.split('/')[-1].split('_')[0] != '1':
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

			
	return check_SA, read_d, stats

def run_triple_minimap2(ref_file, input_file, output_file, longjoin_bandwidth, chaining_bandwidth, threads, minimap2_opt, len_cutoff = 100, fasta_check=False, label=None):
	output_path = output_file[:output_file.rfind('/')]
	check_SA = True
	n = 0
	read_d = {}
	align_stats = {'input_reads': 0, 'unmapped': 0}
	file_format = 'fastq'
	if fasta_check:
		file_format = 'fasta'

	if label is None:
		label = f'Aligning {os.path.basename(input_file)}'
	align_task = log.task(label)

	while check_SA == True:
		n += 1
		if n == 1:
			fastq_file = input_file
		else:
			fastq_file = output_path + f'/{n}_soft.{file_format}'

		align_task.update(f'pass {n}')
		run_minimap2(ref_file, fastq_file, output_path + f'/{n}_align.sam' , longjoin_bandwidth, chaining_bandwidth, threads, minimap2_opt, task=align_task)
		if not os.path.exists(output_path + f'/{n}_align.sam'):
			align_task.fail('no alignment file produced')
			log.error(f'minimap2 did not produce {output_path}/{n}_align.sam')
		check_SA, read_out, pass_stats = soft_clipped(output_path + f'/{n}_align.sam', output_path + f'/{n+1}_soft.{file_format}', len_cutoff, fasta_check)
		if n == 1:
			# Only the first pass sees the original reads; later passes
			# realign soft-clipped fragments of reads already counted.
			align_stats = pass_stats
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

	pass_word = 'pass' if n == 1 else 'passes'
	align_task.done(f'{n} {pass_word}, {log.count(len(read_d))} reads')

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

	write_align_stats(output_file, align_stats)
	
	# remove tmp files
	
	for fn in os.listdir(output_path):
		if fn[-10:] == '_align.sam' or fn[-11:] == '_soft.fasta':
			os.system(f'rm {output_path}/{fn}')

