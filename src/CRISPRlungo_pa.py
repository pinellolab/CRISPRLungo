# 20261111_pa
"""Prime assembly (PA) outcome classification.

Ports the read classification used for PA out-out amplicons in Levesque et al.,
Nature 2026, into CRISPRlungo. Every read is aligned against five alleles -
wild type, precise PA, dual flap without donor, 5' flap only and 3' flap only -
the best-matching allele is chosen by alignment identity, and the read is then
classified from the mutations that fall in the windows around that allele's
cleavage, flap and fragment-junction positions.

The four non-wild-type alleles are derived from the reference, the two cleavage
positions and the donor, so the only thing the user has to supply is --pa_donor.

Nothing here touches the main CRISPRlungo output: PA writes into
<output_dir>/pa_results and <output_dir>/pa_align only.
"""
import os
import shutil  #20261111_compat
import csv
from collections import OrderedDict, Counter

import pysam

import CRISPRlungo_log as log
import CRISPRlungo_mutation_analysis as mutation_analysis
from CRISPRlungo_minimap import run_triple_minimap2


REF_ORDER = ['WT', 'PA', 'DUAL', 'F5', 'F3']

REF_LABEL = {'WT': 'wild type',
			 'PA': 'precise PA',
			 'DUAL': 'dual flap, no donor',
			 'F5': "5' flap only",
			 'F3': "3' flap only"}

PATTERN_OF_REF = {'WT': 'none', 'PA': 'full_donor', 'DUAL': 'dual_flap',
				  'F5': '5p_flap', 'F3': '3p_flap'}

CATEGORIES = ['Precise_PA',
			  'Flap_install_wo_indels',
			  'Flap_integration_wDel',
			  'Flap_integration_wIns',
			  'Flap_install_with_Large_deletion',
			  'Large_deletion',
			  'Unedited',
			  'Indel_at_nick_site',
			  'Others',
			  'Not_analysed']

FLAP_PATTERNS = ['full_donor', 'dual_flap', '5p_flap', '3p_flap', 'none']


# --------------------------------------------------------------------------
# inputs

def read_fasta(path):
	records, name, seq = [], None, []
	for line in open(path):
		line = line.strip()
		if not line:
			continue
		if line.startswith('>'):
			if name is not None:
				records.append((name, ''.join(seq).upper()))
			name, seq = line[1:].split()[0], []
		else:
			seq.append(line)
	if name is not None:
		records.append((name, ''.join(seq).upper()))
	return records


def parse_flap_len(spec, original_target=1):
	"""--pa_flap_len takes '32', or '32,25' when the two targets write different flaps.

	Two values are read in the order the targets were given on the command line
	(target, additional_target). CRISPRlungo sorts the targets so that
	cv_pos < cv_pos_2 and records in original_target which one was given first,
	so the pair is swapped here when that happened.
	"""
	parts = [i.strip() for i in str(spec).split(',') if i.strip() != '']
	if len(parts) == 1:
		values = [parts[0], parts[0]]
	elif len(parts) == 2:
		values = [parts[0], parts[1]]
	else:
		log.error("--pa_flap_len takes one value, or two separated by a comma "
				  "(the flap at target, the flap at additional_target).")
	try:
		values = [int(i) for i in values]
	except ValueError:
		log.error('--pa_flap_len values must be whole numbers.')
	if min(values) <= 0:
		log.error('--pa_flap_len values must be positive.')
	if original_target == 2:
		values = values[::-1]
	return values[0], values[1]


# --------------------------------------------------------------------------
# the five alleles

def build_alleles(ref_seq, cv_pos, cv_pos_2, donor_records, flap_5, flap_3, pa_override=None):
	"""Return {name: (sequence, sites)}.

	Positions follow CRISPRlungo: 0-based, and the cut lies after that index.
	`sites` are the positions whose windows decide the classification on that
	allele: the two genome junctions, the flap/donor-body boundaries, and the
	joins between donor fragments.
	"""
	donor = ''.join(seq for _, seq in donor_records)
	left, right = ref_seq[:cv_pos + 1], ref_seq[cv_pos_2 + 1:]
	f5seq, f3seq = donor[:flap_5], donor[-flap_3:]

	alleles = OrderedDict()
	alleles['WT'] = (ref_seq, {'cleavage': [cv_pos, cv_pos_2], 'flap': [], 'junction': []})

	if pa_override is not None:
		pa_seq = pa_override
		start = pa_seq.find(donor)
		if start < 0:
			log.error('--pa_reference does not contain the --pa_donor sequence.')
		base = start - 1
	else:
		pa_seq = left + donor + right
		base = cv_pos
	junctions, offset = [], 0
	for _, seq in donor_records[:-1]:
		offset += len(seq)
		junctions.append(base + offset)
	alleles['PA'] = (pa_seq, {'cleavage': [base, base + len(donor)],
							  'flap': [base + flap_5, base + len(donor) - flap_3],
							  'junction': junctions})

	alleles['DUAL'] = (left + f5seq + f3seq + right,
					   {'cleavage': [cv_pos, cv_pos + flap_5 + flap_3],
						'flap': [cv_pos + flap_5], 'junction': []})
	alleles['F5'] = (left + f5seq + right,
					 {'cleavage': [cv_pos, cv_pos + flap_5], 'flap': [], 'junction': []})
	alleles['F3'] = (left + f3seq + right,
					 {'cleavage': [cv_pos, cv_pos + flap_3], 'flap': [], 'junction': []})
	return alleles, donor, (f5seq, f3seq)


def validate(ref_seq, donor, flap_5, flap_3):
	if flap_5 + flap_3 >= len(donor):
		log.error('--pa_flap_len: the two flaps together are longer than the donor.')
	if donor[:flap_5] in ref_seq:
		log.warn("PA: the 5' flap sequence also occurs in the reference; "
				 "calls at that junction may be ambiguous.")
	if donor[-flap_3:] in ref_seq:
		log.warn("PA: the 3' flap sequence also occurs in the reference; "
				 "calls at that junction may be ambiguous.")


# --------------------------------------------------------------------------
# reading a run back

def parse_mutation_info(text):
	"""results/read_classification.txt column 3 -> [(start, end, kind, length)]."""
	out = []
	if not text or text == 'None':
		return out
	for token in text.split(','):
		token = token.strip()
		if not token or ':' not in token:
			continue
		span, detail = token.split(':', 1)
		if '_' not in span:
			continue
		a, b = span.split('_', 1)
		try:
			start, end = int(a), int(b)
		except ValueError:
			continue
		pieces = detail.split('_')
		kind = pieces[0]
		try:
			length = int(pieces[1])
		except (IndexError, ValueError):
			length = end - start + 1
		out.append((start, end, kind, length))
	return out


def base_read_id(read_id):
	for marker in ('_alignStrandInfo_', '_AI_'):
		if marker in read_id:
			read_id = read_id.split(marker)[0]
			break
	if '_split_' in read_id:
		read_id = read_id.split('_split_')[0]
	return read_id


def load_align_fraction(path):
	"""Per read, the fraction of its bases that this allele accounts for.

	A read is assigned to the allele that explains the most of it. Without this an
	allele far shorter than the read wins by default, because everything it cannot
	explain is simply clipped away and costs nothing.
	"""
	aligned, length = {}, {}
	handle = pysam.AlignmentFile(path, 'r', check_sq=False)
	for read in handle:
		if read.is_unmapped or read.is_secondary:
			continue
		rid = base_read_id(read.query_name)
		aligned[rid] = aligned.get(rid, 0) + (read.query_alignment_length or 0)
		length[rid] = max(length.get(rid, 0), read.infer_read_length() or 0)
	handle.close()
	return {rid: (aligned[rid] / length[rid]) for rid in aligned if length.get(rid)}


def load_run(path):
	table = OrderedDict()
	with open(path) as handle:
		reader = csv.reader(handle, delimiter='\t')
		next(reader, None)
		for row in reader:
			if len(row) < 3:
				continue
			rid = base_read_id(row[0])
			table.setdefault(rid, [])
			table[rid].extend(parse_mutation_info(row[2]))
	return table


# --------------------------------------------------------------------------
# scoring

def overlaps(mut, site, half):
	"""A site is the index the cut sits after, so its window covers site and site+1."""
	return not (mut[1] < site - half or mut[0] > site + 1 + half)


def window_identity(muts, cleavage):
	"""The published identity: matches / (matches + substitutions + insertions).

	Deletions are deliberately left out of the denominator, which is what lets a
	read that lost the donor body still match the allele it came from.
	"""
	lo, hi = min(cleavage), max(cleavage) + 1
	span = max(1, hi - lo + 1)
	sub = ins = dele = 0
	for start, end, kind, length in muts:
		if end < lo or start > hi:
			continue
		if kind == 'Sub':
			sub += length
		elif kind in ('Ins', 'LargeIns'):
			ins += length
		elif kind in ('Del', 'LargeDel'):
			dele += min(length, span)
	matches = max(0, span - sub - dele)
	denom = matches + sub + ins
	return (matches / denom) if denom else 0.0


def classify_read(per_ref, alleles, half, large_cut, covered=None):
	"""per_ref = {allele name: mutation list}, for the alleles that kept this read."""
	if not per_ref:
		return dict(category='Not_analysed', used_ref='.', flap_pattern='.',
					identity={}, del_len=0, ins_len=0)

	identity = {name: window_identity(muts, alleles[name][1]['cleavage'])
				for name, muts in per_ref.items()}

	# Which allele the read came from. The published identity
	# (matches / (matches + sub + ins)) leaves deletions out, so with five alleles it
	# returns 1.0 for several of them at once and cannot separate them. Instead: the
	# allele that explains the most of the read wins, and among alleles that explain it
	# equally well, the one needing the least editing to do so. Identity is still
	# reported per allele in pa_read_classification.txt and breaks the remaining ties.
	covered = covered or {}
	burden = {name: sum(m[3] for m in muts) for name, muts in per_ref.items()}
	def rank(name):
		return (-round(covered.get(name, 0.0), 2), burden[name], -identity[name])
	best = None
	for name in REF_ORDER:
		if name not in per_ref:
			continue
		if best is None or rank(name) < rank(best):
			best = name

	muts = per_ref[best]
	sites = (alleles[best][1]['cleavage'] + alleles[best][1]['flap']
			 + alleles[best][1]['junction'])

	in_window = [m for m in muts if any(overlaps(m, s, half) for s in sites)]
	del_len = sum(m[3] for m in in_window if m[2] in ('Del', 'LargeDel'))
	ins_len = sum(m[3] for m in in_window if m[2] in ('Ins', 'LargeIns'))
	large_del = any(m[2] == 'LargeDel' and m[3] > large_cut for m in in_window)

	if not in_window:
		category = 'Unedited' if best == 'WT' else (
			'Precise_PA' if best == 'PA' else 'Flap_install_wo_indels')
	elif best == 'WT':
		category = 'Large_deletion' if large_del else (
			'Indel_at_nick_site' if (del_len or ins_len) else 'Others')
	elif large_del:
		category = 'Flap_install_with_Large_deletion'
	elif del_len > ins_len:
		category = 'Flap_integration_wDel'
	elif ins_len > del_len:
		category = 'Flap_integration_wIns'
	else:
		category = 'Others'

	return dict(category=category, used_ref=best, flap_pattern=PATTERN_OF_REF[best],
				identity=identity, del_len=del_len, ins_len=ins_len)


def near_edge(muts, sites, half, slack=2):
	for mut in muts:
		for site in sites:
			for edge in (site - half, site + 1 + half):
				if abs(mut[0] - edge) <= slack or abs(mut[1] - edge) <= slack:
					return True
	return False


# --------------------------------------------------------------------------
# one analysis pass per allele

def analyse_allele(name, seq, sites, ref_name, pa_dir, align_dir, treated_file, control_file,
				   window, threads, args, range_align_end, longjoin_bandwidth,
				   chaining_bandwidth, treated_sam=None, control_sam=None, fasta_path=None):
	# Laid out like a normal CRISPRlungo run: <root>/ref_seq beside <root>/results.
	# analysis_function_without_control() hands the output directory to
	# make_visualization_sam(), which looks for ../ref_seq/ref_wo_umi.fasta, so a
	# flat directory per allele crashes the no-control path.
	run_root = os.path.join(pa_dir, 'reference_' + name.lower())  #20261111_pa
	out_dir = os.path.join(run_root, 'results')  #20261111_pa
	ref_dir = os.path.join(run_root, 'ref_seq')  #20261111_pa
	os.makedirs(out_dir, exist_ok=True)  #20261111_pa
	os.makedirs(ref_dir, exist_ok=True)  #20261111_pa

	allele_fasta = os.path.join(ref_dir, 'ref_wo_umi.fasta')  #20261111_pa
	with open(allele_fasta, 'w') as handle:  #20261111_pa
		handle.write(f'>{ref_name}_{name}\n{seq}\n')  #20261111_pa
	if fasta_path is None:
		fasta_path = allele_fasta

	# named as CRISPRlungoAllele expects (<root>/ref_seq/ref_wo_umi.mmi), and built
	# for every allele so that Allele can be pointed at any per-allele directory
	index = os.path.join(ref_dir, 'ref_wo_umi.mmi')  #20261111_compat
	os.system(f'minimap2 -d {index} {allele_fasta} > /dev/null 2>&1')  #20261111_compat

	if treated_sam is None:
		treated_sam = os.path.join(align_dir, f'{name}_Treated_alignment.sam')
		run_triple_minimap2(index, treated_file, treated_sam, longjoin_bandwidth,
							chaining_bandwidth, threads, args.minimap2_opt,
							len_cutoff=args.align_sa_len_threshold, fasta_check=True)
		if control_file:
			control_sam = os.path.join(align_dir, f'{name}_Control_alignment.sam')
			run_triple_minimap2(index, control_file, control_sam, longjoin_bandwidth,
								chaining_bandwidth, threads, args.minimap2_opt,
								len_cutoff=args.align_sa_len_threshold, fasta_check=True)

	# CRISPRlungoAllele's last step reads <root>/align/Treated_alignment.sam, so
	# give each allele directory that path too. Hard-linked where the filesystem
	# allows it, since these files are large.
	link_dir = os.path.join(run_root, 'align')  #20261111_compat
	os.makedirs(link_dir, exist_ok=True)  #20261111_compat
	for src, dst_name in ((treated_sam, 'Treated_alignment.sam'),  #20261111_compat
						  (control_sam, 'Control_alignment.sam')):  #20261111_compat
		if not src or not os.path.exists(src):  #20261111_compat
			continue  #20261111_compat
		dst = os.path.join(link_dir, dst_name)  #20261111_compat
		if os.path.exists(dst):  #20261111_compat
			os.remove(dst)  #20261111_compat
		try:  #20261111_compat
			os.link(src, dst)  #20261111_compat
		except OSError:  #20261111_compat
			shutil.copyfile(src, dst)  #20261111_compat

	cv_pos, cv_pos_2 = sites['cleavage']

	# The keys CRISPRlungoAllele reads back, written with the same names and in
	# the same order as CRISPRlungo.main(), so `CRISPRlungoAllele <this dir>`
	# works on an allele exactly as it does on a normal run. The coordinates are
	# the allele's own, not the WT reference's.
	with open(os.path.join(out_dir, 'input_summary.txt'), 'w') as handle:  #20261111_compat
		handle.write(f'Target_1 :{args.target}\n')  #20261111_compat
		handle.write(f'Target_2 :{args.additional_target}\n')  #20261111_compat
		handle.write(f'Ref_seq :{seq}\n')  #20261111_compat
		handle.write(f'CleavagePos_1 :{cv_pos}\n')  #20261111_compat
		handle.write('CleavageStrand_1 :1\n')  #20261111_compat
		handle.write(f'CleavagePos_2 :{cv_pos_2}\n')  #20261111_compat
		handle.write('CleavageStrand_2 :-1\n')  #20261111_compat
		handle.write(f'Window : {window}\n')  #20261111_compat
		handle.write('Window_between_cleavage : False\n')  #20261111_compat
		handle.write('induced_mut :False\n')  #20261111_compat
		handle.write(f'cut_pos_in_target :{args.cleavage_pos}\n')  #20261111_compat
		handle.write('original_target :1\n')  #20261111_compat
		handle.write(f'minimum_read_count : {args.min_read_cnt}\n')  #20261111_compat
		handle.write(f'minimum_read_frequency : {args.min_read_freq}\n')  #20261111_compat
		handle.write('induced_sequence_path : None\n')  #20261111_compat
		handle.write(f'largeins_cutlen : {args.largeins_cutlen}\n')  #20261111_compat
		handle.write(f'largedel_cutlen : {args.largedel_cutlen}\n')  #20261111_compat
		handle.write('\n# Run\n')  #20261111_compat
		handle.write(f'PA_allele : {name}\n')  #20261111_compat
		handle.write('PA_note : written by CRISPRlungo_pa; coordinates are this allele\'s own\n')  #20261111_compat

	table_path = os.path.join(out_dir, 'read_classification.txt')
	if control_sam:
		edited_dict, _ctrl, ids, control_cnt, edited_cnt = mutation_analysis.analysis_function_with_control(
			control_sam, treated_sam, fasta_path, out_dir + '/',
			cv_pos, cv_pos_2, window, True, [], [False],
			args.mutation_length_threshold_pval,
			range_align_end=range_align_end, threads=threads,
			largeins_cutlen=args.largeins_cutlen, largedel_cutlen=args.largedel_cutlen,
			use_all_mutations=args.using_all_mutations, length_min=args.length_min,
			allowance_value=0, umi_clustered=False,
			background_max_len=args.largedel_cutlen)
		mutation_analysis.process_mutations(edited_dict, table_path, ids, [], [False],
											args.induced_paritial_similiarity,
											args.largeins_cutlen, args.largedel_cutlen, seq)
		mutation_analysis.write_cnt_file(control_cnt, edited_cnt,
										 os.path.join(out_dir, 'preprocess_count.txt'))
	else:
		mutation_analysis.analysis_function_without_control(
			seq, ref_name + '_' + name, cv_pos, 1, cv_pos_2, -1, window, treated_sam,
			out_dir, True, [], os.path.dirname(os.path.abspath(__file__)), [False], False,  #20261111_pa
			threads=threads, mix_tag=args.mix_tag, write_cnt=True,
			partial_induce_cutoff=args.induced_paritial_similiarity,
			range_align_end=range_align_end,
			largeins_cutlen=args.largeins_cutlen, largedel_cutlen=args.largedel_cutlen)
	return table_path


# --------------------------------------------------------------------------

def run(args, output_dir, ref_seq, ref_name, cv_pos, cv_pos_2, original_target,
		treated_file, control_file, threads, range_align_end,
		longjoin_bandwidth, chaining_bandwidth):

	if cv_pos_2 is False or cv_pos_2 is None:
		log.error('--pa_donor needs two targets; add --additional_target.')

	pa_dir = os.path.join(output_dir, 'pa_results')
	align_dir = os.path.join(output_dir, 'pa_align')
	os.makedirs(pa_dir, exist_ok=True)
	os.makedirs(align_dir, exist_ok=True)

	log.blank()
	log.info('Prime assembly (PA) analysis')
	if not control_file:  #20261111_pa
		log.warn('PA: no --control. Every sequencing error inside a window counts as a '  #20261111_pa
				 'mutation, so Precise_PA is under-reported and Indel_at_nick_site inflated.')  #20261111_pa

	flap_5, flap_3 = parse_flap_len(args.pa_flap_len, original_target)
	donor_records = read_fasta(args.pa_donor)
	if not donor_records:
		log.error(f'No sequence found in {args.pa_donor}.')

	pa_override = None
	if args.pa_reference:
		records = read_fasta(args.pa_reference)
		if len(records) != 1:
			log.error('--pa_reference must contain exactly one sequence.')
		pa_override = records[0][1]

	alleles, donor, _flap_seqs = build_alleles(ref_seq, cv_pos, cv_pos_2, donor_records,
											   flap_5, flap_3, pa_override)
	validate(ref_seq, donor, flap_5, flap_3)

	log.info(f'  Donor           : {len(donor)} bp in {len(donor_records)} fragment(s)')
	log.info(f"  Flap lengths    : 5' {flap_5} nt, 3' {flap_3} nt")
	for name in REF_ORDER:
		seq, sites = alleles[name]
		marks = sites['cleavage'] + sites['flap'] + sites['junction']
		log.info(f'  {name:<4} {REF_LABEL[name]:<22} {len(seq):>6} bp   sites {marks}')

	with open(os.path.join(pa_dir, 'pa_coordinates.txt'), 'w') as handle:
		handle.write(f'donor_file\t{args.pa_donor}\n')
		handle.write(f'donor_fragments\t{len(donor_records)}\n')
		handle.write(f'donor_length\t{len(donor)}\n')
		handle.write(f'flap_5\t{flap_5}\n')
		handle.write(f'flap_3\t{flap_3}\n')
		handle.write(f'pa_window\t{args.pa_window}\n')
		handle.write(f'pa_site_window\t{args.pa_site_window}\n')
		handle.write('\nreference\tlength\tcleavage\tflap\tjunction\n')
		for name in REF_ORDER:
			seq, sites = alleles[name]
			handle.write('%s\t%d\t%s\t%s\t%s\n' % (
				name, len(seq),
				','.join(str(i) for i in sites['cleavage']),
				','.join(str(i) for i in sites['flap']) or '-',
				','.join(str(i) for i in sites['junction']) or '-'))

	for name in REF_ORDER:
		seq, _sites = alleles[name]
		with open(os.path.join(pa_dir, f'pa_reference_{name}.fasta'), 'w') as handle:
			handle.write(f'>{ref_name}_{name}\n{seq}\n')

	tables = {}
	coverage = {}
	for name in REF_ORDER:
		seq, sites = alleles[name]
		if name == 'WT':
			table = analyse_allele(
				name, seq, sites, ref_name, pa_dir, align_dir, treated_file, control_file,
				args.pa_window, threads, args, range_align_end, longjoin_bandwidth,
				chaining_bandwidth,
				treated_sam=os.path.join(output_dir, 'align', 'Treated_alignment.sam'),
				control_sam=(os.path.join(output_dir, 'align', 'Control_alignment.sam')
							 if control_file else None),
				fasta_path=os.path.join(output_dir, 'ref_seq', 'ref_wo_umi.fasta'))
		else:
			table = analyse_allele(
				name, seq, sites, ref_name, pa_dir, align_dir, treated_file, control_file,
				args.pa_window, threads, args, range_align_end, longjoin_bandwidth,
				chaining_bandwidth)
		tables[name] = load_run(table)
		sam = (os.path.join(output_dir, 'align', 'Treated_alignment.sam') if name == 'WT'
			   else os.path.join(align_dir, f'{name}_Treated_alignment.sam'))
		coverage[name] = load_align_fraction(sam) if os.path.exists(sam) else {}

	half = args.pa_site_window
	large_cut = args.largedel_cutlen
	counts, patterns, chosen = Counter(), Counter(), Counter()
	edge_cases = 0

	order, seen = [], set()
	for name in REF_ORDER:
		for rid in tables[name]:
			if rid not in seen:
				seen.add(rid)
				order.append(rid)

	rows_path = os.path.join(pa_dir, 'pa_read_classification.txt')
	with open(rows_path, 'w', newline='') as handle:
		writer = csv.writer(handle, delimiter='\t')
		writer.writerow(['Read_id', 'Category', 'Best_reference', 'Flap_pattern',
						 'Deleted_bp', 'Inserted_bp']
						+ [f'Identity_{n}' for n in REF_ORDER] + ['References_kept'])
		for rid in order:
			per_ref = {n: tables[n][rid] for n in REF_ORDER if rid in tables[n]}
			covered = {n: coverage[n][rid] for n in per_ref if rid in coverage.get(n, {})}
			result = classify_read(per_ref, alleles, half, large_cut, covered)
			counts[result['category']] += 1
			patterns[result['flap_pattern']] += 1
			chosen[result['used_ref']] += 1
			best = result['used_ref']
			if best in alleles:
				sites = (alleles[best][1]['cleavage'] + alleles[best][1]['flap']
						 + alleles[best][1]['junction'])
				if near_edge(per_ref[best], sites, half):
					edge_cases += 1
			writer.writerow([rid, result['category'], best, result['flap_pattern'],
							 result['del_len'], result['ins_len']]
							+ [f"{result['identity'][n]:.4f}" if n in result['identity'] else '.'
							   for n in REF_ORDER]
							+ [''.join(n[0] if n != 'F3' else '3' for n in REF_ORDER if n in per_ref)])

	total = sum(counts.values())
	analysed = total - counts['Not_analysed']
	with open(os.path.join(pa_dir, 'pa_outcome_count.txt'), 'w') as handle:
		handle.write('Category\tCount\tPercent_of_analysed\n')
		for name in CATEGORIES:
			pct = (100.0 * counts[name] / analysed) if analysed else 0.0
			handle.write(f'{name}\t{counts[name]}\t{pct:.2f}\n')
		handle.write(f'#analysed_reads\t{analysed}\t\n')
		handle.write(f'#reads_with_a_breakpoint_within_2bp_of_a_window_edge\t{edge_cases}\t\n')
		handle.write('\nFlap_pattern\tCount\tPercent_of_analysed\n')
		for name in FLAP_PATTERNS:
			pct = (100.0 * patterns[name] / analysed) if analysed else 0.0
			handle.write(f'{name}\t{patterns[name]}\t{pct:.2f}\n')
		handle.write('\nBest_reference\tCount\n')
		for name in REF_ORDER:
			handle.write(f'{name}\t{chosen[name]}\n')

	with open(os.path.join(pa_dir, 'pa_summary.txt'), 'w') as handle:
		handle.write('\t'.join(['Analysed'] + CATEGORIES) + '\n')
		handle.write('\t'.join([str(analysed)] + [str(counts[n]) for n in CATEGORIES]) + '\n')

	log.summary_table('PA outcome', OrderedDict((n, counts[n]) for n in CATEGORIES if counts[n]))
	if edge_cases:
		log.warn(f'PA: {edge_cases} reads carry a breakpoint within 2 bp of a classification '
				 f'window edge; --pa_site_window may be too tight for this data.')
	log.info(f'PA results : {pa_dir}')

	# 20261111_pa_report: what the report needs to show PA results in place of the
	# wild-type ones. Nothing here is written to disk; the caller decides.
	labels = pa_labels_by_mutstr(output_dir, pa_dir)
	return {
		'counts': counts, 'patterns': patterns, 'chosen': chosen,
		'analysed': analysed, 'alleles': alleles, 'flap': (flap_5, flap_3),
		'labels': labels, 'dir': pa_dir,
	}


# ---------------------------------------------------------------- reporting

# The report shows the unedited reads first, under the same name the standard
# outcome table uses for them.
DISPLAY_ORDER = ['Unedited', 'Precise_PA', 'Flap_install_wo_indels',
	'Flap_integration_wDel', 'Flap_integration_wIns',
	'Flap_install_with_Large_deletion', 'Large_deletion',
	'Indel_at_nick_site', 'Others']
DISPLAY_NAME = {'Unedited': 'WT'}


def display_name(category):
	return DISPLAY_NAME.get(category, category.replace('_', ' '))


def pa_labels_by_mutstr(output_dir, pa_dir):
	"""Mutation pattern -> the PA category most of its reads fall in.

	An allele-plot row is one mutation pattern, not one read, and the reads behind
	it need not share a PA category, so the row is labelled with the category that
	holds most of them."""
	cls = os.path.join(output_dir, 'results', 'read_classification.txt')
	pa = os.path.join(pa_dir, 'pa_read_classification.txt')
	if not (os.path.exists(cls) and os.path.exists(pa)):
		return {}

	category = {}
	with open(pa) as handle:
		next(handle, None)
		for line in handle:
			f = line.rstrip('\n').split('\t')
			if len(f) > 1:
				category[base_read_id(f[0])] = f[1]

	tally = {}
	with open(cls) as handle:
		next(handle, None)
		for line in handle:
			f = line.rstrip('\n').split('\t')
			if len(f) < 3:
				continue
			cat = category.get(base_read_id(f[0]))
			if cat is None:
				continue
			bucket = tally.setdefault(f[2], {})
			bucket[cat] = bucket.get(cat, 0) + 1

	labels = {}
	for mutstr, cats in tally.items():
		best = max(cats, key=lambda c: cats[c])
		labels[mutstr] = display_name(best)
	return labels


def pa_table_html(pa, edited_reads_cnt):
	"""The heading, the explanation and the outcome table, in place of the
	standard one. The heading sits above the table it describes."""
	counts, analysed = pa['counts'], pa['analysed']
	f5, f3 = pa['flap']
	allele = pa['alleles']['PA']
	cleavage = ' and '.join(str(x) for x in allele[1]['cleavage'])
	flaps = ' and '.join(str(x) for x in allele[1]['flap'])

	def pct(v):
		return f'{v} ({(100.0 * v / analysed) if analysed else 0.0:.1f} %)'

	note = ('Every read is compared with the wild-type, precise-PA, dual-flap, '
			'5&prime;-flap and 3&prime;-flap alleles and assigned to the one that explains '
			'the most of it. The table, the pie charts and the allele-plot labels below are '
			'that classification; outcome names follow Levesque et al., Nature 2026. '
			f'Flaps 5&prime; {f5} nt and 3&prime; {f3} nt; on the precise-PA allele the '
			f'cleavage sites are at {cleavage} and the flap boundaries at {flaps}. '
			f'Percentages are of the {analysed:,} reads at least one of the five alleles '
			'could place &mdash; another allele places reads the wild-type analysis drops, so '
			'this is not the same total as a run without PA. The per-position insertion, '
			'deletion and substitution plots further down stay on the wild-type reference, '
			'since PA has no equivalent of them. The standard files in results/ are unchanged; '
			'the PA tables are in pa_results/.')

	# One table per block, the way the web version lays it out: read counts, then
	# the PA outcomes, then what the standard classification made of the reads PA
	# could not place.
	def table(headers, values):
		width = 100.0 / len(headers)
		cols = ''.join(f'<col style="width: {width}%;">' for _ in headers)
		head = ''.join(f'<th>{h}</th>' for h in headers)
		body = ''.join(f'<td>{v}</td>' for v in values)
		return ('<table class="table table-bordered table-striped w-100 result-table mb-3">'
				f'<colgroup>{cols}</colgroup><thead class="thead-light"><tr>{head}</tr></thead>'
				f'<tbody><tr>{body}</tr></tbody></table>')

	html = table(['Total reads', 'Unmapped', 'Reads used in analysis'],
				 [edited_reads_cnt.get('all_reads', '-'),
				  edited_reads_cnt.get('unmapped', '-'), analysed])
	html += table([display_name(c) for c in DISPLAY_ORDER],
				  [pct(counts.get(c, 0)) for c in DISPLAY_ORDER])

	return (
		'<div id="paHeading" class="mb-2"><h2 class="h5 m-0">Prime assembly analysis</h2>'
		f'<p class="small text-secondary mt-2 mb-2">{note}</p></div>'
		f'<div id="cnt-table-container" class="mb-3">{html}</div>')


def pa_pie_html(pa):
	"""The two pie charts, over the PA categories: every outcome, and unedited
	against edited."""
	import plotly.graph_objects as go
	counts, analysed = pa['counts'], pa['analysed']
	shown = [c for c in DISPLAY_ORDER if counts.get(c, 0) > 0]
	colours = ['#CECFD1', '#4357A7', '#6BBD45', '#FABF37', '#98D9ED',
			   '#BF71AE', '#ED2024', '#8c6d4a', '#000000']
	fig1 = go.Figure(data=[go.Pie(
		labels=[display_name(c) for c in shown],
		values=[counts[c] for c in shown],
		marker=dict(colors=['#CECFD1' if c == 'Unedited' else colours[i % len(colours)]
							for i, c in enumerate(shown)]))])
	fig1.update_layout(title='Prime assembly outcomes')

	unedited = counts.get('Unedited', 0)
	fig2 = go.Figure(data=[go.Pie(labels=['Unedited', 'Edited'],
								  values=[unedited, max(0, analysed - unedited)],
								  marker=dict(colors=['#CECFD1', '#4357A7']))])
	fig2.update_layout(title='Edited reads (PA)')
	return ([fig1.to_html(full_html=False, include_plotlyjs='cdn')],
			[fig2.to_html(full_html=False, include_plotlyjs='cdn')])


def cli():
	"""CRISPRlungoPA: run the PA stage against a finished CRISPRlungo run directory."""
	import argparse
	parser = argparse.ArgumentParser(prog='CRISPRlungoPA',
		description='Prime assembly outcome classification for a finished CRISPRlungo run')
	parser.add_argument('run_dir', help='A CRISPRlungo output directory')
	parser.add_argument('--pa_donor', required=True)
	parser.add_argument('--pa_flap_len', type=str, default='32')
	parser.add_argument('--pa_window', type=int, default=2)
	parser.add_argument('--pa_site_window', type=int, default=1)
	parser.add_argument('--pa_reference', type=str, default=None)
	parser.add_argument('--largeins_cutlen', type=int, default=50)
	parser.add_argument('--largedel_cutlen', type=int, default=50)
	parser.add_argument('--length_min', type=float, default=10)
	parser.add_argument('--mutation_length_threshold_pval', type=float, default=-1)
	parser.add_argument('--induced_paritial_similiarity', type=float, default=0.8)
	parser.add_argument('--align_sa_len_threshold', type=int, default=100)
	parser.add_argument('--range_both_end_region', type=int, default=100)
	parser.add_argument('--minimap2_opt', type=str, default='')
	parser.add_argument('--using_all_mutations', action='store_true')
	parser.add_argument('--mix_tag', action='store_true')
	parser.add_argument('-t', '--threads', type=int, default=1)
	args = parser.parse_args()

	summary = {}
	for line in open(os.path.join(args.run_dir, 'results', 'input_summary.txt')):
		if ':' not in line:
			continue
		key, value = line.split(':', 1)
		summary.setdefault(key.strip(), value.strip())

	ref_seq = summary['Ref_seq'].upper()
	cv_pos = int(summary['CleavagePos_1'])
	cv_pos_2 = summary['CleavagePos_2']
	cv_pos_2 = False if cv_pos_2.lower() == 'false' else int(cv_pos_2)
	original_target = int(summary.get('original_target', 1))
	ref_name = summary.get('reference_name', 'Ref')
	treated = summary.get('treated_file_used')
	control = summary.get('control_file_used')
	control = None if control in (None, 'None', '') else control
	longjoin = int(len(ref_seq) * 0.3)
	chaining = longjoin if longjoin > 500 else 500

	run(args, args.run_dir, ref_seq, ref_name, cv_pos, cv_pos_2, original_target,
		treated, control, args.threads, args.range_both_end_region, longjoin, chaining)
