"""Simulate a 50 bp deletion replaced by an 80 bp insert carrying ~35 substitutions.

The question is whether CRISPRlungo still reports that edit as one allele rather
than shredding it into a pile of separate mutations. Every read name carries its
true class, so the answer can be counted rather than eyeballed.
"""
import random, sys

def fa(path):
	return open(path).read().split('\n', 1)[1].replace('\n', '')

AMP = fa('/home/claude/work/ref/amp_WT.fasta')   # 426 bp
CV = 80                                          # 0-based, cut after this index
DEL_LEN = 50
INS_LEN = 80
N_SUB = 35

rng = random.Random(20261111)
BASES = 'ACGT'

# The replacement: 80 bp built from the 80 bp that would sit at the cut, with 35
# of its bases changed, so it is recognisably related to the locus but heavily
# diverged - the hard case for an aligner.
template = AMP[CV + 1: CV + 1 + INS_LEN]
insert = list(template)
for pos in rng.sample(range(INS_LEN), N_SUB):
	insert[pos] = rng.choice([b for b in BASES if b != insert[pos]])
INSERT = ''.join(insert)

LEFT = AMP[:CV + 1]
RIGHT = AMP[CV + 1 + DEL_LEN:]
EDITED = LEFT + INSERT + RIGHT          # 50 bp out, 80 bp in

ALLELES = {
	'Unedited': AMP,
	'Replaced': EDITED,
	'Del_only': LEFT + RIGHT,                              # the 50 bp deletion alone
	'Indel_8':  AMP[:CV - 2] + AMP[CV + 6:],               # a plain 8 bp deletion
}
MIX = {'Unedited': 2400, 'Replaced': 1200, 'Del_only': 250, 'Indel_8': 150}

SUB, INS, DEL = 0.015, 0.010, 0.015     # ~96 % read identity, Plasmidsaurus-like

def corrupt(seq, rng):
	out = []
	for base in seq:
		r = rng.random()
		if r < DEL:
			continue
		if r < DEL + INS:
			out.append(rng.choice(BASES)); out.append(base); continue
		if r < DEL + INS + SUB:
			out.append(rng.choice([b for b in BASES if b != base])); continue
		out.append(base)
	return ''.join(out)

def write(path, counts, seed, tag):
	r = random.Random(seed)
	n = 0
	with open(path, 'w') as handle:
		for name, want in counts.items():
			for _ in range(want):
				read = corrupt(ALLELES[name], r)
				if r.random() < 0.5:
					read = read[::-1].translate(str.maketrans('ACGT', 'TGCA'))
				n += 1
				handle.write(f'@{tag}{n}_TRUTH_{name}\n{read}\n+\n{"I" * len(read)}\n')
	return n

if __name__ == '__main__':
	with open('induced.fasta', 'w') as h:
		h.write('>replacement_allele\n' + EDITED + '\n')
	with open('insert_only.fasta', 'w') as h:
		h.write('>insert80\n' + INSERT + '\n')
	n1 = write('treated.fastq', MIX, 1, 'T')
	n2 = write('control.fastq', {'Unedited': 3000}, 2, 'C')
	ident = sum(1 for a, b in zip(INSERT, template) if a == b) / INS_LEN
	print(f'amplicon {len(AMP)} bp, cut after {CV}')
	print(f'deleted {DEL_LEN} bp, inserted {INS_LEN} bp with {N_SUB} substitutions '
		  f'({ident*100:.0f} % identity to the sequence it replaces)')
	print(f'edited allele {len(EDITED)} bp (net +{len(EDITED)-len(AMP)} bp)')
	print(f'treated {n1} reads, control {n2} reads')
