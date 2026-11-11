"""Simulate a PA out-out amplicon library with a known allele mixture.

Every read name carries its true class, so PA mode can be scored read by read.
"""
import random, sys

def fa(path):
	return open(path).read().split('\n', 1)[1].replace('\n', '')

AMP = fa('/home/claude/work/ref/amp_WT.fasta')      # 426 bp out-out amplicon
DONOR = fa('/home/claude/work/ref/donor.fasta')     # 2483 bp donor
CV, CV2, FLAP = 80, 265, 32                         # 0-based, cut after index

LEFT, RIGHT = AMP[:CV + 1], AMP[CV2 + 1:]
F5, F3 = DONOR[:FLAP], DONOR[-FLAP:]

ALLELES = {
	'Unedited':   AMP,
	'Precise_PA': LEFT + DONOR + RIGHT,
	'Dual_flap':  LEFT + F5 + F3 + RIGHT,
	'Flap_5p':    LEFT + F5 + RIGHT,
	'Flap_3p':    LEFT + F3 + RIGHT,
	'Nick_indel': AMP[:CV - 2] + AMP[CV + 6:],            # 8 bp deletion over the left nick
	'Nick_ins':   AMP[:CV + 1] + 'GATCAGTC' + AMP[CV + 1:],  # 8 bp insertion at the left nick
}

MIX = {'Unedited': 0.55, 'Precise_PA': 0.10, 'Dual_flap': 0.15,
	   'Flap_5p': 0.06, 'Flap_3p': 0.06, 'Nick_indel': 0.04, 'Nick_ins': 0.04}

SUB, INS, DEL = 0.015, 0.010, 0.015      # ~96 % identity, Plasmidsaurus-like
BASES = 'ACGT'

def corrupt(seq, rng):
	out = []
	for base in seq:
		r = rng.random()
		if r < DEL:
			continue
		if r < DEL + INS:
			out.append(rng.choice(BASES))
			out.append(base)
			continue
		if r < DEL + INS + SUB:
			out.append(rng.choice([b for b in BASES if b != base]))
			continue
		out.append(base)
	return ''.join(out)

def write(path, counts, seed, tag):
	rng = random.Random(seed)
	n = 0
	with open(path, 'w') as handle:
		for name, want in counts.items():
			seq = ALLELES[name]
			for i in range(want):
				read = corrupt(seq, rng)
				if rng.random() < 0.5:
					read = read[::-1].translate(str.maketrans('ACGT', 'TGCA'))
				n += 1
				handle.write(f'@{tag}{n}_TRUTH_{name}\n{read}\n+\n{"I" * len(read)}\n')
	return n

if __name__ == '__main__':
	total = int(sys.argv[1]) if len(sys.argv) > 1 else 4000
	counts = {k: int(round(total * v)) for k, v in MIX.items()}
	n = write('pa_sim.fastq', counts, 20261111, 'r')
	c = write('pa_sim_control.fastq', {'Unedited': total}, 909090, 'c')
	print('treated reads:', n, counts)
	print('control reads:', c)
	for k, v in ALLELES.items():
		print(f'  {k:12s} {len(v):5d} bp')
