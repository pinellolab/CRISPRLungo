#!/usr/bin/env bash
# 20261111 — 50 bp deleted and replaced by an 80 bp insert carrying 35 substitutions.
# Does CRISPRlungo still report that edit? Simulated reads, truth in every read name.
set -e

python simulate_replacement.py        # writes treated.fastq, control.fastq, induced.fasta

# 1. defaults: detection is fine, "Precise" is not
CRISPRlungo amp_WT.fasta treated.fastq out_default TATTGAGAAGCCGAAACCTC \
    --control control.fastq --induced_sequence_path induced.fasta --window 60

# 2. the setting that matters: the replacement is 19 induced patterns, and an ONT
#    read rarely carries all of them, so the match threshold has to come down
CRISPRlungo amp_WT.fasta treated.fastq out_sim05 TATTGAGAAGCCGAAACCTC \
    --control control.fastq --induced_sequence_path induced.fasta --window 60 \
    --induced_paritial_similiarity 0.5

python score2.py out_default
python score2.py out_sim05
