#!/usr/bin/env bash
# 20261111 — CRISPRlungoAllele / CRISPRlungoBatch compatibility check for PA mode.
# Run from the directory holding ref/ and fq/.
set -e
SRC=${SRC:-src}            # CRISPRlungo source directory
PY=${PY:-python}

# 1. the same sample with and without PA
$PY $SRC/CRISPRlungo.py ref/amp_WT.fasta fq/s20.fastq out_pa \
    TATTGAGAAGCCGAAACCTC --additional_target TGGGAAGTGCTGCTGACCGC \
    --control fq/ctrl19.fastq --whole_window_between_targets --allele_plot_lines 10 \
    --pa_donor ref/donor.fasta
$PY $SRC/CRISPRlungo.py ref/amp_WT.fasta fq/s20.fastq out_nopa \
    TATTGAGAAGCCGAAACCTC --additional_target TGGGAAGTGCTGCTGACCGC \
    --control fq/ctrl19.fastq --whole_window_between_targets --allele_plot_lines 10

# 2. the main output must not move
for f in mutation_summary_count.txt read_classification.txt Allele_table.txt mutation_patter_count.txt; do
    cmp out_pa/results/$f out_nopa/results/$f && echo "SAME $f"
done

# 3. CRISPRlungoAllele on the main directory, then on every PA allele
$PY $SRC/CRISPRlungoAllele.py out_pa   -i categories.txt
$PY $SRC/CRISPRlungoAllele.py out_nopa -i categories.txt
for f in $(cd out_pa/custom_results && find . -name '*.txt'); do
    cmp out_pa/custom_results/$f out_nopa/custom_results/$f && echo "SAME $f"
done
for a in wt pa dual f5 f3; do
    $PY $SRC/CRISPRlungoAllele.py out_pa/pa_results/reference_$a -i categories.txt
done

# 4. a mixed batch, donor path written relative to the batch file
$PY $SRC/CRISPRlungoBatch.py -b batch.txt -o batch_out -t 4
cut -f1,5,6,28,29,35 batch_out/CRISPRlungoBatch_summary.txt
