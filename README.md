# CRISPRlungo

![Logo](https://github.com/pinellolab/CRISPR-LongReads/blob/main/docs/title_icon.png)

**CRISPRlungo** is a software pipeline designed to analyze genome editing outcomes using long-read sequencing data.  
It supports multiple CRISPR platforms (base editors, prime editors) and is compatible with various sequencing methods such as amplicon sequencing, UMI-tagged long-read sequencing, and nCATS.


---

## Pipeline Overview
1. Align sequencing reads to the reference genome, filter out low-quality reads, and remove chimeric reads.
2. If UMIs are used, cluster UMI-tagged reads and generate consensus sequences.
3. If control samples are provided, perform background error filtering using statistical analysis.
4. Quantify small indels, large indels, and inversions.
5. Use the submodule **CRISPRlungoAllele** to classify allele groups and identify PCR-induced chimeric reads.


---

## What can CRISPRlungo do?
- Filtering of low-quality reads  
- Identification and removal of chimeric reads generated during library preparation  
- Alignment using an optimized pipeline to detect structural variants and inversions  
- UMI extraction, clustering, and consensus read generation (if UMIs are present)  
- Background error estimation and filtering using control data (if available)  
- Quantification of small insertions/deletions (indels), large indels, inversions, and sequence integrations  
- Detection and quantification of intended mutations when a reference for the edited sequence is provided  
- Visualization of:
  - Indel size and position distributions  
  - Substitution patterns and their positions  
  - Allele frequency and edit spectrum  


---

## Installation
```bash
git clone https://github.com/pinellolab/CRISPRlungo
cd CRISPRlungo
conda env create -n {env_name} -f environment.yml  # use mamba instead of conda if preferred
conda activate {env_name} # use mamba instead of conda if preferred
pip install -e .
CRISPRlungo -h
```


---
## Usage

CRISPRlungo has **4 usage modes**:

### 1. Default
```bash
CRISPRlungo {reference.FASTA} {sequencing_result.FASTQ} {Output directory} {target sequence}
```

### 2. Background filtering
```bash
CRISPRlungo {reference.FASTA} {sequencing_result.FASTQ} {Output directory} {target sequence} --control {control.FASTQ}
```

### 3. UMI option
The UMI context must be annotated in the reference FASTA using parentheses ( and ):

>PD1
GGTGCTGAAGAAAGTTGTCGGT...AGGTTT(VVVVTTVVVVTTVVVVTTVVVV)TTTGGGACACCG...

```bash
CRISPRlungo --umi {reference.FASTA} {sequencing_result.FASTQ} {Output directory} {target sequence}
```

### 4. UMI + background filtering
```bash
CRISPRlungo --umi {reference.FASTA} {sequencing_result.FASTQ} {Output directory} {target sequence} --control {control.FASTQ}
```

---

## Example Runs

### Example 1: Background error filter
```bash
cd data
CRISPRlungo PD1.fasta --control Nanopore_umi_Run_test_control_wo_chi.fastq \
    Nanopore_umi_Run_test_wo_chi.fastq regular_output ggcgccctggccagtcgtct
```
* Generates regular_output/
* Results available at: regular_output/combined_graphs.html

### Example 2: UMI + background error filter
```bash
cd data
CRISPRlungo --umi PD1_umi.fasta --control Nanopore_umi_Run_test_control_wo_chi.fastq \
    Nanopore_umi_Run_test_wo_chi.fastq umi_output ggcgccctggccagtcgtct
```
* Generates umi_output/
* Results available at: umi_output/combined_graphs.html


---

## Paramters
* --umi : Enable UMI mode
* --control : Control FASTQ for background error filtering
* --cleavage_pos : Cleavage position relative to target (default: 16)
* --additional_target : Add extra target sequences
* --window : Window size around cleavage site
* --whole_window_between_targets : Include entire region between two targets
* --induced_sequence_path : Desired induced sequence reference file
* --integration_file : FASTA file with potential integration sequences
* --merge_substitution : Treat consecutive substitutions as one event
* --min_read_cnt : Minimum read count filter
* --min_read_freq : Minimum frequency filter
* --mix_tag : Handle multiple mutations (False: prioritize indels, True: complex)
* --min_mut_freq_no_control_refmut : Minimum mutation frequency (no cotrol/reference)
* --induced_paritial_similiarity : Similarity threshold for partial induced classification
* --range_both_end_region : Range filter for short fragments
* --align_sa_len_threshold : Minimum alignment length for soft-clipped/supplementary reads
* --p_value_threshold : p-value cutoff
* --mut_freq_threshold : Mutation frequency threshold
* -c/--clust_cutoff : Minimum UMI cluster size
* --just_visualization : Only visualization (skip consensus generation)
* --allele_plot_window : Window size for allele plots (default: window+10)
* --allele_plot_lines : Number of representative sequences in allele plots
* --show_all_between_allele : Show all sequences spanning between targets
* -t/--threads : Number of threads

---
## Result Files
* combined_graphs.html – summary of all plots/results
* Input_summary.txt – options used in the run
* Allele_table.txt – allele classification with mutations, counts, frequencies
* Edit_pattern_count.txt – mutation pattern counts/frequencies
* Mutation_pattern_p_values.txt – statistical test results for patterns
* Mutation_summary_count.txt – summary of mutation patterns
* Preprocess_count.txt – preprocessing read filtering summary
* Read_classification.txt – mutation classification per read

---
## Subanalysis Tool - CRISPRlungoAllele
**CRISPRlungoAllele** performs post-analysis of CRISPRlungo results by classifying alleles into multiple groups based on mutation type.

This enables:
* Identification of mutation patterns within alleles
* Detection of similar/homologous sequences
* Estimation/removal of PCR chimeric products

## Usage
```bash
CRISPRlungoAllele {analysis_file_dir} {custom_category_file}
```
* analysis_file_dir : output directory from CRISPRlungo run
* custom_category_file : tab-separated file with mutation conditions

### Example:
```bash
 cd data
 CRISPRlungo Allele_ref.fasta Allele_Cas_treated.fastq Allele_example --control Allele_Mock.fastq --threads 8 --induced_sequence_path induced.fasta actgaaatctgtaagcaggt
 CRISPRlungoAllele -i Allele_mutations.txt Allele_example
```

---
## Custom Category File Format
* Column 1 : Category name (or AND for multi-condition)
* Column 2-3 : Start and end positions
* Column 4 : Mutation type (WT, Sub, Deletion)
* Column 5 : Threshold (e.g., 403/426/WT/0.8)
* Column 6 : Include outside window? (TRUE/FALSE)
* Column 7 : Condition required? (TRUE/FALSE)
* Column 8 : Frequency calculation within group (TRUE/FALSE)
* Column 9 : [reserved]
* Column 10 : Number of top alleles to show


---
## CRISPRlungoAllele Options
* --min_read_cnt : Minimum read count
* --min_read_freq : Minimum read frequency
* --allele_plot_window : Window size for allele plots (default: window+10)
* --show_all_between_allele : Show all sequences spanning between targets


---
## CRISPRlungoAllele Result Files

Located in custom_results/ inside the analysis directory:
* allele_plot.png
* custom_mutation_allele_plot_input.txt
* custom_mutation_pattern_count.txt
Each category subfolder contains:
* deletion_and_insertions_per_position.png
* mutation_pattern_count.txt
* Mutation_pie_chart.png
* mutation_summary_count.txt
* pattern_pie_chart.png
* percent_of_alleles_pie_chart.png
* read_classification.txt



