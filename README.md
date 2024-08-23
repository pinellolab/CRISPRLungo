# CRISPR-LongReads


### Requirements

- biopython
- pysam
- scipy
- vsearch 
- matplotlib
- plotly
- pandas
- numpy
- editdistance
- pyspoa
- kaleido  (conda install conda-forge::python-kaleido)


### Example

1. REGULAR mode

```bash
    python3 CRISPRlungo.py regular data/PD1.fasta data/Nanopore_umi_Run_test_control_wo_chi.fastq data/Nanopore_umi_Run_test_wo_chi.fastq regular_output --target ggcgccctggccagtcgtct
```

2. UMI mode

```bash
    python3 CRISPRlungo.py umi data/PD1_umi.fasta data/Nanopore_umi_Run_test_wo_chi.fastq output_umi/ ggcgccctggccagtcgtct
```


![CRISPRlungo_image](https://github.com/user-attachments/assets/1b9a03b2-21d6-4450-b8c3-b7ac8afe57e6)

