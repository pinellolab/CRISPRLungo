# CRISPR-LongReads


![image](https://github.com/user-attachments/assets/842d4741-c261-4bce-8e36-709f7ae3c76f)



### Requirements
- python3 >= 3.8
- Minimap2 = 2.28
- biopython = 1.85
- pysam
- scipy
- vsearch = 2.30.0
- matplotlib
- plotly >= 6.0.1
- pandas
- numpy
- editdistance
- pyspoa
- kaleido  (conda install conda-forge::python-kaleido)


### Install

After download repository,

``` bash
     cd CRISPR-LongReads
     conda env create -n {env_name} -f environment.yml # can use mamba instead of conda
     pip install -e .
```

### Example

1. REGULAR mode

```bash
    CRISPRlungo PD1.fasta --control Nanopore_umi_Run_test_control_wo_chi.fastq Nanopore_umi_Run_test_wo_chi.fastq regular_output ggcgccctggccagtcgtct
```

2. UMI mode

```bash
    CRISPRlungo --umi PD1_umi.fasta --control Nanopore_umi_Run_test_control_wo_chi.fastq Nanopore_umi_Run_test_wo_chi.fastq regular_output ggcgccctggccagtcgtct
```



