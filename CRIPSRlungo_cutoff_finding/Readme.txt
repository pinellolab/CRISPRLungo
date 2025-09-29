This is the code used for generating figure S10
It can also be used for setting a personalized threshold based on user untreated data

Required packages:
1)library(MASS)
2)library(ggplot2)
3)library(scales)
4)library(knitr)
5)library(kableExtra)

Input:
1)File of allele lengths extracted from allele_table.txt (in CRISPRlungo output)

Output: 
1)Graphs with allele length on X-axis and count in Y-axis (trans - Y-axis is log scaled)
2)P-value of allele in control distribution being larger than threshold
3)Table of P-value of allele being exactly a certain length from 1-10

There is an example file of allele lengths taken from gangbao_BCL11A_CONTROL_ins_lengths