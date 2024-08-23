#!/usr/bin/env python
# coding: utf-8

# In[9]:


import argparse, sys
import time, os
from subprocess import Popen, PIPE


# In[23]:


import statistical_integration_clean as regular_py
import Post_process_final as visual
import CRISPRlungo_umi


# In[25]:





# In[ ]:


def main():

    start_time = time.time()

    if len(sys.argv) == 1:
        print('usage: CRISPRlungo [regular | umi]')
        sys.exit()

    parser = argparse.ArgumentParser(description='Analysis CRISPR mutaion pattern for Long-read sequencing')
    subparsers = parser.add_subparsers(dest='command')

    regular_parser = subparsers.add_parser('regular', help='Analyze CRISPR mutation for Long-read sequencing data without UMI')
    regular_parser.add_argument('ref', type=str,  help='Reference FASTA file')
    regular_parser.add_argument('control', type=str,  help='Control FASTQ file of Long-read sequencing data')
    regular_parser.add_argument('treated', type=str,  help='Treated FASTQ file of Long-read sequencing data')
    regular_parser.add_argument('output', type=str, help='Output file name')
    regular_parser.add_argument('-t', '--threads', type=int, default=1, help='Output file name')
    regular_parser.add_argument('--target', type=str, default=None, help='Target sequence')

    umi_parser = subparsers.add_parser('umi', help='Analyze CRISPR mutation for Long-read sequencing data with UMI')
    umi_parser.add_argument('ref', type=str,  help='Reference FASTA file')
    umi_parser.add_argument('input_file', type=str,  help='FASTQ file of Long-read sequencing data')
    umi_parser.add_argument('output_dir', type=str,  help='Directory for output files')
    umi_parser.add_argument('target', type=str,  help='Target sequence without PAM')
    umi_parser.add_argument('-t', '--threads', type=int, default=1, help='Output file name')
    umi_parser.add_argument('-c', '--clust_cutoff', type=int, default=5, help='The minimum of UMI cluster size')
    umi_parser.add_argument('--cleavage_pos', type=int, default=16, help='Cleavage position in target sequence')
    umi_parser.add_argument('--window', type=int, default=10, help='The range for mutation analysis around cleavage site')

    
    args = parser.parse_args()


    def run_minimap2(ref_file, input_file, output_file, longjoin_bandwidth, chaining_bandwidth, threads):
        align_st_time = time.time()
        print(f'minimap2 aligning {input_file} ...\r', end='')
        p = Popen(f'minimap2 -ax map-ont -t {threads} -r {chaining_bandwidth},{longjoin_bandwidth} {ref_file} {input_file} -o {output_file}', shell=True, stderr=PIPE, stdout=PIPE).communicate()
        if p[1].decode('utf-8').find('ERROR') != -1:
            print(p[1].decode('utf-8'))
            sys.exit()
        else:
            print(f'minimap2 aligning {input_file} ... Done {time.time() - align_st_time} s')
        
    def create_dir(dir_name):
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)
    
    if args.command == 'regular':
        
        print('Start to analysis using REGULAR mode ...')

               
        output_dir = os.path.dirname(args.output)
        if output_dir == '':
            output_dir = './'
        threads = args.threads
        
        # Run Miniseq2 to align received files
        ref_file = args.ref
        ref_seq = ''.join(open(ref_file).readlines()[1:]).replace('\n','').upper()
        ref_len = len(ref_seq)
        if args.target != None:
            target = args.target.upper()
            if ref_seq.find(target) != -1:
                cv_pos = ref_seq.find(target) + 16
            elif ref_seq.find(target.translate(target.maketrans('ATGC','TACG'))[::-1]) != -1:
                cv_pos = ref_seq.find(target.translate(target.maketrans('ATGC','TACG'))[::-1]) + 3
        longjoin_bandwidth = int(ref_len * 0.3)
        if longjoin_bandwidth < 500:
            chaining_bandwidth = longjoin_bandwidth
        else:
            chaining_bandwidth = 500

        ref_indexed_file = os.path.splitext(ref_file)[0] + '.mmi'
        Popen(f'minimap2 -d {ref_indexed_file} {ref_file}', shell=True, stderr=PIPE, stdout=PIPE).communicate()

        run_minimap2(ref_indexed_file, args.control,  output_dir + '/Control_alignment.sam', longjoin_bandwidth, chaining_bandwidth, threads)
        run_minimap2(ref_indexed_file, args.treated,  output_dir + '/Treated_alignment.sam', longjoin_bandwidth, chaining_bandwidth, threads)


        # Run Statistical anlaysis.py

        edited_dictionary, controll_dictionary, List_of_valid_IDs = regular_py.analysis_function(output_dir + '/Control_alignment.sam', output_dir + '/Treated_alignment.sam', ref_file)
        regular_py.process_mutations(edited_dictionary, args.output + '.tsv', List_of_valid_IDs)


    elif args.command == 'umi':

        print('Start to analyze using UMI mode ...')
        
        ref_file = args.ref
        output_dir = args.output_dir
        input_file = args.input_file
        threads = args.threads

        # Get UMI position information
        umi_pos, umi_len, ref_name, ref_seq = CRISPRlungo_umi.input_organize(ref_file, output_dir,input_file)
        ref_len = len(ref_seq)

        # Make output folders
        create_dir(output_dir)
    
        ref_dir = f'{output_dir}/ref_seq'
        create_dir(ref_dir)
        
        with open(ref_dir + '/ref.fasta', 'w') as f:
        	f.write(f'>{ref_name}\n{ref_seq}')
        
        create_dir(output_dir + '/align')

        # Align recieved FASTQ file
        longjoin_bandwidth = int(ref_len * 0.3)
        if longjoin_bandwidth < 500:
            chaining_bandwidth = longjoin_bandwidth
        else:
            chaining_bandwidth = 500
            
        p = Popen(f'minimap2 -d {ref_dir}/ref.mmi {ref_dir}/ref.fasta', shell=True, stderr=PIPE, stdout=PIPE).communicate()
        run_minimap2(ref_dir + "/ref.mmi", args.input_file,  output_dir + '/align/fastq.sam', longjoin_bandwidth, chaining_bandwidth, threads)

        # Get UMI sequence
        create_dir(f'{output_dir}/demultiplexing')

        CRISPRlungo_umi.extract_index_umi(output_dir + '/align/fastq.sam', ref_file, output_dir, umi_pos[0], umi_pos[1])

        # Clustering UMI sequence

        create_dir(f'{output_dir}/clustering')
        create_dir(f'{output_dir}/consensus')

        index_info = 'result,NNNNNNNNNN' #For index demultiplexing, not yet supported.

        index_info_list = index_info.split(',')
        index_names = []

        for i in range(int(len(index_info_list)/2)):
            
            i = index_info_list[2*i]
            index_names.append(i)
            
            create_dir(f'{output_dir}/clustering/{i}')
            create_dir(f'{output_dir}/clustering/{i}/1st_clusters')
            create_dir(f'{output_dir}/clustering/{i}/2nd_clusters')
            create_dir(f'{output_dir}/clustering/{i}/medaka_input')
            create_dir(f'{output_dir}/consensus/{i}')
            
        create_dir(f'{output_dir}/results')

        for fn in index_names:
            
            CRISPRlungo_umi.clustering_umi( f'{output_dir}/demultiplexing/umi_{fn}.fasta', umi_len, f'{output_dir}/clustering/{fn}', f'{output_dir}/consensus/{fn}', args.clust_cutoff, threads)

            run_minimap2(ref_dir + "/ref.mmi", f'{output_dir}/consensus/{fn}/consensus.fasta', f'{output_dir}/align/consensus_result.sam', longjoin_bandwidth, chaining_bandwidth, threads)

            # Mutation Analysis
            print('Mutation analysis...')
            create_dir(f'{output_dir}/results/{fn}')
            cv_pos = CRISPRlungo_umi.mutation_analysis(ref_seq, args.target.upper(), args.cleavage_pos, args.window, f'{output_dir}/align/consensus_{fn}.sam', f'{output_dir}/results/{fn}', threads=threads)

    #visualization

    if args.command == 'regular':
        visual.regular_statistic_plot(args.output + '.tsv', args.treated, output_dir)
        visual.regular_accuracy_plot(ref_seq, read_per_position, output_dir)
        read_per_position = visual.visualization_preprocess(output_dir + '/Treated_alignment.sam', ref_file)
        tsv_file = args.output + '.tsv'

    else:
        tsv_file = output_dir + '/results/result/read_classification.txt'
        read_per_position = visual.visualization_preprocess(output_dir + '/align/consensus_result.sam', ref_file)

    visual.base_proportion(read_per_position, output_dir, ref_seq)
    visual.mutation_pie_chart(tsv_file, output_dir)
    visual.indel_per_position(tsv_file, ref_seq, output_dir)
    visual.Insertion_length(tsv_file, output_dir)
    visual.Deletion_length(tsv_file, output_dir)
    visual.Deletion_count_length(tsv_file, output_dir)
    if args.target != None:
        visual.allele_plot(ref_seq, cv_pos, tsv_file, output_dir)
    #LD_tornado()
      

        


# In[ ]:


if __name__=='__main__':
    main()

