import time, pysam, sys, os
from subprocess import Popen, PIPE


def reverse_complementary(s):
	return s.translate(s.maketrans('ATGCatgc','TACGtacg'))[::-1]

def run_minimap2(ref_file, input_file, output_file, longjoin_bandwidth, chaining_bandwidth, threads):
		align_st_time = time.time()
		print(f'minimap2 aligning {input_file} ...\r', end='')
		p = Popen(f'minimap2 -ax map-ont -t {threads} -p 0.5 -r {chaining_bandwidth},{longjoin_bandwidth} {ref_file} {input_file} -o {output_file}', shell=True, stderr=PIPE, stdout=PIPE).communicate()
		if p[1].decode('utf-8').find('ERROR') != -1:
			print(p[1].decode('utf-8'))
			sys.exit()
		else:
			print(f'minimap2 aligning {input_file} ... Done {round(time.time() - align_st_time, 2)} s')			
	
def soft_clipped(align_file_1, output_file, cut_len_threshold, fasta_check=False):

    check_SA = False
    read_d = []
    with pysam.AlignmentFile(align_file_1) as f1, open(output_file, 'w') as fw:
        for line in f1:
            if line.is_unmapped or line.is_secondary or line.is_supplementary:
                continue
            if line.mapping_quality < 30 and align_file_1.split('/')[-1].split('_')[0] != '1':
                continue
            clipped_seq = [[],[]]
            cigar = line.cigar
            read_d.append(line)
            if 'alignStrandInfo' not in line.query_name:
                line.query_name += f'_alignStrandInfo'
            
            if line.is_reverse:
                strand = -1
            else:
                strand = 1
            line.query_name += f'_{strand}'
            
            if fasta_check == False:
                if cigar[0][0] == 4 and cigar[0][1] > cut_len_threshold:
                    fw.write(f'@{line.query_name}\n{line.query_sequence[:cigar[0][1]]}\n+\n{line.qual[:cigar[0][1]]}\n')
                    check_SA = True
                if cigar[-1][0] == 4 and cigar[-1][1] > cut_len_threshold:
                    fw.write(f'@{line.query_name}\n{line.query_sequence[-cigar[-1][1]:]}\n+\n{line.qual[-cigar[-1][1]:]}\n')
                    check_SA = True
            else:
                if cigar[0][0] == 4 and cigar[0][1] > cut_len_threshold:
                    fw.write(f'>{line.query_name}\n{line.query_sequence[:cigar[0][1]]}\n')
                    check_SA = True
                if cigar[-1][0] == 4 and cigar[-1][1] > cut_len_threshold:
                    fw.write(f'>{line.query_name}\n{line.query_sequence[-cigar[-1][1]:]}\n')
                    check_SA = True

            
    return check_SA, read_d

def run_triple_minimap2(ref_file, input_file, output_file, longjoin_bandwidth, chaining_bandwidth, threads, len_cutoff = 200, fasta_check=False):
    output_path = output_file[:output_file.rfind('/')]
    check_SA = True
    n = 0
    read_d = {}
    file_format = 'fastq'
    if fasta_check:
        file_format = 'fasta'
    while check_SA == True:
        n += 1
        if n == 1:
            fastq_file = input_file
        else:
            fastq_file = output_path + f'/{n}_soft.{file_format}'

        run_minimap2(ref_file, fastq_file, output_path + f'/{n}_align.sam' , longjoin_bandwidth, chaining_bandwidth, threads)
        check_SA, read_out = soft_clipped(output_path + f'/{n}_align.sam', output_path + f'/{n+1}_soft.{file_format}', len_cutoff, fasta_check)
        for read in read_out:
            query_name = read.query_name
            strand_info = query_name.split('_alignStrandInfo_')[1].split('_')
            strand = 1
            for x in strand_info:
                strand *= int(x)
            if strand == -1:
                read.is_reverse = True
            else:
                read.is_reverse = False
            query_name = query_name[:query_name.find('_alignStrandInfo_')]
            if n != 1:
                read.is_supplementary = True
            if query_name not in read_d:
                read_d[query_name] = [read]
            else:
                ori_seq = read_d[query_name][0].query_sequence
                ori_strand = -1 if read_d[query_name][0].is_reverse else 1
                if ori_strand * strand == -1:
                    part_st = ori_seq.find(reverse_complementary(read.query_alignment_sequence))
                    part_ed = ori_seq.find(reverse_complementary(read.query_alignment_sequence)) + len(read.query_alignment_sequence)
                else:
                    part_st = ori_seq.find(read.query_alignment_sequence)
                    part_ed = ori_seq.find(read.query_alignment_sequence) + len(read.query_alignment_sequence)
                
                if not fasta_check: q = read.query_qualities

                if read.cigar[0][0] == 4:
                    read.query_sequence = read.query_alignment_sequence
                    if not fasta_check: q = q[read.cigar[0][1]:]
                    read.cigar = [(5, part_st)] + read.cigar[1:]
                    
                else:
                    read.cigar = [(5, part_st)] + read.cigar

                if read.cigar[-1][0] == 4:
                    read.query_sequence = read.query_alignment_sequence
                    if not fasta_check: q = q[:-read.cigar[-1][1]]
                    read.cigar = read.cigar[:-1] + [(5, len(ori_seq) - part_ed)]
                else:
                    read.cigar += [(5, len(ori_seq) - part_ed)]
                if not fasta_check: read.query_qualities = q

                read_d[query_name].append(read)
    fw = open(output_file, 'w') 
    with open(output_path + '/1_align.sam') as f:
        for line in f:
            if line[0] == '@':
                fw.write(line)
            else:
                break
    for x, reads in read_d.items():
        read_n = 0
        if len(reads) == 1:
            tags = reads[0].get_tags()
            filtered_tags = [t for t in tags if t[0] != "SA"]
            reads[0].set_tags(filtered_tags)
        else:
            reads[0].set_tag("SA", ('N;'*(len(reads)))[:-1], value_type='Z')
        for read in reads:
            read_n += 1
            fw.write(read.to_string() +f'\t\n')
    fw.close()
    
    # remove tmp files
    
    for fn in os.listdir(output_path):
         if fn[:-10] == '_align.sam' or fn[:-11] == '_soft.fasta':
            os.system(f'rm {output_path}/{fn}')

