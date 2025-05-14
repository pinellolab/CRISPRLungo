import os, sys
from subprocess import Popen, PIPE

def confirm_insertion_seq(mutation_dict, ref_seq, ori_ref_name, possible_ref_path, output_dir, threads, get_ins_len = 20):

    print(f'Analyzing the insertions (>{get_ins_len}) ...')
    if not os.path.exists(output_dir + '/insert_analysis/'):
        os.mkdir(output_dir + '/insert_analysis/')
    fw = open(output_dir + '/insert_analysis/reference_for_ins.fasta', 'w')

    fw.write(f'>{ori_ref_name}\n{ref_seq}\n')
    fw.write(''.join(open(possible_ref_path).readlines()))
    fw.close()

    p = Popen(f'minimap2 -d {output_dir}/insert_analysis/reference_for_ins.mmi {output_dir}/insert_analysis/reference_for_ins.fasta', shell=True, stderr=PIPE, stdout=PIPE).communicate()

    fw = open(output_dir + '/insert_analysis/large_insertion.fasta', 'w')

    for read_id, muts in mutation_dict.items():
        for mut_n, mut in enumerate(muts[0]):
            if mut[0] == 'insertion' and len(mut[3]) > get_ins_len:
                fw.write(f'>{read_id}_{mut_n}_{mut[5]}_{mut[6]}\n{mut[3]}\n')

    
    fw.close()

    p = Popen(f'minimap2 -ax map-ont -t {threads} {output_dir}/insert_analysis/reference_for_ins.mmi {output_dir}/insert_analysis/large_insertion.fasta -o {output_dir}/insert_analysis/large_insertion.sam', shell=True, stderr=PIPE, stdout=PIPE).communicate()
    if p[1].decode('utf-8').find('ERROR') != -1:
        print(p[1].decode('utf-8'))
        sys.exit()
    
    with open(f'{output_dir}/insert_analysis/large_insertion.sam') as f:
        for line in f:
            if line[0] == '@':
                continue
            line_sp = line.split()
            if line_sp[5] == '*' or int(line_sp[4]) < 40:
                continue
            if line_sp[1] not in ['0', '16']:
                continue
            read_id = line_sp[0].split('_')
            read_id = [int(i) for i in read_id]
            reads = [line_sp]
            inversion_check = False
            query_seq = line_sp[9]
            query_seq_len = len(query_seq)

            if line.find('SA:Z:') != -1:
                for i in line_sp:
                    if i[:5] == 'SA:Z:':
                        SA_tag = i  
                        break

                while len(reads) < len(SA_tag.split(';')):
                    
                    line_sec = next(f).split('\t')
                    flag = int(line_sec[1])
                    if flag & 0x800:
                        reads.append(line_sec)

            muts = []

            for line_sp in reads:
                flag = int(line_sp[1])
                s = ''
                align_len = 0
                clip_len = [0, 0]
                front_clip = True
                for i in line_sp[5]:
                    if i in '0123456789':
                        s += i
                    else:
                        s = int(s)
                        if i in 'MD':
                            align_len += s
                        if i in 'SH':
                            if front_clip:
                                clip_len[0] = s
                            else:
                                clip_len[1] = s
                        front_clip = False
                        s = ''
                if align_len < 20:
                    continue
                ref_name = line_sp[2]
                strand = 1
                if flag & 0x10:
                    clip_len = clip_len[::-1]
                    strand = -1
                if ref_name == ori_ref_name and strand == -1:
                    inversion_check = True
                    ref_name = 'inversion'
                ref_start = int(line_sp[3])
                ref_end = ref_start + align_len
                muts.append((ref_name, read_id[2], read_id[3], ref_start, ref_end, strand))
            
            before_mut = list(mutation_dict[read_id[0]][0][read_id[1]])
            if inversion_check:
                before_mut[0] == 'inversion'
            before_mut.append(tuple(muts))

            mutation_dict[read_id[0]][0][read_id[1]] = before_mut


    return mutation_dict


                    
                    

    

