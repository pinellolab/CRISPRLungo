import os, sys, edlib
from subprocess import Popen, PIPE
from CRISPRlungo_minimap import *
import CRISPRlungo_mutation_analysis as mutation_analysis

def rev_comp(s): return s.translate(s.maketrans('ATGC','TACG'))[::-1]

def align_edlib(seq_info, ref_dict, ori_ref_name, max_dist = 0.05):
    res_dict = {}
    insert_seq = seq_info[3]
    insert_seq_len = len(insert_seq)
    insert_seq_r = rev_comp(insert_seq)
    for ref_id, ref_seq in ref_dict.items():
        res_dict[ref_id] = [False, ] #Edit distance, pos_st, pos_ed, strand
        if insert_seq_len <= len(ref_seq):
            align_res = edlib.align(insert_seq, ref_seq, mode="HW", task="locations")
            if align_res['editDistance'] < insert_seq_len * max_dist:
                res_dict[ref_id] = [align_res['editDistance'], align_res['locations'][0][0], align_res['locations'][0][1], '+']
            align_res = edlib.align(insert_seq_r, ref_seq, mode="HW", task="locations")
            if align_res['editDistance'] < insert_seq_len * max_dist:
                if res_dict[ref_id][0] == False or align_res['editDistance'] < res_dict[ref_id][0]:
                    res_dict[ref_id] = [align_res['editDistance'], align_res['locations'][0][0], align_res['locations'][0][1], '-']
        else:
            align_res = edlib.align(ref_seq, insert_seq, mode="HW", task="locations")
            if align_res['locations'][0][1] - align_res['locations'][0][0] > 20:
                trimmed_insert = insert_seq[align_res['locations'][0][0]: align_res['locations'][0][1]]
                align_res = edlib.align(trimmed_insert, ref_seq, mode="HW", task="locations")
                if align_res['editDistance'] < insert_seq_len * max_dist:
                    res_dict[ref_id] = [align_res['editDistance'], align_res['locations'][0][0], align_res['locations'][0][1], '+'] 
            align_res = edlib.align(ref_seq, insert_seq_r, mode="HW", task="locations")
            if align_res['locations'][0][1] - align_res['locations'][0][0] > 20:
                trimmed_insert = insert_seq_r[align_res['locations'][0][0]: align_res['locations'][0][1]]
                align_res = edlib.align(trimmed_insert, ref_seq, mode="HW", task="locations")
                if align_res['editDistance'] < insert_seq_len * max_dist:
                    if res_dict[ref_id][0] == False or align_res['editDistance'] < res_dict[ref_id][0]:
                        res_dict[ref_id] = [align_res['editDistance'], align_res['locations'][0][0], align_res['locations'][0][1], '-'] 

    res = sorted(res_dict.items(), key = lambda x: x[1][0])
    muts = []
    inversion_check = False
    for i in res:
        if i[1][0] == False:
            continue
        i = list(i)
        if i[0] == ori_ref_name and i[1][3] == '-':
            i[0] = 'inversion'
            inversion_check = True
        muts.append((i[0], seq_info[5], seq_info[6], i[1][1], i[1][2], i[1][3]))
    return muts, inversion_check


    

def confirm_insertion_seq(mutation_dict, ref_seq, ori_ref_name, possible_ref_path, output_dir, threads, anchor_validation_only=False, get_ins_len = 20, minimap_len_limit=100):

    print(f'Analyzing the insertions (>{get_ins_len}) ...')
    if not os.path.exists(output_dir + '/insert_analysis/'):
        os.mkdir(output_dir + '/insert_analysis/')

    ref_dict = {}
    fw = open(output_dir + '/insert_analysis/reference_for_ins.fasta', 'w')
    fw.write(f'>{ori_ref_name}\n{ref_seq}\n')
    ref_dict[ori_ref_name] = ref_seq.upper()
    if possible_ref_path != False:
        f_ref = open(possible_ref_path).readlines()
        for i in range(int(len(f_ref)/2)):
            insref_name = f_ref[i*2].strip()
            insref_seq = f_ref[i*2+1].strip()
            ref_dict[insref_name] = insref_seq
            fw.write(f'{insref_name}\n{insref_seq}\n')
    fw.close()

    #p = Popen(f'minimap2 -d {output_dir}/insert_analysis/reference_for_ins.mmi {output_dir}/insert_analysis/reference_for_ins.fasta', shell=True, stderr=PIPE, stdout=PIPE).communicate()

    fw = open(output_dir + '/insert_analysis/large_insertion.fasta', 'w')

    for read_id, muts in mutation_dict.items():
        for mut_n, mut in enumerate(muts[0]):
            if mut[0] == 'insertion' and len(mut[3]) > get_ins_len:
                if len(mut[3]) < minimap_len_limit:
                    align_res, inv_check = align_edlib(mut, ref_dict, ori_ref_name)
                    if align_res != []:
                        before_mut = list(mutation_dict[read_id][0][mut_n])
                        if inv_check:
                            before_mut[0] == 'inversion'
                        before_mut.append(tuple(align_res))
                        mutation_dict[read_id][0][mut_n] = before_mut
                else:
                    fw.write(f'>{read_id}_{mut_n}_{mut[5]}_{mut[6]}\n{mut[3]}\n')

    
    fw.close()

    p = Popen(f'minimap2 -ax map-ont -t {threads} {output_dir}/insert_analysis/reference_for_ins.fasta {output_dir}/insert_analysis/large_insertion.fasta -o {output_dir}/insert_analysis/large_insertion.sam', shell=True, stderr=PIPE, stdout=PIPE).communicate()
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
                sig_mut = False
                for i in line_sp[5]:
                    if i in '0123456789':
                        s += i
                    else:
                        s = int(s)
                        if i in 'DI' and s > 20:
                            sig_mut = True
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
                muts.append((ref_name, read_id[2], read_id[3], ref_start, ref_end, strand, sig_mut))
            
            before_mut = list(mutation_dict[read_id[0]][0][read_id[1]])
            if inversion_check:
                before_mut[0] == 'inversion'
            before_mut.append(tuple(muts))

            mutation_dict[read_id[0]][0][read_id[1]] = before_mut


    return mutation_dict

def confirm_induced_ins_with_simulation(mutation_dict, 
                                        List_of_valid_IDs,
                                        ref_seq, 
                                        induced_sequence_path,
                                        simulation_mod,
                                        output_dir, 
                                        threads, 
                                        cv_pos,
                                        cv_pos_2,
                                        window,
                                        induced_mutations,
                                        anchor_validation_only=False, 
                                        get_ins_len = 20,
                                        minimap_len_limit=100):
    
    induce_seq = ''.join(open(induced_sequence_path).readlines()[1:]).replace('\n','').split('>')[0].upper()
    fw = open(f'{output_dir}/ref_seq/induce_ref.fasta', 'w')
    fw.write(f'>induce_seq\n{induce_seq}\n')

    if cv_pos_2:
        cv_pos_2 = len(induce_seq) - (len(ref_seq) - cv_pos_2)
    else:
        cv_pos_2 = len(induce_seq) - (len(ref_seq) - cv_pos)


    print('Generate induced sequence control simulation file using badread...')
    p = Popen(f'minimap2 -d {output_dir}/ref_seq/induce_ref.mmi {output_dir}/ref_seq/induce_ref.fasta', shell=True, stderr=PIPE, stdout=PIPE).communicate()

    try:
        p = Popen(["badread", "--help"], stdout=PIPE, stderr=PIPE)
        p.communicate()
    except FileNotFoundError:
        print('ERROR!: need to install badread')
        sys.exit()

    fw = open(f'{output_dir}/simulation/induce_simul.fasta', 'w')

    badread_identity = False
    badread_model = False

    if simulation_mod == 'pacbio':
        badread_model = 'pacbio2021' 
        badread_identity = '30,3'
    elif simulation_mod == 'ont_sup':
        badread_model = 'nanopore2023' 
        badread_identity = '98.5,99.5,1'
    elif simulation_mod == 'ont_hac':
        badread_model = 'nanopore2023' 
        badread_identity = '97.5,99,2.5'
    elif simulation_mod == 'ont_fast':
        badread_model = 'nanopore2023' 
        badread_identity = '95,99,5'
    else:
        print('ERROR: --large_induced_insertion value should be among pacbio, ont_sup, ont_hac, ont_fast')
        sys.exit()
    
    for i in range(10000): fw.write(f'>{i}\n{induce_seq}\n')

    cmd = ["badread", "simulate",
        "--reference", f'{output_dir}/simulation/induce_simul.fasta',
        "--quantity", "1x",
        "--error_model", badread_model,
        "--qscore_model", badread_model,
        "--identity", badread_identity]
    
    """
    with open(f'{output_dir}/simulation/induce_simul.fastq', "w") as outfile:
        proc = Popen(cmd, stdout=outfile, stderr=PIPE, text=True)

        print('')
        for line in proc.stderr:
            print(line.strip(), end="\r")

        proc.wait()
    
    if proc.returncode == 0:
        print("badread simulation completed.")
    else:
        print("badread error:")
        print(proc.stderr.decode())"""


    longjoin_bandwidth = int(len(induce_seq) * 0.3)
    if longjoin_bandwidth < 500:
        chaining_bandwidth = longjoin_bandwidth
    else:
        chaining_bandwidth = 500

    
    run_triple_minimap2(f'{output_dir}/ref_seq/induce_ref.mmi ', f'{output_dir}/simulation/induce_simul.fastq',  
                        output_dir + '/align/induced_control_alignment.sam', 
                        longjoin_bandwidth, chaining_bandwidth, 
                        threads, len_cutoff=100, 
                        fasta_check=True)
    
    run_triple_minimap2(f'{output_dir}/ref_seq/induce_ref.mmi ', f'{output_dir}/align/Treated_separated.fasta',  
                        output_dir + '/align/induced_treated_alignment.sam', 
                        longjoin_bandwidth, chaining_bandwidth, 
                        threads, len_cutoff=100, 
                        fasta_check=True)
    

    ind_edited_dictionary, controll_dictionary, ind_List_of_valid_IDs, control_reads_cnt, edited_reads_cnt = mutation_analysis.analysis_function(
				output_dir + '/align/induced_control_alignment.sam', 
				output_dir + '/align/induced_treated_alignment.sam', 
				f'{output_dir}/simulation/induce_simul.fasta', 
				output_dir + '/simulation/induced_ins/',
				cv_pos, 
				cv_pos_2, 
				window, 
				False,
				[], 
				[],
                -1)

    print(cv_pos)
    print(cv_pos_2)
    reverse_List_of_valid_IDs = {v: k for k, v in List_of_valid_IDs.items()} 

    for x, y in ind_edited_dictionary.items():
        if y[0] != []:
            continue
        read_id = ind_List_of_valid_IDs[x]
        print(read_id) 
        if read_id not in reverse_List_of_valid_IDs:
            print('None')
        input()

    









    



    

