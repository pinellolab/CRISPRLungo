import sys, argparse, os
import CRISPRlungo_visualization as visual
import CRISPRlungo_regular as regular_py
from CRISPRlungo_minimap import *
from itertools import combinations



def main():

    def convert_postion(pos):
        try:
            pos = int(pos)
            return pos
        except:
            if 'cv_pos1' in pos or 'cv_pos' == pos:
                if pos == 'cv_pos1' or pos == 'cv_pos':
                    return cv_pos
                elif '-' in pos:
                    return cv_pos - int(pos.split('-')[1])
                elif '+' in pos:
                    return cv_pos + int(pos.split('+')[1])
            elif 'cv_pos2' in pos:
                if pos == 'cv_pos2':
                    return cv_pos2
                elif '-' in pos:
                    return cv_pos2 - int(pos.split('-')[1])
                elif '+' in pos:
                    return cv_pos2 + int(pos.split('+')[1])
        return False
    
    def check_position(mut, start_pos, end_pos):
        mut = mut.split(':')[0].split('_')
        if int(mut[1]) < start_pos or int(mut[0]) > end_pos:
            return False
        return True
    
    def create_dir(dir_name):
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)

    def marker_find(mut_set):
        all_vals = set().union(*mut_set.values())
        for r in range(1, len(all_vals)+1):
            for combo in combinations(all_vals, r):
                marker_set = set(combo)
                sigs = {}
                for g, vals in mut_set.items():
                    sigs[g] = tuple(sorted(set(vals) & marker_set))
                unique_check = len(set(sigs.values())) == len(sigs)
                nonempty_check = all(len(sig) > 0 for sig in sigs.values())
                if unique_check:
                    return sigs
        return False
        
    def marker_find_mild(mut_set):
        all_vals = set().union(*mut_set.values())
        print(all_vals)
        for r in range(1, len(all_vals)+1):
            for combo in combinations(all_vals, r):
                marker_set = set(combo)
                sigs = {}
                for g, vals in mut_set.items():
                    sigs[g] = tuple(sorted(set(vals) & marker_set))
            return sigs
        return False

    def set_custom_mutation(input_file, cv_pos, cv_pos2, induced_mut, analysis_res_dir):
        mutation_category = []
        additional_analysis_cateogry = {}

        if input_file[input_file.rfind('.')+1:].upper() in ['FASTA', 'FNA', 'FA']:
            create_dir(analysis_res_dir + '/custom_results/align')
            run_triple_minimap2(f'{analysis_res_dir}/ref_seq/ref_wo_umi.mmi', args.induced_sequence_path,  f'{analysis_res_dir}/custom_results/align/induced_mutation_reference.sam', longjoin_bandwidth, chaining_bandwidth, 1, len_cutoff=args.align_sa_len_threshold)
            induced_mutations, induced_mutation_str = regular_py.get_induced_mutation(output_dir + '/induced_mutation_reference.sam', ref_file, cv_pos, cv_pos_2, args.window, window_between, range_align_end)
        else:
            with open(input_file) as f:
                mut_type = False
                for line in f:
                    if line[0] == '#' or line.strip() == '':
                        continue
                    
                    line = line.strip().split()

                    if line[0] == 'AND':
                        if not mut_type:
                            print('ERROR: "AND" can not be used in first line')
                            sys.exit()
                        pass
                    else:
                        mutation_category.append([line[0], []])
                        mut_type = line[0]
                
                    start_pos = convert_postion(line[1])

                    if not start_pos:
                        print('ERROR: not proper format at start_pos') 
                        sys.exit()
                    end_pos = convert_postion(line[2])
                    if not end_pos:
                        print('ERROR: not proper format at end_pos') 
                        sys.exit()
                    mut = []
                    mutation_marker = True
                    if line[5].upper() == 'FALSE':
                        mutation_marker = False
                    OrNot_marker = True
                    if line[6].upper() == 'FALSE':
                        OrNot_marker = False
                    Include_marker = True
                    if line[8].upper() == 'FALSE':
                        Include_marker = False

                
                    if line[7] == 'TRUE':
                        additional_analysis_cateogry[line[0]] = {}
                        for i in ['WT', 'Del', 'Ins', 'Sub', 'LargeDel', 'LargeIns', 'Inv', 'Complex']:
                            additional_analysis_cateogry[line[0]][i] = 0

                    if line[3].upper() == 'WT':
                        mutation_category[-1][1].append(['WT', [start_pos, end_pos, float(line[4])], mutation_marker, OrNot_marker, Include_marker])
                    elif line[3].upper() == 'DESIRED':
                        for i in induced_mut:
                            if check_position(i, start_pos, end_pos):
                                mut.append(i)
                        mutation_category[-1][1].append([mut, [0, 0, float(line[4])], mutation_marker, OrNot_marker, Include_marker])
                    elif 'DEL' in line[3].upper():
                        mutation_category[-1][1].append(['DEL', [start_pos, end_pos, float(line[4])], mutation_marker, OrNot_marker, Include_marker])
                    elif 'SUB' in line[3].upper():
                        mutation_category[-1][1].append(['SUB', [start_pos, 0, line[4], line[3].split('_')[1], line[3].split('_')[2]], mutation_marker, OrNot_marker, Include_marker])
                
        print(mutation_category)

        return mutation_category, additional_analysis_cateogry
    
    def set_mutation_reference(reference_file, cv_pos, cv_pos2, induced_mut, analysis_res_dir, mut_range_reference):
        
        mutation_category_dict = {}
        additional_analysis_cateogry = {}

        cv_pos_ed = cv_pos
        if cv_pos2 != False:
            cv_pos_ed = cv_pos2

        f = open(reference_file).readlines()
        ref_seq = ''.join(open(f'{analysis_res_dir}/ref_seq/ref_wo_umi.fasta').readlines()[1:]).replace('\n','')
        window_front_st = 100
        window_front_ed = mut_range_reference
        if window_front_ed > cv_pos - 100:
            window_front_ed = cv_pos - 100
        window_back_st = len(ref_seq) -  mut_range_reference
        window_back_ed = len(ref_seq) - 100
        if window_back_st < cv_pos_ed + 100:
            window_back_st = cv_pos_ed + 100

        front_mut = []
        back_mut = []

        induced_mut_front = {}
        induced_mut_back = {}


        for i in range(int(len(f)/2)):
            read_id = f[2*i].replace('>', '').strip()
            seq = f[2*i+1].strip()
            create_dir(f'{analysis_res_dir}/custom_results/align')
            fw = open(f'{analysis_res_dir}/custom_results/align/{read_id}.fasta', 'w')
            fw.write(f'>{read_id}\n{seq}\n')
            fw.close()
            longjoin = len(ref_seq)*0.3
            if longjoin < 500:
                chain = longjoin
            else:
                chain = 500
            run_triple_minimap2(f'{analysis_res_dir}/ref_seq/ref_wo_umi.mmi', f'{analysis_res_dir}/custom_results/align/{read_id}.fasta',  f'{analysis_res_dir}/custom_results/align/{read_id}.sam', longjoin, chain, 1, len_cutoff=100)
            mutations, mutation_str = regular_py.get_induced_mutation(f'{analysis_res_dir}/custom_results/align/{read_id}.sam', f'{analysis_res_dir}/ref_seq/ref_wo_umi.fasta', cv_pos, False, window, False, 100, largedel_cutlen, largeins_cutlen, full_window=True)

            mutation_side = [[],[]]

            for mut in mutations:
                if mut[0] in ['substitution']:
                    if window_front_st < mut[1] < window_front_ed:
                        mutation_side[0].append(mut)
                        front_mut.append(mut)
                    elif window_back_st < mut[1] < window_back_ed:
                        mutation_side[1].append(mut)
                        back_mut.append(mut)
                if mut[0] == 'deletion':
                    if window_front_st < mut[1] and window_front_ed > mut[1] + mut[2]:
                        mutation_side[0].append(mut)
                        front_mut.append(mut)
                    elif window_back_st < mut[1] and window_back_ed > mut[1] + mut[2]:
                        mutation_side[1].append(mut)
                        back_mut.append(mut)

            induced_mut_front[read_id] = sorted(mutation_side[0], key= lambda x: x[1])
            induced_mut_back[read_id] = sorted(mutation_side[1], key= lambda x: x[1], reverse=True)


            mutation_category_dict[read_id] = []

            additional_analysis_cateogry[read_id] = {}
            for i in ['WT', 'Del', 'Ins', 'Sub', 'LargeDel', 'LargeIns', 'Inv', 'Complex']:
                additional_analysis_cateogry[read_id][i] = 0

        marker_front = marker_find(induced_mut_front)
        marker_back = marker_find(induced_mut_back)

        if marker_front == False and marker_back == False:
            print('ERROR:: can not find proper mutation markers')
            sys.exit()
        
        n = 0
        confirmed_table = False
        read_ids = induced_mut_back.keys()

        if marker_front:
            for read_id, muts in marker_front.items():
                for mut in muts:
                    if mut[0] == 'deletion':
                        mutation_category_dict[read_id].append(['DEL', [mut[1], mut[1] + mut[2], 0.8], True, True, False])
                    elif mut[0] == 'substitution':
                        mutation_category_dict[read_id].append(['SUB', [mut[1], 0, 1, mut[3], mut[4]], True, True, False])
        else:
            marker_front_mild = marker_find_mild(induced_mut_front)
            if marker_front_mild:
                for read_id, muts in marker_front_mild.items():
                    for mut in muts:
                        if mut[0] == 'deletion':
                            mutation_category_dict[read_id].append(['DEL', [mut[1], mut[1] + mut[2], 0.8], True, True, False])
                        elif mut[0] == 'substitution':
                            mutation_category_dict[read_id].append(['SUB', [mut[1], 0, 1, mut[3], mut[4]], True, True, False])
            

        if marker_back:
            for read_id, muts in marker_back.items():
                for mut in muts:
                    if mut[0] == 'deletion':
                        mutation_category_dict[read_id].append(['DEL', [mut[1], mut[1] + mut[2], 0.8], True, True, False])
                    elif mut[0] == 'substitution':
                        mutation_category_dict[read_id].append(['SUB', [mut[1], 0, 1, mut[3], mut[4]], True, True, False]) 
        else:   
            marker_back_mild = marker_find_mild(induced_mut_back)
            if marker_back_mild:
                for read_id, muts in marker_back_mild.items():
                    for mut in muts:
                        if mut[0] == 'deletion':
                            mutation_category_dict[read_id].append(['DEL', [mut[1], mut[1] + mut[2], 0.8], True, True, False])
                        elif mut[0] == 'substitution':
                            mutation_category_dict[read_id].append(['SUB', [mut[1], 0, 1, mut[3], mut[4]], True, True, False])

        updated_mutation_category_dict = {}
        for read_id in read_ids:
            updated_mutation_category_dict[read_id] = []
        
        for read_id, muts in mutation_category_dict.items():
            for mut in muts:
                updated_mutation_category_dict[read_id].append(mut)
                if mut[0] == 'SUB':
                    for read_id_2, muts_2 in updated_mutation_category_dict.items():
                        if read_id == read_id_2: continue
                        range_check = True
                        for mut_2 in mutation_category_dict[read_id_2]:
                            if mut_2[0] in ['SUB', 'WT'] and mut[1][0] == mut_2[1][0]:
                                range_check = False
                                break
                            elif mut_2[0] == 'DEL' and mut_2[1][0] <= mut[1][1] <= mut_2[1][1]:
                                range_check = False
                                break

                        for mut_2 in muts_2:
                            if mut_2[0] in ['SUB', 'WT'] and mut[1][0] == mut_2[1][0]:
                                range_check = False
                                break
                            elif mut_2[0] == 'DEL' and mut_2[1][0] <= mut[1][1] <= mut_2[1][1]:
                                range_check = False
                                break
                        if range_check:
                            updated_mutation_category_dict[read_id_2].append(['WT', [mut[1][0], mut[1][0], 1], True, True, False])

                if mut[0] == 'DEL':
                    for read_id_2, muts_2 in updated_mutation_category_dict.items():
                        if read_id == read_id_2: continue
                        range_check = True
                        for mut_2 in mutation_category_dict[read_id_2]:
                            if mut_2[0] in ['SUB', 'WT'] and mut[1][0] <= mut_2[1][0] <= mut[1][1]:
                                range_check = False
                                break
                            elif mut_2[0] == 'DEL' and (mut_2[1][0] > mut[1][1]  or mut[1][0] > mut_2[1][1]):
                                range_check = False
                                break

                        for mut_2 in muts_2:
                            if mut_2[0] in ['SUB', 'WT'] and mut[1][0] <= mut_2[1][0] <= mut[1][1]:
                                range_check = False
                                break
                            elif mut_2[0] == 'DEL' and (mut_2[1][0] > mut[1][1]  or mut[1][0] > mut_2[1][1]):
                                range_check = False
                                break
                        if range_check:
                            updated_mutation_category_dict[read_id_2].append(['WT', [mut[1][0], mut[1][1], 0.8], True, True, False])



        print(updated_mutation_category_dict)
        mutation_category = []
        for read_id, mut in updated_mutation_category_dict.items():
            mutation_category.append([read_id, mut])
    
        return mutation_category, additional_analysis_cateogry
                
    

    parser = argparse.ArgumentParser(description='Re-analyze CRISPRlungo results with Customized Mutation Category')

    parser.add_argument('analysis_file_dir', type=str,  help='directory for CRISPRlungo ouptut')

    parser.add_argument('-i', '--custom_category_file', type=str,  help='file for custom group conditions or allele reference fasta file')
    parser.add_argument('-r', '--reference_fasta', type=str, help='FASTA file for allele reference. The mutation should be located out of the window range.')

    parser.add_argument('--mut_range_reference', type=int, default=1000, help='The range from both end of reference. The mutations only in this range are used for allele separation')
    parser.add_argument('--mut_num_reference', type=str, default='Auto', help='The range from both end of reference. The mutations only in this range are used for allele separation')

    parser.add_argument('--min_read_cnt', type=int, default=0, help='After counting based on mutation pattern, reads with counts less than the value are removed')
    parser.add_argument('--min_read_freq', type=float, default=0, help='After counting based on mutation pattern, reads with frequency less than the value are removed')
    parser.add_argument('--allele_plot_window', type=int, default=20, help='Window for allele plot, [default: window + 10]')
    parser.add_argument('--show_all_between_allele', action='store_true', help='Draw all sequences between two targets in an allele plot')

    args = parser.parse_args()

    if args.custom_category_file == None and args.reference_fasta == None:
        print('ERROR: Please input custom_category_file (-i) or reference_fasta (-r)')
        sys.exit()
    

    analysis_res_dir = args.analysis_file_dir
    create_dir(analysis_res_dir + '/custom_results/')

    f = open(analysis_res_dir + '/results/input_summary.txt').readlines()
    input_opt = {}
    for i in f:
        i = i.lower().strip().split(':')
        input_opt[i[0].replace(' ', '')] = ':'.join(i[1:]).replace(' ', '')

    ref_seq = input_opt['ref_seq'].upper()
    cv_pos = int(input_opt['cleavagepos_1'])
    strand = int(input_opt['cleavagestrand_1'])
    target = input_opt['target_1'].upper()
    cv_pos2 = input_opt['cleavagepos_2']
    window = int(input_opt['window'])
    window_between = input_opt['window_between_cleavage']
    largeins_cutlen = int(input_opt['largeins_cutlen'])
    largedel_cutlen = int(input_opt['largedel_cutlen'])

    if cv_pos2.upper() == 'FALSE':
        cv_pos2 = False
        strand_2 = False
        target_2 = False
    else:
        cv_pos2 = int(cv_pos2)
        strand_2 = input_opt['cleavagestrand_2']
        target_2 = input_opt['target_2'].upper()
    induced_mut = input_opt['induced_mut'].split(',')

    if induced_mut[0].upper() == 'FALSE':
        induced_mut = False
    cleavage_pos = int(input_opt['cleavagepos_1'])
    original_target = input_opt['target_1']


    if args.custom_category_file != None:
        mutation_category, additional_analysis_category = set_custom_mutation(args.custom_category_file, cv_pos, cv_pos2, induced_mut, analysis_res_dir)
    elif args.reference_fasta != None:
        mutation_category, additional_analysis_category = set_mutation_reference(args.reference_fasta, cv_pos, cv_pos2, induced_mut, analysis_res_dir, args.mut_range_reference)
    else:
        print('ERROR: input mutation information files (--custom_category_file or --reference_fasta)')
        sys.exit()
    for i in mutation_category: print(i)

    f = open(analysis_res_dir + '/results/read_classification.txt').readlines()

    classification_cnt = {'all': 0 }
    classification_fw = {}
    classification_read = {}
    for i in mutation_category:
        mut_name = i[0]
        if mut_name not in classification_cnt:
            classification_cnt[mut_name] = 0
            classification_read[mut_name] = {}
            create_dir(analysis_res_dir + '/custom_results/' + mut_name)
            classification_fw[mut_name] = open(analysis_res_dir + '/custom_results/' + mut_name + '/read_classification.txt', 'w')
            classification_fw[mut_name].write(f[0])
    classification_cnt['Others'] = 0
    classification_read['Others'] = {}
    create_dir(analysis_res_dir + '/custom_results/Others')
    classification_fw['Others'] = open(analysis_res_dir + '/custom_results/Others/read_classification.txt', 'w')
    classification_fw['Others'].write(f[0])
    total_cnt = 0

    for line in f[1:]:
        if line.strip() == '':
            continue
        line_sp = line.strip().split('\t')
        custom_mut_type = 'Others'
        filtered_mut_info = line_sp[2].split(',')
        additional_mut_info = line_sp[5].split(',')

        if additional_mut_info == ['[]'] or additional_mut_info == ['None']:
            additional_mut_info = []
        if filtered_mut_info == ['None']:
            filtered_mut_info = []

        whole_mut_info = filtered_mut_info + additional_mut_info

        if whole_mut_info == ['None']:
            whole_mut_info = []

        for category in mutation_category:
            check = True
            obtained_mut = []
            ##print(category[0])
            for condition in category[1]:
                if condition[2]:
                    used_mut_list = whole_mut_info
                else:
                    used_mut_list = filtered_mut_info
                True_or_Not = condition[3]
                cutoff = condition[1][2]
                if condition[0] == 'WT':
                    mut_len = 0
                    for m in used_mut_list:
                        if check_position(m, condition[1][0], condition[1][1]) == True_or_Not:
                            if m.find('Del') != -1:
                                del_st = int(m.split(':')[0].split('_')[0])
                                del_ed = int(m.split(':')[0].split('_')[1])
                                if del_st < condition[1][0]:
                                    del_st = condition[1][0]
                                if del_ed > condition[1][1]:
                                    del_ed = condition[1][1]
                                mut_len += del_ed - del_st + 1
                            elif m.find('Ins') != -1:
                                mut_len += 1
                            elif m.find('Sub') != -1:
                                mut_len += int(m.split('_')[2])
                    if mut_len >= (condition[1][1] - condition[1][0] + 1) * condition[1][2]:
                        check = False
                        break
                elif type(condition[0]) == list:
                    all_desired = len(condition[0])
                    counted_desired = 0
                    for m in used_mut_list:
                        if m in condition[0]:
                            counted_desired += 1
                            if m not in filtered_mut_info:
                                obtained_mut.append(m)
                    ##print(counted_desired, all_desired*cutoff)
                    if (counted_desired < all_desired * cutoff) == True_or_Not:
                        check = False
                        break
                elif condition[0] == 'DEL':
                    len_threshold = (condition[1][1] - condition[1][0] + 1) * cutoff
                    del_len = 0
                    for m in used_mut_list:
                        if m.find('Del') != -1 and check_position(m, condition[1][0], condition[1][1]):
                            del_st = int(m.split(':')[0].split('_')[0])
                            del_ed = int(m.split(':')[0].split('_')[1])
                            if del_st < condition[1][0]:
                                del_st = condition[1][0]
                            if del_ed > condition[1][1]:
                                del_ed = condition[1][1]
                            del_len += del_ed - del_st + 1
                            if m not in filtered_mut_info:
                                obtained_mut.append(m)
                    ##print(condition[1][0], condition[1][1])
                    ##print('del_len', del_len, len_threshold)
                    if (del_len < len_threshold) == True_or_Not:
                        check = False
                        break
                elif condition[0] == 'SUB':
                    sub_pos = condition[1][0]
                    sub_pat = condition[1][4]
                    get_sub = False

                    for m in used_mut_list:
                        if m.find('Sub') != -1:
                            sub_len = int(m.split('_')[2])
                            mut_nt = m.split('>')[1]
                            mut_st = int(m.split('_')[0])
                            for x in range(sub_len):
                                if mut_st + x == sub_pos and mut_nt[x] == sub_pat:
                                    get_sub = True
                                    break
                                if m not in filtered_mut_info:
                                    obtained_mut.append(m)
                            if get_sub == True:
                                break
                    if get_sub != True_or_Not:
                        check = False
                        break
            
            if check:
                custom_mut_type = category[0]
                if condition[2] and condition[4]:
                    filtered_mut_info += obtained_mut
                    filtered_mut_info = sorted(filtered_mut_info, key= lambda x: x.split('_')[0])
                break
            
                
        ##print(line)
        ##print(custom_mut_type)
        ##input()


        if custom_mut_type in additional_analysis_category:
            mut_classification = line_sp[1].replace('WithLargeMut', '')
            additional_analysis_category[custom_mut_type][mut_classification] += 1
            if line_sp[4] == '-':
                line_sp[4] = mut_classification
            else:
                line_sp[4] = mut_classification + '_' + line_sp[4]
            final_mut_info = f'{line_sp[3]}\t{line_sp[4]}\t' + ','.join(filtered_mut_info)
        else:
            final_mut_info = f'{line_sp[3]}\t{line_sp[4]}\t' + ','.join(filtered_mut_info)
        classification_cnt[custom_mut_type] += 1
        if final_mut_info not in classification_read[custom_mut_type]:
            classification_read[custom_mut_type][final_mut_info] = 1
        else:
            classification_read[custom_mut_type][final_mut_info] += 1
        total_cnt += 1
        classification_fw[custom_mut_type].write(line)

    for fw in classification_fw.values():
        fw.close()

    filted_classification_cnt = {}
    filted_total_cnt = 0
    for annot, dict in classification_read.items():
        filted_classification_cnt[annot] = {}
        for x, cnt in dict.items():
            if cnt < args.min_read_cnt or cnt < total_cnt * (args.min_read_freq / 100):
                continue
            filted_classification_cnt[annot][x] = cnt
            filted_total_cnt += cnt

    print('Total cnt', total_cnt)
    print('Filted cnt', 
    filted_total_cnt)

    fw = open(analysis_res_dir + '/custom_results/custom_mutation_patter_count.txt', 'w')

    for anno, res_dict in filted_classification_cnt.items():
        group_all_cnt = sum(res_dict.values())
        for mut_pat, cnt in sorted(res_dict.items(), key=lambda x: x[1], reverse=True):
            per = round(cnt*100/group_all_cnt, 2)
            fw.write(f'{anno}\t{cnt}\t{per}\t{mut_pat}\n')
    fw.close()

    fw = open(analysis_res_dir + '/custom_results/custom_mutation_allele_plot_input.txt', 'w')
    fw.write('Pat\tcount\t(%)\t\t\n')
    allele_line= 3
    plot_line = 0
    
    for anno, res_dict in filted_classification_cnt.items():
        n = 0
        group_all_cnt = sum(res_dict.values())
        plot_line += 1
        for mut_pat, cnt in sorted(res_dict.items(), key = lambda x: x[1], reverse=True):
            n += 1
            if anno == 'WT' or mut_pat.split('\t')[2] == '':
                mut_pat = '\t'.join(mut_pat.split('\t')[:2] + ['None'])
            per = round(cnt*100/group_all_cnt, 2)
            fw.write(f'{anno}\t{cnt}\t{per}\t{mut_pat}\n')
            plot_line += 1
            if n == allele_line:
                break


    fw.close()

    print(classification_cnt)
    
    print('Drawing graphs ...\r', end='')
    if induced_mut:
        induced_mut = ','.join(induced_mut)
    visual.custom_mutation_pie_chart(classification_cnt, analysis_res_dir + '/custom_results')

    visual.allele_plot(ref_seq, cv_pos, cv_pos2, strand, strand_2, analysis_res_dir + '/custom_results/custom_mutation_allele_plot_input.txt', analysis_res_dir + '/custom_results/', int(input_opt['cut_pos_in_target']), target, target_2, 1, args.min_read_cnt, args.min_read_freq, args.allele_plot_window, plot_line, induced_mut, args.show_all_between_allele, group_separate=True)
    min_read_cnt = 0
    min_read_freq = 0
    induced_sequence_path = False
    ref_seq = ''
    for i in open(analysis_res_dir + '/results/input_summary.txt').readlines():
        i = i.strip().split(':')
        x = i[0].strip()
        if x == 'minimum_read_count':
            min_read_cnt = int(i[1])
        elif x == 'minimum_read_frequency':
            min_read_freq = int(i[1])
        elif x == 'induced_sequence_path':
            induced_sequence_path = i[1]
        elif x == 'Ref_seq':
            ref_seq = i[1]

    for i in mutation_category:
        plots = {}
        tsv_file = f'{analysis_res_dir}/custom_results/{i[0]}/read_classification.txt'
        read_cnt_file = f'{analysis_res_dir}/custom_results/{i[0]}/mutation_patter_count.txt'
        graph_output_dir = f'{analysis_res_dir}/custom_results/{i[0]}'
        visual.write_read_count(tsv_file,  f'{analysis_res_dir}/results/preprocess_count.txt', read_cnt_file, f'{analysis_res_dir}/custom_results/{i[0]}/mutation_summary_count.txt', min_read_cnt, min_read_freq, induced_sequence_path, out_print=False)
        read_per_position = visual.visualization_preprocess_regular(analysis_res_dir + '/align/Treated_alignment.sam', analysis_res_dir + '/ref_seq/ref_wo_umi.fasta')
        plots['mutation_pie'], plots['pattern_pie'], plots['allele_pie'] = visual.mutation_pie_chart(read_cnt_file, graph_output_dir)
        plots['indel_per_pos'] = visual.indel_per_position(read_cnt_file, ref_seq, graph_output_dir)

if __name__=='__main__':
	main()



        


                    
                    



