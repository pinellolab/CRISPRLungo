import sys, argparse
import CRISPRlungo_visualization as visual


def run_customized_analysis():

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
        if int(mut[1]) <= start_pos or int(mut[0]) >= end_pos:
            return False
        return True

    def set_custom_mutation(input_file, cv_pos, cv_pos2, induced_mut):
        mutation_category = []
        additional_analysis_cateogry = {}

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
                    for i in ['WT', 'Del', 'Ins', 'Sub', 'LargeDel', 'LargeIns', 'Inv']:
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
                    mutation_category[-1][1].append(['SUB', [start_pos, 0, line[4]], mutation_marker, OrNot_marker, Include_marker])
            

        return mutation_category, additional_analysis_cateogry
    

    parser = argparse.ArgumentParser(description='Re-analyze CRISPRlungo results with Customized Mutation Category')

    parser.add_argument('analysis_file_dir', type=str,  help='directory for CRISPRlungo ouptut')
    parser.add_argument('custom_category_file', type=str,  help='directory for CRISPRlungo ouptut')

    parser.add_argument('--min_read_cnt', type=int, default=0, help='After counting based on mutation pattern, reads with counts less than the value are removed')
    parser.add_argument('--min_read_freq', type=float, default=0, help='After counting based on mutation pattern, reads with frequency less than the value are removed')
    parser.add_argument('--allele_plot_window', type=int, default=20, help='Window for allele plot, [default: window + 10]')
    parser.add_argument('--show_all_between_allele', action='store_true', help='Draw all sequences between two targets in an allele plot')

    args = parser.parse_args()

    analysis_res_dir = args.analysis_file_dir

    f = open(analysis_res_dir + '/results/input_summary.txt').readlines()

    ref_seq = f[2].split(':')[1]
    cv_pos = int(f[3].strip().split(':')[1])
    strand = f[4].strip().split(':')[1]
    target = f[0].strip().split(':')[1]
    cv_pos2 = f[5].strip().split(':')[1]
    if cv_pos2 == 'False':
        cv_pos2 = False
        strand_2 = False
        target_2 = False
    else:
        cv_pos2 = int(cv_pos2)
        strand_2 = f[6].strip().split(':')[1]
        target_2 = f[1].strip().split(':')[1]
    induced_mut = f[7].strip().split(' :')[1].split(',')
    if induced_mut[0] == 'False':
        induced_mut = False
    cleavage_pos = int(f[8].strip().split(' :')[1])
    original_target = int(f[9].strip().split(' :')[1])
        
    mutation_category, additional_analysis_category = set_custom_mutation(args.custom_category_file, cv_pos, cv_pos2, induced_mut)

    f = open(analysis_res_dir + '/results/read_classification.txt').readlines()

    classification_cnt = {'all': 0 }
    classification_read = {}
    for i in mutation_category:
        mut_name = i[0]
        if mut_name not in classification_cnt:
            classification_cnt[mut_name] = 0
            classification_read[mut_name] = {}
    classification_cnt['Others'] = 0
    classification_read['Others'] = {}
    total_cnt = 0

    for line in f[1:]:
        if line.strip() == '':
            continue
        line_sp = line.strip().split('\t')
        custom_mut_type = 'Others'
        filtered_mut_info = line_sp[2].split(',')

        if filtered_mut_info == ['None']:
            filtered_mut_info = []
        whole_mut_info = filtered_mut_info + line_sp[5].split(',')

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
                    if mut_len > (condition[1][1] - condition[1][0] + 1) * condition[1][2]:
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
                    sub_pat = condition[1][2]
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
            mut_classification = line_sp[1]
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
    print('Filted cnt', filted_total_cnt)

    fw = open(analysis_res_dir + '/custom_results/custom_mutation_patter_count.txt', 'w')

    for anno, res_dict in filted_classification_cnt.items():
        for mut_pat, cnt in sorted(res_dict.items(), key=lambda x: x[1], reverse=True):
            per = round(cnt*100/filted_total_cnt, 2)
            fw.write(f'{anno}\t{cnt}\t{per}\t{mut_pat}\n')
    fw.close()

    fw = open(analysis_res_dir + '/custom_results/custom_mutation_allele_plot_input.txt', 'w')
    fw.write('Pat\tcount\t(%)\t\t\n')
    allele_line= 10
    for anno, res_dict in filted_classification_cnt.items():
        n = 0
        for mut_pat, cnt in sorted(res_dict.items(), key = lambda x: x[1], reverse=True):
            n += 1
            if anno == 'WT' or mut_pat.split('\t')[2] == '':
                mut_pat = '\t'.join(mut_pat.split('\t')[:2] + ['None'])
            per = round(cnt*100/filted_total_cnt, 2)
            fw.write(f'{anno}\t{cnt}\t{per}\t{mut_pat}\n')
            if n == allele_line:
                break
    fw.close()

    print(classification_cnt)
    
    if induced_mut:
        induced_mut = ','.join(induced_mut)
    visual.custom_mutation_pie_chart(classification_cnt, analysis_res_dir + '/custom_results')
    visual.allele_plot(ref_seq, cv_pos, cv_pos2, strand, strand_2, analysis_res_dir + '/custom_results/custom_mutation_allele_plot_input.txt', analysis_res_dir + '/custom_results/', cleavage_pos, target, target_2, original_target, args.min_read_cnt, args.min_read_freq, args.allele_plot_window, allele_line*len(classification_read), induced_mut, args.show_all_between_allele, group_separate=True)

if __name__=='__main__':
	run_customized_analysis()



        


                    
                    



