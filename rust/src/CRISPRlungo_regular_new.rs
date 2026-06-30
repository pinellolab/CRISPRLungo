use fishers_exact::fishers_exact;
use statrs::distribution::{ChiSquared, ContinuousCDF};
use std::collections::HashSet;
use std::{collections::HashMap, path::absolute};

pub struct AnalysisResult {
    pub read_id: String,
    pub classification: String,
    pub mutation_info: String,
    pub precise_info: String,
}

#[derive(Debug, Eq, PartialEq, Hash, Clone)]
//0: Substitutions, 1: Insertions, 2: Deletions
pub struct MutationKeys {
    pub mutation: u8,
    pub start_ref: u32,
    pub end_ref: u32,
    pub mutation_len: u32,
    pub reference_seq: String,
    pub mutated_seq: String,
}

fn reverse_complementary(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|base| match base {
            'A' => 'T',
            'T' => 'A',
            'G' => 'C',
            'C' => 'G',
            _ => base,
        })
        .collect()
}

fn cigar_aligned_length(cigar: &str) -> u32 {
    let mut aligned_length = 0;
    let cigar_chars = cigar.chars();
    let mut sub_len = 0;

    for c in cigar_chars {
        if c.is_digit(10) {
            sub_len = sub_len * 10 + c.to_digit(10).unwrap() as i32;
        } else {
            if c == 'M' || c == 'D' {
                aligned_length += sub_len;
            }
            sub_len = 0;
        }
    }
    aligned_length as u32
}

fn check_in_window(
    mutation: String,
    window_filter: bool,
    cv_pos: i32,
    cv_pos_2: Option<i32>,
    window: i32,
    check_window_between_targets: bool,
) -> bool {
    if !window_filter {
        return true;
    }

    let window_st = cv_pos - window;
    let window_ed = cv_pos + window;

    let replaced = mutation.replace(":", "_");

    let mut_split: Vec<&str> = replaced.split("_").collect();
    let mut_st = mut_split[0].parse::<i32>().expect("Failed to parse mut_st");
    let mut mut_ed;
    if mutation.contains("Del") {
        mut_ed = mut_split[1].parse::<i32>().expect("Failed to parse mut_st") + 1;
    } else {
        mut_ed = mut_split[1].parse::<i32>().expect("Failed to parse mut_st");
    }

    let mut check_position_in_window = false;

    if !(mut_ed < window_st || mut_st > window_ed) {
        check_position_in_window = true;
    }

    if let Some(val) = cv_pos_2 {
        let window_st2 = if check_window_between_targets {
            cv_pos - window
        } else {
            val - window
        };
        let window_ed2 = val + window;

        if !(mut_ed < window_st2 || mut_st > window_ed2) {
            check_position_in_window = true;
        }
    }

    check_position_in_window
}

fn cigar_to_tup(cigar: &str) -> Vec<(u8, u32)> {
    let mut cigar_tup: Vec<(u8, u32)> = Vec::new();
    let cigar_chars = cigar.chars();
    let mut sub_len: u32 = 0;

    for c in cigar_chars {
        if c.is_digit(10) {
            sub_len = sub_len * 10 + c.to_digit(10).unwrap() as u32;
        } else {
            if c == 'M' {
                cigar_tup.push((0, sub_len));
            } else if c == 'I' {
                cigar_tup.push((1, sub_len));
            } else if c == 'D' {
                cigar_tup.push((2, sub_len));
            } else if c == 'S' {
                cigar_tup.push((4, sub_len));
            } else if c == 'H' {
                cigar_tup.push((5, sub_len));
            }
            sub_len = 0;
        }
    }
    cigar_tup
}

fn cigar_mutation_keys(
    ref_pos: u32,
    cigar_tup: Vec<(u8, u32)>,
    query_seq: &str,
    ref_char_vec: Vec<char>,
    overlap_len: u32,
    ref_start_pos: u32,
) -> Vec<String> {
    let mut mutations_in_read: Vec<String> = Vec::new();
    let mut seq_pos: u32 = 0;
    let mut ref_pos = ref_pos;

    let using_query_seq_vec: Vec<char> = query_seq.chars().collect();

    for (op, len) in cigar_tup {
        if op == 0 {
            let mut tmp_sub_vec: Vec<(u8, u32, u32, String, String)> = Vec::new();
            for x in 0..len {
                let index: usize = (ref_pos + x) as usize;
                let query_index: usize = (seq_pos + x) as usize;
                if using_query_seq_vec[(seq_pos + x) as usize] != ref_char_vec[index]
                    && ref_start_pos + overlap_len <= ref_pos + x
                {
                    let formated = format!(
                        "{}_{}:Sub_{}>{}",
                        ref_pos + x,
                        ref_pos + x,
                        ref_char_vec[index].to_string(),
                        using_query_seq_vec[query_index].to_string()
                    );
                    mutations_in_read.push(formated);
                }
            }
            /*if let Some((u_x, u_y, u_z, str_x, str_y))= tmp_sub_vec.last_mut() {
                        if *u_y + *u_z == ref_pos + x {
                            *u_z += 1;
                            *str_x += &(ref_char_vec[index].to_string());
                            *str_y += &(using_query_seq_vec[query_index].to_string());
                        } else {
                            tmp_sub_vec.push((0, ref_pos+x, 1, ref_char_vec[index].to_string(), using_query_seq_vec[query_index].to_string()));
                        }
                    } else {
                        tmp_sub_vec.push((0, ref_pos+x, 1, ref_char_vec[index].to_string(), using_query_seq_vec[query_index].to_string()));
                    }
                }
            }
            for tmp_sub in tmp_sub_vec {

                let formated = format!(
                    "{}_{}:Sub_{}>{}",
                    tmp_sub.1,
                    tmp_sub.1 + tmp_sub.2 - 1,
                    tmp_sub.3.to_string(),
                    tmp_sub.4.to_string()
                );

                mutations_in_read.push( formated );
            }*/
            ref_pos += len;
            seq_pos += len;
        } else if op == 1 {
            // R1 fix: only the recording is gated by the overlap guard; the
            // query position must ALWAYS advance for an insertion, otherwise
            // every subsequent coordinate in the read desyncs.
            if ref_start_pos + overlap_len < ref_pos {
                let query_index = seq_pos as usize;
                let formated = format!(
                    "{}_{}:Ins_{}",
                    ref_pos,
                    ref_pos + 1,
                    using_query_seq_vec[query_index..query_index + len as usize]
                        .iter()
                        .copied()
                        .collect::<String>()
                );
                mutations_in_read.push(formated);
            }
            seq_pos += len;
        } else if op == 2 {
            // R1 fix: reference position must ALWAYS advance for a deletion.
            if ref_start_pos + overlap_len < ref_pos {
                let formated = format!("{}_{}:Del_{}", ref_pos, ref_pos + len - 1, len);
                mutations_in_read.push(formated);
            }
            ref_pos += len;
        } else if op == 3 {
            // R1 fix: N (skipped reference region) consumes reference only.
            ref_pos += len;
        } else if op == 4 {
            seq_pos += len;
        }
    }
    mutations_in_read
}

fn analyze_sa_reads(
    partial_reads: Vec<(u32, u32, String, String, u16)>,
    query_seq: String,
    ref_char_vec: Vec<char>,
    ref_len: u32,
    query_name: &str,
) -> Vec<Vec<String>> {
    let mut query_char_vector: Vec<char> = query_seq.chars().collect();
    let mut partial_read_info: Vec<(u32, u32, u32, u32, i8, i8, String, String)> = Vec::new();

    for read in partial_reads {
        let (pos_refst, pos_refed, partial_seq, cigar, flag) = read;
        let mut cigar_tup = cigar_to_tup(&cigar);
        let mut pos_in_ref = false;

        let partial_seq_vec: Vec<char> = partial_seq.chars().collect();
        let query_seq_len: u32 = query_seq.len() as u32;
        let mut align_strand_in_refseq: i8 = 2;
        let mut align_strand_to_seq: i8 = 2;
        let mut pos_in_seq: u32 = 0;
        let mut pos_in_seq_ed: u32 = 0;
        let mut pos_find_check = false;

        let align_len = cigar_aligned_length(&cigar);

        if let Some(index) = query_seq.find(&partial_seq) {
            if cigar_tup[0].0 == 4 || cigar_tup[0].0 == 5 {
                pos_in_seq = cigar_tup[0].1;
                pos_find_check = true;
            } else {
                pos_in_seq = 0;
                pos_find_check = true;
            }
            if let Some((op, val)) = cigar_tup.last() {
                if *op == 4 || *op == 5 {
                    pos_in_seq_ed = query_seq_len - val - 1;
                } else {
                    pos_in_seq_ed = query_seq_len - 1;
                }
            } else {
                println!("Warning: truncated cigar");
            }
            align_strand_to_seq = 1;
        } else if let Some(index) = query_seq.find(&reverse_complementary(&partial_seq)) {
            if cigar_tup[0].0 == 4 || cigar_tup[0].0 == 5 {
                pos_in_seq_ed = query_seq_len - cigar_tup[0].1 - 1;
            } else {
                pos_in_seq_ed = query_seq_len - 1;
            }
            if let Some((op, val)) = cigar_tup.last() {
                if *op == 4 || *op == 5 {
                    pos_in_seq = *val;
                    pos_find_check = true;
                } else {
                    pos_in_seq = 0;
                    pos_find_check = true;
                }
            } else {
                println!("Warning: truncated cigar");
            }
            align_strand_to_seq = -1;
        }

        if pos_find_check == false {
            continue;
        }

        if pos_refst < 100 || pos_refed > ref_len - 100 {
            //pass
        } else {
            continue;
        }

        if pos_refst < 100 && pos_refed > ref_len - 100 {
            align_strand_in_refseq = 3;
        } else if flag & 0x10 != 0 {
            align_strand_in_refseq = 1;
        } else {
            align_strand_in_refseq = -1;
        }

        partial_read_info.push((
            pos_refst,
            pos_refed,
            pos_in_seq,
            pos_in_seq_ed,
            align_strand_in_refseq,
            align_strand_to_seq,
            cigar,
            partial_seq,
        ));
    }

    partial_read_info.sort_by(|a, b| a.2.cmp(&b.2));
    partial_read_info.push((2, 2, 2, 2, 2, 2, "".to_string(), "".to_string()));

    let mut mutations_in_reads: Vec<Vec<String>> = Vec::new();

    for i in 0..partial_read_info.len() - 1 {
        let mut mutations_in_read: Vec<String> = Vec::new();
        let mut info = &partial_read_info[i..i + 2];

        if info[1].4 == 3 {
            continue;
        }

        if info[0].4 == 3 {
            let mut ref_pos = info[0].0;
            let mut seq_pos = info[0].2;
            let using_query_seq = info[0].7.clone();

            //if info[0].5 == -1{
            //    let using_query_seq = reverse_complementary(&query_seq);
            //    seq_pos = (query_seq.len() as u32 - info[0].3) as u32;
            //}

            mutations_in_reads.push(cigar_mutation_keys(
                ref_pos,
                cigar_to_tup(&info[0].6),
                &using_query_seq,
                ref_char_vec.clone(),
                0,
                0,
            ));
            continue;
        }

        if info[0].5 == info[1].5
            && info[0].4 == info[1].4
            && info[0].5 == 1
            && info[0].1 < info[1].1
        {
            let mut overlap_len = 0;
            if info[1].0 < info[0].1 {
                overlap_len = info[0].1 - info[1].0;
            }

            if info[0].1 < info[1].0 - 100 {
                let formated = format!(
                    "{}_{}:Del_{}",
                    info[0].1,
                    info[1].0 - 1,
                    info[1].0 - info[0].1,
                );
                mutations_in_read.push(formated);
            }
            let mut_len = (info[0].3 as i32 - info[1].2 as i32).abs() as u32;
            if mut_len - 1 > 20 {
                let mut ins_end_pos = info[1].0 + overlap_len;
                if ins_end_pos <= info[0].1 {
                    ins_end_pos = info[0].1 + 1;
                }
                let mut ins_seq = "".to_string();
                if query_seq.contains(&info[0].7) {
                    let start_idx = (info[0].3 as usize) + 1;
                    let diff = (info[0].3 as i32 - info[1].2 as i32).abs() as usize;
                    let overlap = overlap_len as usize;

                    let end_idx = start_idx + diff - 1 + overlap;

                    ins_seq = query_seq[start_idx..end_idx].to_string();
                } else {
                    let pos0 = info[0].3 as i32;
                    let pos1 = info[1].2 as i32;
                    let overlap = overlap_len as usize;

                    let start_idx = (info[0].3 as usize + 1).saturating_sub(overlap);
                    let diff = (pos0 - pos1).abs() as usize;
                    let end_idx = (info[0].3 as usize) + 1 + diff - 1;

                    ins_seq = reverse_complementary(&query_seq[start_idx..end_idx]);
                }

                let formated = format!("{}_{}:Ins_{}", info[0].1, ins_end_pos, ins_seq);
                mutations_in_read.push(formated);
            }

            let mut info_n = 0;

            for sub_info in info {
                let mut overlap_st = 0;
                if info_n == 1 {
                    overlap_st = overlap_len;
                }
                info_n += 1;
                let ref_pos: u32 = sub_info.0;
                let seq_pos: u32 = 0;
                for tmp_sub in cigar_mutation_keys(
                    ref_pos,
                    cigar_to_tup(&sub_info.6),
                    &sub_info.7,
                    ref_char_vec.clone(),
                    overlap_len,
                    sub_info.0,
                ) {
                    mutations_in_read.push(tmp_sub);
                }
            }
            mutations_in_reads.push(mutations_in_read);
        } else if info[0].5 == info[1].5
            && info[0].4 == info[1].4
            && info[0].5 == -1
            && info[1].1 < info[0].1
        {
            let mut overlap_len = 0;
            if info[1].1 > info[0].0 {
                overlap_len = info[1].1 - info[0].0;
            }
            if info[1].1 < info[0].0 - 100 {
                let formated = format!(
                    "{}_{}:Del_{}",
                    info[1].1,
                    info[1].1 - 1,
                    info[0].0 - info[1].1 - 1,
                );
                mutations_in_read.push(formated);
            }
            let mut mut_len = (info[1].2 as i32 - info[0].3 as i32).abs() as u32;
            if (info[1].2 as i32 - info[0].3 as i32).abs() - 1 > 20 {
                let mut ins_end_pos = info[0].0 + overlap_len;
                if ins_end_pos <= info[1].1 {
                    ins_end_pos = info[1].1 + 1;
                }
                let mut ins_seq = "".to_string();
                let pos0 = info[0].3 as i32;
                let pos1 = info[1].2 as i32;
                let overlap = overlap_len as i32;
                let diff = (pos1 - pos0).abs();

                let start_idx = ((pos0 + 1 - overlap).max(0)) as usize;
                let end_idx = (pos0 + 1 + diff - 1 + overlap) as usize;

                if query_seq.contains(&info[0].7) {
                    ins_seq = query_char_vector[start_idx..end_idx]
                        .iter()
                        .collect::<String>();
                } else {
                    let temp_seq = query_char_vector[start_idx..end_idx]
                        .iter()
                        .collect::<String>();
                    ins_seq = reverse_complementary(&temp_seq);
                }

                let formated = format!("{}_{}:Ins_{}", info[1].1, ins_end_pos, ins_seq);
                mutations_in_read.push(formated);
            }
            let query_seq_len = query_seq.len() as u32;

            let mut info_n = 0;
            for sub_info in info {
                let ref_pos = sub_info.0;
                let mut seq_pos = 0;
                let mut overlap_st = 0;
                if info_n == 1 {
                    overlap_st = overlap_len;
                }
                info_n += 1;
                let used_query_seq = sub_info.7.clone();
                for tmp_sub in cigar_mutation_keys(
                    ref_pos,
                    cigar_to_tup(&sub_info.6),
                    &used_query_seq,
                    ref_char_vec.clone(),
                    overlap_st,
                    sub_info.0,
                ) {
                    mutations_in_read.push(tmp_sub);
                }
            }
            mutations_in_reads.push(mutations_in_read)
        }
    }
    mutations_in_reads
}

fn quant_unique_indels(
    align_res: Vec<String>,
    ref_seq: &str,
    induced_filter: bool,
    cv_pos: i32,
    cv_pos_2: Option<i32>,
    window: i32,
    check_window_between_targets: bool,
) -> (
    HashMap<String, i32>,
    HashMap<String, i32>,
    HashMap<String, Vec<Vec<String>>>,
) {
    let mut align_summary_cnt: HashMap<String, i32> = HashMap::new();
    let mut mutations: HashMap<String, i32> = HashMap::new();
    let mut align_read_to_mut: HashMap<String, Vec<Vec<String>>> = HashMap::new();

    let ref_len: u32 = ref_seq.len() as u32;
    let ref_char_vec: Vec<char> = ref_seq.chars().collect();

    align_summary_cnt.insert("All_reads".to_string(), 0);
    align_summary_cnt.insert("Unmapped".to_string(), 0);
    align_summary_cnt.insert("Low_quality".to_string(), 0);
    align_summary_cnt.insert("Short".to_string(), 0);
    align_summary_cnt.insert("Used".to_string(), 0);

    let mut align_iter = align_res.iter();

    while let Some(read) = align_iter.next() {
        let fields: Vec<&str> = read.trim().split('\t').collect();

        // check Falg
        if (fields.len() < 9) {
            continue;
        }

        let flag = match fields[1].parse::<u16>() {
            Ok(val) => val,
            Err(_) => { eprintln!("Skipping malformed line (flag)"); continue; }
        };

        if (flag & 0x100) != 0 || (flag & 0x800) != 0 {
            continue;
        }

        let mapq = match fields[4].parse::<i32>() {
            Ok(num) => num,
            Err(_) => { eprintln!("Skipping malformed line (MAPQ)"); continue; }
        };

        align_summary_cnt
            .get_mut("All_reads")
            .map(|count| *count += 1);

        if mapq <= 30 {
            if let Some(count) = align_summary_cnt.get_mut("Low_quality") {
                *count += 1;
            }
            continue;
        }

        let read_id: String = fields[0].to_string();

        let start_pos = match fields[3].parse::<u32>() {
            Ok(num) => num - 1,
            Err(_) => { eprintln!("Skipping malformed line (POS)"); continue; }
        };

        let cigar = fields[5].to_string();

        let tags: Vec<String> = if fields.len() > 11 {
            fields[11..].iter().map(|s| s.to_string()).collect()
        } else {
            Vec::new()
        };

        let mut sa_n: usize = 0;
        if let Some(sa_field) = fields.iter().find(|s| s.starts_with("SA:")) {
            sa_n = sa_field.matches(';').count() + 1;
        }
        //let sa_n: usize = tags.iter().filter(|tag| tag.starts_with("SA:")).count() + 1;

        let sa_check = sa_n > 0;
        let query_seq = fields[9].to_string();

        if sa_check {
            let mut sa_reads: Vec<(u32, u32, String, String, u16)> = Vec::new();
            sa_reads.push((
                start_pos,
                start_pos + cigar_aligned_length(&cigar),
                query_seq.clone(),
                cigar,
                flag,
            ));
            let mut next_cnt = 0;
            let query_name = fields[0].to_string();

            while sa_reads.len() < sa_n {
                if let Some(next_read) = align_iter.next() {
                    let next_fields: Vec<&str> = next_read.trim().split('\t').collect();
                    let next_flag = match next_fields[1].parse::<u16>() {
                        Ok(num) => num,
                        Err(_) => {
                            eprintln!("FLAG parsing error");
                            panic!("Stopping the software due to parsing failure");
                        }
                    };
                    if (next_flag & 0x800) != 0 {
                        let next_start_pos = match next_fields[3].parse::<u32>() {
                            Ok(num) => num - 1,
                            Err(_) => {
                                eprintln!("POS parsing error");
                                panic!("Stopping the software due to parsing failure");
                            }
                        };

                        let next_cigar = next_fields[5].to_string();
                        let next_query_seq = next_fields[9].to_string();
                        sa_reads.push((
                            next_start_pos,
                            next_start_pos + cigar_aligned_length(&next_cigar),
                            next_query_seq,
                            next_cigar,
                            next_flag,
                        ));
                    }
                } else {
                    break;
                }
            }

            let sa_reads_muts: Vec<Vec<String>> = analyze_sa_reads(
                sa_reads.clone(),
                query_seq.clone(),
                ref_char_vec.clone(),
                ref_len,
                &query_name,
            );
            let sa_reads_muts_len = sa_reads_muts.len();
            if sa_reads_muts_len == 0 {
                align_summary_cnt.get_mut("Short").map(|count| *count += 1);
                continue;
            }

            let mut n = -1;
            for sa_mutations in sa_reads_muts {
                n += 1;
                align_summary_cnt.get_mut("Used").map(|count| *count += 1);
                let mut mutation_keys_vec: Vec<String> = Vec::new();
                let mut out_mutation_keys_vec: Vec<String> = Vec::new();
                for tmp_mut in sa_mutations {
                    if (!check_in_window(
                        tmp_mut.clone(),
                        induced_filter,
                        cv_pos,
                        cv_pos_2,
                        window,
                        check_window_between_targets,
                    )) {
                        out_mutation_keys_vec.push(tmp_mut);
                        continue;
                    }
                    mutations
                        .entry(tmp_mut.clone())
                        .and_modify(|e| *e += 1)
                        .or_insert(1);
                    mutation_keys_vec.push(tmp_mut);
                }
                if sa_reads_muts_len > 1 {
                    align_read_to_mut.insert(
                        format!("{}_s{}", read_id, n).to_string(),
                        vec![mutation_keys_vec.clone(), out_mutation_keys_vec.clone()],
                    );
                } else {
                    align_read_to_mut.insert(
                        read_id.to_string(),
                        vec![mutation_keys_vec.clone(), out_mutation_keys_vec.clone()],
                    );
                }
            }
            // Run analyze_sa_reads func
        } else {
            let end_pos = start_pos + cigar_aligned_length(&cigar);
            if start_pos > 100 || end_pos < ref_len - 100 {
                align_summary_cnt.get_mut("Short").map(|count| *count += 1);
                continue;
            }
            let mut query_pos = 0;
            let mut mutation_keys_vec: Vec<String> = Vec::new();
            let mut out_mutation_keys_vec: Vec<String> = Vec::new();
            align_summary_cnt.get_mut("Used").map(|count| *count += 1);

            for tmp_mut in cigar_mutation_keys(
                start_pos,
                cigar_to_tup(&cigar),
                &query_seq,
                ref_char_vec.clone(),
                0,
                0,
            ) {
                if (!check_in_window(
                    tmp_mut.clone(),
                    induced_filter,
                    cv_pos,
                    cv_pos_2,
                    window,
                    check_window_between_targets,
                )) {
                    out_mutation_keys_vec.push(tmp_mut);
                    continue;
                }
                mutations
                    .entry(tmp_mut.clone())
                    .and_modify(|e| *e += 1)
                    .or_insert(1);
                mutation_keys_vec.push(tmp_mut);
            }
            align_read_to_mut.insert(
                read_id.to_string(),
                vec![mutation_keys_vec.clone(), out_mutation_keys_vec.clone()],
            );
        }
    }

    (align_summary_cnt, mutations, align_read_to_mut)
}

fn custom_round(num: u32, margin: f64, length: u32) -> u32 {
    let step = (length as f64 * margin).floor();
    if step == 0.0 {
        return num;
    }
    let quotient = num as f64 / step;
    let rounded = quotient.round();
    let closest_value = rounded * step;
    closest_value as u32
}

fn pool(mut mut_dict: HashMap<String, i32>, allowance: f64) -> HashMap<String, u32> {
    let mut outputdict: HashMap<String, u32> = HashMap::new();

    for (key, count) in mut_dict.drain() {
        let new_key = key.clone();
        let replaced = new_key.replace(":", "_"); // 임시값을 replaced에 저장
        let mut_split: Vec<&str> = replaced.split("_").collect();
        let mut mutation_len = 0;
        if (mut_split[2] == "Sub" || mut_split[2] == "Ins") {
            mutation_len = mut_split[3].len() as i32;
        } else {
            mutation_len = mut_split[3].parse::<i32>().expect("Failed to parse mut_st");
        }
        if (mutation_len as f64 * allowance).floor() == 0.0 {
            *outputdict.entry(new_key).or_insert(0) += count as u32;
        } else {
            let mut new_key = key.clone();
            *outputdict.entry(new_key).or_insert(0) += count as u32;
        }
    }
    outputdict
}

fn chi2_contingency(table: &[Vec<f64>]) -> f64 {
    let n_rows = table.len();
    let n_cols = table[0].len();

    let mut row_sums = vec![0.0; n_rows];
    let mut col_sums = vec![0.0; n_cols];
    let mut total = 0.0;
    for (i, row) in table.iter().enumerate() {
        for (j, &value) in row.iter().enumerate() {
            row_sums[i] += value;
            col_sums[j] += value;
            total += value;
        }
    }

    let mut chi2_stat = 0.0;
    for i in 0..n_rows {
        for j in 0..n_cols {
            let expected = row_sums[i] * col_sums[j] / total;
            if expected != 0.0 {
                let diff = (table[i][j] - expected).abs();
                let corrected_diff = if n_rows == 2 && n_cols == 2 {
                    (diff - 0.5).max(0.0)
                } else {
                    diff
                };
                chi2_stat += corrected_diff.powi(2) / expected;
            }
        }
    }

    let dof = ((n_rows - 1) * (n_cols - 1)) as u32;

    let chi2_dist = ChiSquared::new(dof as f64).expect("Can not make ChiSquared distribution");
    let p_value = 1.0 - chi2_dist.cdf(chi2_stat);

    p_value
}

fn calculate_pvalue(
    control_mut_cnt: u32,
    control_all_cnt: u32,
    treated_mut_cnt: u32,
    treated_all_cnt: u32,
    key: &String,
) -> f64 {
    let p_value: f64 = if control_mut_cnt as f64 / control_all_cnt as f64
        > treated_mut_cnt as f64 / treated_all_cnt as f64
    {
        2.0
    } else if control_mut_cnt == 0 && treated_mut_cnt != 0 {
        let mut mut_len: u32 = 0;
        if key.contains("Del") {
            mut_len = key.rsplit('_').next().unwrap().parse::<u32>().unwrap();
        } else if key.contains("Ins") {
            mut_len = key.rsplit('_').next().unwrap().len() as u32;
        } else if key.contains("Sub") {
            mut_len = key.rsplit('>').next().unwrap().len() as u32;
        }
        if mut_len < 5 && (treated_mut_cnt as f64 / treated_all_cnt as f64) < 0.001 {
            2.0
        } else {
            0.0
        }
    } else if control_mut_cnt <= 5 {
        fishers_exact(&[
            control_mut_cnt,
            control_all_cnt - control_mut_cnt,
            treated_mut_cnt,
            treated_all_cnt - treated_mut_cnt,
        ])
        .unwrap()
        .two_tail_pvalue
    } else {
        let table = vec![
            vec![
                control_mut_cnt as f64,
                (control_all_cnt - control_mut_cnt) as f64,
            ],
            vec![
                treated_mut_cnt as f64,
                (treated_all_cnt - treated_mut_cnt) as f64,
            ],
        ];
        chi2_contingency(&table)
    };
    p_value
}

fn significant_mutations(
    control1: HashMap<String, u32>,
    treated1: HashMap<String, u32>,
    control_used_cnt: i32,
    treated_used_cnt: i32,
    p_limit_value: f64,
    wo_control: bool,
    alpha: f64,
    length_min: i32,
) -> HashMap<String, (f64, u32)> {
    fn lengths_from_control(control_dict: &HashMap<String, u32>) -> (Vec<u32>, Vec<u32>) {
        let mut del_lens: Vec<u32> = Vec::new();
        let mut ins_lens: Vec<u32> = Vec::new();
        let mut mut_len: u32;

        for (key, value) in control_dict.iter() {
            if key.contains("Del") {
                mut_len = key.rsplit('_').next().unwrap().parse::<u32>().unwrap();
                del_lens.push(mut_len);
            } else if key.contains("Ins") {
                mut_len = key.rsplit('_').next().unwrap().len() as u32;
                ins_lens.push(mut_len);
            }
        }
        del_lens.sort();
        ins_lens.sort();
        (del_lens, ins_lens)
    }

    fn probability_based_on_length(
        mutation: &String,
        ins_sorted: &Vec<u32>,
        del_sorted: &Vec<u32>,
    ) -> (f64, i32) {
        let mut mut_len: u32 = 1;
        let mut arr: Vec<u32>;
        let mut pval: f64;

        if mutation.contains("Sub") {
            return (1.0, 1);
        }
        if mutation.contains("Del") {
            mut_len = mutation.rsplit('_').next().unwrap().parse::<u32>().unwrap();
            arr = del_sorted.clone();
        } else if mutation.contains("Ins") {
            mut_len = mutation.rsplit('_').next().unwrap().len() as u32;
            arr = ins_sorted.clone();
        } else {
            return (0.0, 1);
        }

        let idx = match arr.binary_search(&(mut_len - 1)) {
            Ok(pos) => pos + 1,
            Err(pos) => pos,
        };

        let pval: f64 = 1.0 - (idx as f64 / arr.len() as f64);
        (pval, mut_len as i32)
    }

    fn probability_based_on_count(
        mutation: &String,
        ctrl_counts: &HashMap<String, u32>,
        edit_counts: &HashMap<String, u32>,
        total_ctrl: i32,
        total_edit: i32,
    ) -> f64 {
        let a = ctrl_counts.get(mutation).cloned().unwrap_or(0) as f64;
        let c = edit_counts.get(mutation).cloned().unwrap_or(0) as f64;
        let b = total_ctrl as f64 - a;
        let d: f64 = total_edit as f64 - c;

        if a / (total_ctrl as f64) >= c / (total_edit as f64) {
            return 1.0;
        }

        let row1 = a + b;
        let row2: f64 = c + d;
        let col1 = a + c;
        let col2 = b + d;
        let N = row1 + row2;
        let exp: [f64; 4] = [
            row1 * col1 / N,
            row1 * col2 / N,
            row2 * col1 / N,
            row2 * col2 / N,
        ];

        if exp.iter().any(|&e| e <= 5.0) {
            return fishers_exact(&[a as u32, b as u32, c as u32, d as u32])
                .unwrap()
                .two_tail_pvalue;
        } else {
            let table = vec![vec![a, b], vec![c, d]];
            return chi2_contingency(&table);
        }
    }

    let mut significant_keys: HashMap<String, (f64, u32)> = HashMap::new();
    let total = treated1.len() as i32;
    let interval = total / 100;
    let (ctrl_del_sroted, ctrl_ins_sorted) = lengths_from_control(&control1);
    let mut records: Vec<(&String, f64, i32, i32, f64, f64, f64, bool)> = Vec::new(); //mutations, p_len, cnt_edit, cnt_con, p_cnt, p_double, bh, bh_flase, sig
    let mut len_df: i32 = 0;

    for (idx, (mut_key, cnt_edit)) in treated1.iter().enumerate() {
        let (p_len, mut_len) =
            probability_based_on_length(mut_key, &ctrl_ins_sorted, &ctrl_del_sroted);
        let cnt_ctrl = control1.get(mut_key).cloned().unwrap_or(0) as i32;
        let cnt_edit = treated1.get(mut_key).cloned().unwrap_or(0) as i32;
        let p_cnt = probability_based_on_count(
            mut_key,
            &control1,
            &treated1,
            control_used_cnt,
            treated_used_cnt,
        );

        let mut double_min: f64;

        if (cnt_ctrl == 0 && mut_len > length_min) || (wo_control) {
            double_min = -1.0;
        } else {
            double_min = p_cnt.clone();
            len_df += 1;
        }

        records.push((
            mut_key, p_len, cnt_ctrl, cnt_edit, p_cnt, double_min, 0.0, false,
        ));
    }

    records.sort_by(|a, b| a.5.partial_cmp(&b.5).unwrap());

    let mut n = 0;
    for record in &mut records {
        if record.5 == -1.0 {
            continue;
        }
        n += 1;
        record.6 = (n as f64 / len_df as f64) * alpha;
        record.7 = record.5 <= record.6;
    }

    if n != 0 {
        let mut p_cutoff: f64 = 0.0;
        for i in &records {
            if i.5 != -1.0 && i.7 == true {
                if i.5 > p_cutoff {
                    p_cutoff = i.5;
                }
            }
        }
        for i in &records {
            if (i.5 != -1.0 && i.5 <= p_cutoff) || (i.5 == -1.0) {
                significant_keys.insert(i.0.to_string(), (i.5, i.3 as u32));
            }
        }
    } else {
        for i in &records {
            if (i.5 == -1.0) {
                significant_keys.insert(i.0.to_string(), (i.5, i.3 as u32));
            }
        }
    }

    significant_keys
}

fn create_reads_results(
    read_to_mut: HashMap<String, Vec<MutationKeys>>,
    significant_keys: Vec<MutationKeys>,
    induced_mutations: Vec<MutationKeys>,
    cv_pos: i32,
    cv_pos_2: Option<i32>,
    window: i32,
) -> Vec<(String, Vec<MutationKeys>, i32)> {
    let mut read_to_mut_final: Vec<(String, Vec<MutationKeys>, i32)> = Vec::new();

    for (read_id, mutations) in read_to_mut {
        let mut filtered_mutations: Vec<MutationKeys> = Vec::new();
        let mut induced_count: i32 = 0;

        for mutation in &mutations {
            if significant_keys.contains(mutation) {
                filtered_mutations.push(mutation.clone());
            }
            if induced_mutations.contains(mutation) {
                induced_count += 1;
            }
        }

        println!("{:?}", filtered_mutations);
        read_to_mut_final.push((read_id, filtered_mutations, induced_count));
    }

    read_to_mut_final
}

fn classify_mutations(
    read_to_mut_final: Vec<(String, Vec<MutationKeys>, i32)>,
    induced_mutations_len: i32,
) -> (HashMap<String, i32>, HashMap<i32, i32>) {
    let mut final_string_result: HashMap<String, i32> = HashMap::new();
    let mut final_cnt: HashMap<i32, i32> = HashMap::new();

    for key in 0..=8 {
        final_cnt.insert(key, 0);
    }

    for (read_id, mutations, induced_count) in read_to_mut_final {
        let mut mutations_string = "".to_string();
        let mut mutation_pat: Vec<i8> = Vec::new(); // 0 Sub, 1 smIns, 2 smDel, 3 LargeIns, 4 LargeDel
        let mutations_cnt = mutations.len() as i32;

        for mutation in mutations {
            if mutation.mutation == 0 {
                let formatted = format!(
                    "{}_{}:Sub_{}_{}>{}",
                    mutation.start_ref,
                    mutation.end_ref,
                    mutation.mutation_len,
                    mutation.reference_seq,
                    mutation.mutated_seq
                );
                mutation_pat.push(0);
                mutations_string.push_str(&formatted);
            } else if mutation.mutation == 1 {
                if mutation.mutation_len > 20 {
                    let formatted = format!(
                        "{}_{}:LargeIns_{}{}",
                        mutation.start_ref,
                        mutation.end_ref,
                        mutation.mutation_len,
                        mutation.mutated_seq
                    );
                    mutation_pat.push(3);
                    mutations_string.push_str(&formatted);
                } else {
                    let formatted = format!(
                        "{}_{}:Ins_{}{}",
                        mutation.start_ref,
                        mutation.end_ref,
                        mutation.mutation_len,
                        mutation.mutated_seq
                    );
                    mutation_pat.push(1);
                    mutations_string.push_str(&formatted);
                }
            } else if mutation.mutation == 2 {
                if mutation.mutation_len > 50 {
                    let formatted = format!(
                        "{}_{}:LargeDel_{}",
                        mutation.start_ref, mutation.end_ref, mutation.mutation_len
                    );
                    mutation_pat.push(4);
                    mutations_string.push_str(&formatted);
                } else {
                    let formatted = format!(
                        "{}_{}:Del_{}",
                        mutation.start_ref, mutation.end_ref, mutation.mutation_len
                    );
                    mutation_pat.push(2);
                    mutations_string.push_str(&formatted);
                }
            }
        }

        if induced_mutations_len != 0 {
            if (induced_count as f64) > (induced_mutations_len as f64 * 0.8) {
                if induced_count == induced_mutations_len {
                    if induced_count == mutations_cnt {
                        mutations_string.push_str(";1");
                    } else {
                        mutations_string.push_str(";2");
                    }
                } else {
                    if induced_count == mutations_cnt {
                        mutations_string.push_str(";3");
                    } else {
                        mutations_string.push_str(";4");
                    }
                }
            } else {
                mutations_string.push_str(";0");
            }
        }

        let mut final_pattern = 0; //1, Sub, 2 smins, 3 smdel, 4 LI, 5 LD, 6 Complex, 7 Inv

        if mutation_pat.is_empty() {
            final_pattern = 0;
        } else {
            if mutation_pat.contains(&4) || mutation_pat.contains(&3) {
                if mutation_pat.contains(&4) && mutation_pat.contains(&3) {
                    final_pattern = 6;
                } else if mutation_pat.contains(&4) {
                    final_pattern = 5;
                } else if mutation_pat.contains(&3) {
                    final_pattern = 4;
                }
            } else if mutation_pat.contains(&1) && mutation_pat.contains(&2) {
                final_pattern = 6;
            } else if mutation_pat.contains(&2) {
                final_pattern = 3;
            } else if mutation_pat.contains(&1) {
                final_pattern = 2;
            } else {
                if mutation_pat.contains(&0) {
                    final_pattern = 1;
                } else {
                    final_pattern = 6;
                }
            }
        }

        if let Some(count) = final_cnt.get_mut(&final_pattern) {
            *count += 1;
        }

        mutations_string = format!("{};{}", final_pattern, mutations_string);

        final_string_result
            .entry(mutations_string)
            .and_modify(|v| *v += 1)
            .or_insert(1);
    }

    (final_string_result, final_cnt)
}

pub fn analysis_function(
    control_align_res: Vec<String>,
    treated_align_res: Vec<String>,
    ref_seq: &str,
    cv_pos: i32,
    cv_pos_2: Option<i32>,
    window: i32,
    whole_window_between_targets: bool,
    induced_align_res: Vec<String>,
    filter1: bool,
    window_filter: bool,
    wo_control: bool,
    alpha: f64,
    length_min: i32,
) -> (
    HashMap<String, i32>,
    HashMap<String, i32>,
    Vec<String>,
    HashMap<String, (f64, u32)>,
    HashMap<String, i32>,
    HashMap<String, Vec<Vec<String>>>,
) {
    //(HashMap<String, i32>, Vec<AnalysisResult>)

    let mut align_summary_cnt: HashMap<String, i32> = HashMap::new();
    let mut analysis_res: HashMap<String, i32> = HashMap::new();

    let (induced_align_cnt, induced_mutations_hashmap, induced_read_to_mut) = quant_unique_indels(
        induced_align_res.clone(),
        ref_seq,
        false,
        cv_pos,
        cv_pos_2,
        window,
        whole_window_between_targets,
    );
    let induced_mutations: Vec<String> = induced_mutations_hashmap.keys().cloned().collect();

    std::mem::drop(induced_align_cnt);
    std::mem::drop(induced_read_to_mut);
    std::mem::drop(induced_mutations_hashmap);
    let (treated_align_cnt, treated_mut, treated_read_to_mut) = quant_unique_indels(
        treated_align_res.clone(),
        ref_seq,
        true,
        cv_pos,
        cv_pos_2,
        window,
        whole_window_between_targets,
    );

    println!("Analyze control alignments");
    let (control_align_cnt, control_mut, control_read_to_mut) = quant_unique_indels(
        control_align_res.clone(),
        ref_seq,
        true,
        cv_pos,
        cv_pos_2,
        window,
        whole_window_between_targets,
    );
    println!("Control reads : {:?}", control_align_cnt);
    std::mem::drop(control_read_to_mut);
    println!("Analyze treated alignments");
    println!("Treated reads : {:?}", treated_align_cnt);

    println!("Pooling... Control");
    let control1: HashMap<String, u32> = pool(control_mut.clone(), 0.00);
    println!("Pooling... Treated");
    let treated1: HashMap<String, u32> = pool(treated_mut.clone(), 0.00);
    println!("Statistically validate the mutation patterns...");
    let mut significant_keys = significant_mutations(
        control1,
        treated1,
        control_align_cnt["Used"],
        treated_align_cnt["Used"],
        0.002,
        wo_control,
        alpha,
        length_min,
    );

    for i in &induced_mutations {
        significant_keys.insert(i.to_string(), (0.0, 0 as u32));
    }

    println!("filtered out unsignificant mutations ...");
    //let read_to_mut_final = create_reads_results(treated_read_to_mut, significant_keys, induced_mutations, cv_pos, cv_pos_2, window);

    //let (final_string_result, final_cnt) = classify_mutations(read_to_mut_final, induced_mutations_len);

    println!("{:?}", treated_align_cnt);

    (
        control_mut,
        treated_mut,
        induced_mutations,
        significant_keys,
        treated_align_cnt,
        treated_read_to_mut,
    )
}
