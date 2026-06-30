use wasm_bindgen::prelude::*;
use serde::{Serialize, Deserialize};
use std::collections::HashMap;

mod CRISPRlungo_regular_new;

#[derive(Serialize, Deserialize)]
pub struct AnalysisResult {
    //pub control_mut: HashMap<String, i32>,
    pub treated_read_to_mut: HashMap<String, Vec<Vec<String>>>,
    pub treated_mut: HashMap<String, i32>,
    pub induced_mutations: Vec<String>,
    pub significant_keys: HashMap<String, (f64, u32)>,
    pub treated_align_cnt: HashMap<String, i32>,
}

#[wasm_bindgen]
pub fn analyze(
    control_sam: &str,
    treated_sam: &str,
    induced_mutation_input: &str,
    ref_seq: &str,
    cv_pos: i32,
    cv_pos_2: Option<i32>,
    window: i32,
    whole_window_between_targets: bool,
    filter1: bool,
    window_filter: bool,
    alpha: f64,
    length_min: i32
) -> JsValue {
    console_error_panic_hook::set_once();

    let mut wo_control;
    if (control_sam ==  "") {
        wo_control = true;
    } else {
        wo_control = false;
    }
    let control_sam_lines: Vec<String> = control_sam.lines().map(|s| s.to_string()).collect();
    let treated_sam_lines: Vec<String> = treated_sam.lines().map(|s| s.to_string()).collect();

    let induced_mutation_vec: Vec<String> = induced_mutation_input.lines().map(|s| s.to_string()).collect();

    let ref_seq_upper = ref_seq.to_uppercase();

    let (control_mut, treated_mut, induced_mutations, significant_keys, treated_align_cnt, treated_read_to_mut) =
        CRISPRlungo_regular_new::analysis_function(
            control_sam_lines,
            treated_sam_lines,
            &ref_seq_upper,
            cv_pos,
            cv_pos_2,
            window,
            whole_window_between_targets,
            induced_mutation_vec,
            filter1,
            window_filter,
            wo_control,
            alpha,
            length_min
        );

    let result: AnalysisResult = AnalysisResult {
        treated_read_to_mut,
        treated_mut,
        significant_keys,
        induced_mutations,
        treated_align_cnt,
    };

    // json_compatible(): serialize Rust maps as JS Objects (not Map),
    // matching the previous JsValue::from_serde output so the JS consumer is unchanged.
    result.serialize(&serde_wasm_bindgen::Serializer::json_compatible()).unwrap()
}
