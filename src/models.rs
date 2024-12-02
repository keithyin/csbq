use std::{collections::HashMap, fs};

use crate::{LocusInfo, PlpState};

pub trait TModel {
    fn get_prob_with_ref_called_base(&self, ref_ctx_enc: u8, op: PlpState) -> f32;

    fn get_prob_with_ref_base(&self, _ref_ctx_enc: u8) -> f32 {
        let num_choices = 2_u32.pow(self.get_ctx_len() as u32);

        1.0 / num_choices as f32
    }

    fn init_locus_info(&self, ref_seq: &[u8]) -> Vec<LocusInfo> {
        let ctx_len = self.get_ctx_len();

        return match ctx_len {
            1 => ref_seq
                .iter()
                .enumerate()
                .into_iter()
                .map(|(pos, &base)| LocusInfo::new(pos, &vec![base]))
                .collect(),

            2 => ref_seq
                .iter()
                .enumerate()
                .into_iter()
                .map(|(pos, &base)| {
                    let next_base = if (pos + 1) >= ref_seq.len() {
                        'A' as u8
                    } else {
                        unsafe { *ref_seq.get_unchecked(pos + 1) }
                    };
                    LocusInfo::new(pos, &vec![base, next_base])
                })
                .collect(),

            3 => ref_seq
                .iter()
                .enumerate()
                .into_iter()
                .map(|(pos, &base)| {
                    let pre_base = if pos == 0 {
                        'A' as u8
                    } else {
                        unsafe { *ref_seq.get_unchecked(pos - 1) }
                    };
                    let next_base = if (pos + 1) >= ref_seq.len() {
                        'A' as u8
                    } else {
                        unsafe { *ref_seq.get_unchecked(pos + 1) }
                    };

                    LocusInfo::new(pos, &vec![pre_base, base, next_base])
                })
                .collect(),

            a => panic!("not a valid ctx len. {}", a),
        };
    }

    // max
    fn get_ctx_len(&self) -> u8;

    ///
    fn shift_position_for_insertion(&self) -> i64 {
        return match self.get_ctx_len() {
            1 => 1,
            2 => 1,
            3 => 1,
            a => panic!("not a valid ctx len. {}", a),
        };
    }

    fn locus_posterior_prob(&self, locus: &LocusInfo) -> f32 {
        if locus.plp_infos.len() == 0 {
            return 1e-10;
        }

        let join_prob = |ref_ctx_enc: u8, plp_states: &Vec<PlpState>| -> f32 {
            let v = plp_states
                .iter()
                .map(|plp_state| {
                    self.get_prob_with_ref_called_base(ref_ctx_enc, *plp_state)
                        .ln()
                })
                .sum::<f32>()
                + self.get_prob_with_ref_base(locus.base_ctx_enc).ln();
            v.exp()
        };

        let numerator = join_prob(locus.base_ctx_enc, &locus.plp_infos);

        let denominator = locus
            .get_other_choices()
            .into_iter()
            .map(|ref_ctx_enc| join_prob(ref_ctx_enc, &locus.plp_infos))
            .sum::<f32>();

        numerator / denominator
    }

}

// pub struct UnigramModel {
//     refcalled2prob: HashMap<String, f32>,
// }

// impl UnigramModel {
//     pub fn new(fname: &str) -> Self {
//         let refcalled2prob = fs::read_to_string(fname)
//             .unwrap()
//             .split("\n")
//             .into_iter()
//             .map(|line| {
//                 let (left, right) = line.trim().split_once(" ").unwrap();
//                 (left.to_string(), right.parse::<f32>().unwrap())
//             })
//             .collect::<HashMap<_, _>>();

//         Self { refcalled2prob }
//     }

//     pub fn get_prob(&self, key: &str) -> f32 {
//         *self.refcalled2prob.get(key).unwrap()
//     }
// }

// impl TModel for UnigramModel {
//     fn get_prob_with_ref_called_base(&self, ref_ctx: &Vec<u8>, called_base: u8) -> f32 {
//         let key = unsafe { String::from_utf8_unchecked(vec![ref_ctx[0], called_base]) };
//         self.get_prob(&key)
//     }

//     fn get_prob_with_ref_base(&self, ref_ctx_enc: u8) -> f32 {
//         let key = unsafe { String::from_utf8_unchecked(vec![ref_ctx[0]]) };
//         self.get_prob(&key)
//     }

//     fn get_ctx_len(&self) -> u8 {
//         1
//     }

//     fn shift_position_for_insertion(&self) -> usize {
//         1
//     }
// }

// /// AA 0.9  real A -> called A prob
// /// AC 0.02
// /// AG 0.02
// /// AT 0.02
// /// A- 0.02 real A -> deletion  prob
// /// A+ 0.02 real A -> insertion prob
// /// CA  ..
// /// ..
// /// A  0.25 prior
// /// C  0.25 ..
// /// G  0.25 ..
// /// T  0.25 ..
// pub struct DefaultUnigramModel {}

// impl TModel for DefaultUnigramModel {
//     fn get_prob_with_ref_called_base(&self, ref_ctx: &Vec<u8>, called_base: u8) -> f32 {
//         return if ref_ctx[0] == called_base { 0.9 } else { 0.02 };
//     }

//     fn get_ctx_len(&self) -> u8 {
//         1
//     }

//     fn shift_position_for_insertion(&self) -> usize {
//         1
//     }
// }

/// 0->0=0.5 means AAA -> A 0.5
/// 0->1=0.3 means AAA -> C 0.3
/// 0=0.1 means AAA -> deletion 0.3
/// 0->+0=0.4 means AAA -> ins A 0.4
/// 1->0=0.4 means AAC-> C 0.4
pub struct TriGramModel {
    refcalled2prob: HashMap<String, f32>,
}

impl TriGramModel {
    pub fn new(fname: &str) -> Self {
        let refcalled2prob = fs::read_to_string(fname)
            .unwrap()
            .split("\n")
            .into_iter()
            .map(|line| {
                let line = line.replace(" ", "");
                let (left, right) = line.trim().split_once("=").unwrap();
                (left.to_string(), right.parse::<f32>().unwrap())
            })
            .collect::<HashMap<_, _>>();

        Self {
            refcalled2prob: refcalled2prob,
        }
    }

    pub fn new_empty() -> Self {
        Self {
            refcalled2prob: HashMap::new(),
        }
    }
}

impl TModel for TriGramModel {
    fn get_prob_with_ref_called_base(&self, ref_ctx_enc: u8, op: PlpState) -> f32 {
        let key = match op {
            PlpState::Eq(base_enc) | PlpState::Diff(base_enc) => {
                format!("{}->{}", ref_ctx_enc, base_enc)
            }
            PlpState::Ins(base_enc) => format!("{}->+{}", ref_ctx_enc, base_enc),
            PlpState::Del => format!("{}", ref_ctx_enc),
        };

        *self.refcalled2prob.get(&key).unwrap()
    }

    fn get_ctx_len(&self) -> u8 {
        3
    }
}
