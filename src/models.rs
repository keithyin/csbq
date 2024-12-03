use std::{collections::HashMap, fs};

use crate::{decode_base_ctx, LocusInfo};

#[derive(Debug, Clone, Copy)]
pub enum CigarOp {
    Eq,
    Diff, // diff and the difference base
    Match,
    Del,
    Ins, // insertion base
}

pub trait TModel: Sync + Send {
    fn get_prob_with_ref_called_base(
        &self,
        ref_ctx_enc: u8,
        called_base_enc: u8,
        cigar_op: CigarOp,
    ) -> f32;

    fn get_prob_with_ref_base(&self, _ref_ctx_enc: u8) -> f32 {
        // let num_choices = 4_u32.pow(self.get_ctx_len() as u32);

        // 1.0 / num_choices as f32
        0.25
    }

    fn init_locus_info(&self, ref_seq: &[u8]) -> Vec<LocusInfo> {
        let ctx_len = self.get_ctx_len();

        return match ctx_len {
            1 => ref_seq
                .iter()
                .enumerate()
                .into_iter()
                .map(|(pos, &base)| LocusInfo::new(pos, &vec![base], base))
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
                    LocusInfo::new(pos, &vec![base, next_base], base)
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

                    LocusInfo::new(pos, &vec![pre_base, base, next_base], base)
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
        let join_prob = |base_ctx_enc, locus_info: &LocusInfo| -> f32 {
            // match
            let mut ln_prob = locus_info
                .match_cnt
                .iter()
                .enumerate()
                .map(|(called_base_enc, &cnt)| {
                    // print!("cnt: {} , ", cnt);
                    cnt as f32
                        * self
                            .get_prob_with_ref_called_base(
                                base_ctx_enc,
                                called_base_enc as u8,
                                CigarOp::Match,
                            )
                            .ln()
                })
                .sum::<f32>();

            // ins
            ln_prob += locus_info
                .ins_cnt
                .iter()
                .enumerate()
                .map(|(called_base_enc, &cnt)| {
                    // print!("cnt: {} , ", cnt);

                    cnt as f32
                        * self
                            .get_prob_with_ref_called_base(
                                base_ctx_enc,
                                called_base_enc as u8,
                                CigarOp::Ins,
                            )
                            .ln()
                })
                .sum::<f32>();

            // del
            // print!("cnt: {} , ", locus_info.del_cnt);

            ln_prob += 2.0 // 2.0 for more del penalty
                * locus_info.del_cnt as f32
                * self
                    .get_prob_with_ref_called_base(base_ctx_enc, locus_info.base_enc, CigarOp::Del)
                    .ln();

            ln_prob += self.get_prob_with_ref_base(locus.base_ctx_enc).ln();
            ln_prob.exp()
        };

        let numerator = join_prob(locus.base_ctx_enc, &locus);
        // println!("numerator:{}", numerator);

        let denominator = self
            .get_tot_choices(locus.base_ctx_enc)
            .into_iter()
            .map(|ref_ctx_enc| join_prob(ref_ctx_enc, &locus))
            .sum::<f32>();

        // println!("denominator:{}", denominator);

        numerator / denominator
    }

    fn extract_ref_base_enc_from_base_ctx_enc(&self, base_ctx_enc: u8) -> u8 {
        match self.get_ctx_len() {
            1 => base_ctx_enc,
            2 | 3 => (base_ctx_enc >> 2) & 0b11,
            a => panic!("not a valid ctx len. {}", a),
        }
    }

    fn get_tot_choices(&self, base_ctx_enc: u8) -> Vec<u8>;
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
/// 0->=0.1 means AAA -> deletion 0.3
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
                if line.starts_with("##") || !line.contains("=") {
                    None
                } else {
                    let (left, right) = line.trim().split_once("=").unwrap();
                    Some((left.to_string(), right.parse::<f32>().unwrap()))
                }
            })
            .filter(|v| v.is_some())
            .map(|v| v.unwrap())
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
    fn get_prob_with_ref_called_base(
        &self,
        ref_ctx_enc: u8,
        called_base_enc: u8,
        cigar_op: CigarOp,
    ) -> f32 {
        let key = match cigar_op {
            CigarOp::Match | CigarOp::Eq | CigarOp::Diff => {
                format!("{}->{}", ref_ctx_enc, called_base_enc)
            }
            CigarOp::Ins => format!("{}->+{}", ref_ctx_enc, called_base_enc),
            CigarOp::Del => format!("{}->", ref_ctx_enc),
        };

        let res = *self
            .refcalled2prob
            .get(&key)
            .expect(&format!("key:{} not found", key));
        // println!("{}={}", key, res);
        res
    }

    fn get_ctx_len(&self) -> u8 {
        3
    }

    fn get_tot_choices(&self, base_ctx_enc: u8) -> Vec<u8> {
        (0..4)
            .into_iter()
            .map(|base| base_ctx_enc & 0b11110011 | (base << 2))
            .collect()
    }
}

#[cfg(test)]
mod test {
    use crate::{decode_base_ctx, LocusInfo};

    use super::{TModel, TriGramModel};

    #[test]
    fn test_init_locus_info() {}

    #[test]
    fn test_posterior() {
        let model = TriGramModel::new("model.param");
        let mut locus_info = LocusInfo::new(0, b"AAA", 'A' as u8);
        locus_info.add_match('A' as u8);
        locus_info.add_del();
        locus_info.add_del();
        println!("{}", model.locus_posterior_prob(&locus_info));
    }

    #[test]
    fn test_get_other_choices() {
        let model = TriGramModel::new("model.param");
        let mut encs = model.get_tot_choices(0);
        encs.sort();
        assert!(encs.len() == 4);

        assert_eq!(
            decode_base_ctx(encs[0], model.get_ctx_len() as usize),
            "AAA"
        );

        assert_eq!(
            decode_base_ctx(encs[1], model.get_ctx_len() as usize),
            "ACA"
        );
        assert_eq!(
            decode_base_ctx(encs[2], model.get_ctx_len() as usize),
            "AGA"
        );
        assert_eq!(
            decode_base_ctx(encs[3], model.get_ctx_len() as usize),
            "ATA"
        );
    }
}
