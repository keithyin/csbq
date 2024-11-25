use std::{collections::HashMap, fs};

use gskits::gsbam::bam_record_ext::BamRecord;
use rust_htslib::bam::{self, ext::BamRecordExtensions, Read};
// use num::Float;

pub mod hsc_csbq;
pub mod smc_csbq;

pub static ALL_BASES: [u8; 4] = ['A' as u8, 'C' as u8, 'G' as u8, 'T' as u8];
pub static INS: u8 = '+' as u8;
pub static DEL: u8 = '-' as u8;

/// AA 0.9  real A -> called A prob
/// AC 0.02
/// AG 0.02
/// AT 0.02
/// A- 0.02 real A -> deletion  prob
/// A+ 0.02 real A -> insertion prob
/// CA  ..
/// ..
/// A  0.25 prior
/// C  0.25 ..
/// G  0.25 ..
/// T  0.25 ..
pub struct Model {
    refcalled2prob: HashMap<String, f32>,
}

impl Model {
    pub fn new(fname: &str) -> Self {
        let refcalled2prob = fs::read_to_string(fname)
            .unwrap()
            .split("\n")
            .into_iter()
            .map(|line| {
                let (left, right) = line.trim().split_once(" ").unwrap();
                (left.to_string(), right.parse::<f32>().unwrap())
            })
            .collect::<HashMap<_, _>>();

        Self { refcalled2prob }
    }

    pub fn get_prob(&self, key: &str) -> f32 {
        *self.refcalled2prob.get(key).unwrap()
    }

    pub fn get_prob_with_ref_called_base(&self, ref_base: u8, called_base: u8) -> f32 {
        let key = unsafe { String::from_utf8_unchecked(vec![ref_base, called_base]) };
        self.get_prob(&key)
    }

    pub fn get_prob_with_ref_base(&self, ref_base: u8) -> f32 {
        let key = unsafe { String::from_utf8_unchecked(vec![ref_base]) };
        self.get_prob(&key)
    }
}

#[derive(Debug, Clone, Copy)]
pub enum PlpState {
    Eq,
    Diff(u8), // diff and the difference base
    Del,
    Ins(u8), // insertion base
}
pub struct LocusInfo {
    pos: usize,
    cur_base: u8,
    plp_infos: Vec<PlpState>,
}

impl LocusInfo {
    pub fn new(pos: usize, base: u8) -> Self {
        Self {
            pos,
            cur_base: base,
            plp_infos: vec![],
        }
    }

    pub fn push_plp_state(&mut self, plp_state: PlpState) {
        self.plp_infos.push(plp_state);
    }

    pub fn get_other_bases(&self) -> Vec<u8> {
        ALL_BASES
            .iter()
            .copied()
            .filter(|base| *base != self.cur_base)
            .collect()
    }
}

pub fn info_collector(
    record: &BamRecord,
    single_contig_locus_info: &mut Vec<LocusInfo>,
    ref_seq: &str,
) {
    if record.is_secondary() || record.is_unmapped() || record.is_supplementary() {
        return;
    }

    let refseq_bytes = ref_seq.as_bytes();

    let ref_start = record.reference_start();
    let ref_end = record.reference_end();

    let mut rpos_cursor = None;
    let mut qpos_cursor = None;

    let query_seq = record.seq().as_bytes();

    for [qpos, rpos] in record.aligned_pairs_full() {
        if qpos.is_some() {
            qpos_cursor = qpos;
        }
        if rpos.is_some() {
            rpos_cursor = rpos;
        }

        if let Some(rpos_cursor_) = rpos_cursor {
            if rpos_cursor_ < ref_start {
                continue;
            }
            if rpos_cursor_ >= ref_end {
                break;
            }
        } else {
            continue;
        }

        if qpos_cursor.is_none() {
            continue;
        }

        let ref_pos_cur_or_pre = rpos_cursor.unwrap() as usize;

        let locus_info = unsafe { single_contig_locus_info.get_unchecked_mut(ref_pos_cur_or_pre) };

        if qpos.is_none() {
            // deletion
            locus_info.push_plp_state(PlpState::Del);
            continue;
        }
        if rpos.is_none() {
            // ins
            locus_info.push_plp_state(PlpState::Ins(unsafe {
                *query_seq.get_unchecked(qpos.unwrap() as usize)
            }));
            continue;
        }

        unsafe {
            if *refseq_bytes.get_unchecked(rpos.unwrap() as usize)
                == *query_seq.get_unchecked(qpos.unwrap() as usize)
            {
                // eq
                locus_info.push_plp_state(PlpState::Eq);
            } else {
                // diff
                locus_info.push_plp_state(PlpState::Diff(
                    *query_seq.get_unchecked(qpos.unwrap() as usize),
                ));
            }
        }
        if rpos_cursor.unwrap() as usize == (ref_end as usize - 1) {
            break;
        }
    }
}

pub fn calibrate_single_contig_use_bayes(
    single_contig_locus_info: &Vec<LocusInfo>,
    model: &Model,
) -> Vec<u8> {
    let qual = single_contig_locus_info
        .iter()
        .map(|locus_info| {
            let numerator = join_prob(locus_info.cur_base, &locus_info.plp_infos, model);
            let denominator = locus_info
                .get_other_bases()
                .into_iter()
                .map(|cur_base| join_prob(cur_base, &locus_info.plp_infos, model))
                .reduce(|acc, v| acc + v)
                .unwrap();
            numerator / denominator
        })
        .map(|baseq| {
            if baseq > (1. - 1e-6) {
                1. - 1e-6
            } else {
                baseq
            }
        })
        .map(|baseq| -10. * (1. - baseq).log10())
        .map(|phreq| phreq as u8)
        .collect::<Vec<_>>();
    qual
}

fn join_prob(cur_base: u8, plp_infos: &Vec<PlpState>, model: &Model) -> f32 {
    let value = plp_infos
        .iter()
        .map(|plp_state| match *plp_state {
            PlpState::Eq => model.get_prob_with_ref_called_base(cur_base, cur_base).ln(),
            PlpState::Diff(called_base) => model
                .get_prob_with_ref_called_base(cur_base, called_base)
                .ln(),
            PlpState::Ins(_) => model.get_prob_with_ref_called_base(cur_base, INS).ln(),
            PlpState::Del => model.get_prob_with_ref_called_base(cur_base, DEL).ln(),
        })
        .reduce(|acc, cur| acc + cur)
        .unwrap()
        + model.get_prob_with_ref_base(cur_base).ln();

    value.exp()
}

#[cfg(test)]
mod test {

    #[test]
    fn say_hello() {
        println!("hello world");
    }
}
