use std::{collections::HashMap, fs};

use crossbeam::channel::Receiver;
use gskits::{gsbam::bam_record_ext::BamRecord, pbar};
use mm2::AlignResult;
use rust_htslib::bam::ext::BamRecordExtensions;
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
    Eq(u8),
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

    pub fn get_pos(&self) -> usize {
        self.pos
    }
}

pub fn cali_worker_for_hsc(
    recv: Receiver<AlignResult>,
    all_contig_locus_info: &mut HashMap<i32, Vec<LocusInfo>>,
    model: &Model,
    use_pbar: bool,
) -> HashMap<i32, Vec<u8>> {
    let pb = if use_pbar {
        Some(pbar::get_spin_pb(
            "collect plp info".to_string(),
            pbar::DEFAULT_INTERVAL,
        ))
    } else {
        None
    };

    for align_res in recv {
        pb.as_ref().map(|pb_| pb_.inc(1));
        for record in align_res.records {
            let tid = record.tid();
            let single_contig_locus_info = all_contig_locus_info.get_mut(&tid).unwrap();
            collect_plp_info_from_record(&record, single_contig_locus_info);
        }
    }
    pb.as_ref().map(|pb_| pb_.finish());

    let pb = if use_pbar {
        Some(pbar::get_spin_pb(
            " do calibration".to_string(),
            pbar::DEFAULT_INTERVAL,
        ))
    } else {
        None
    };
    let calibrated_qual: HashMap<i32, Vec<u8>> = all_contig_locus_info
        .iter()
        .map(|(tid, single_contig_locus_info)| {
            pb.as_ref().map(|pb_| pb_.inc(1));
            let qual = calibrate_single_contig_use_bayes(single_contig_locus_info, model);
            (*tid, qual)
        })
        .collect::<HashMap<_, _>>();
    pb.as_ref().map(|pb_| pb_.finish());

    calibrated_qual
}

pub fn collect_plp_info_from_record(
    record: &BamRecord,
    single_contig_locus_info: &mut Vec<LocusInfo>,
) {
    if record.is_secondary() || record.is_unmapped() || record.is_supplementary() {
        return;
    }

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
            if locus_info.cur_base == *query_seq.get_unchecked(qpos.unwrap() as usize) {
                // eq
                locus_info.push_plp_state(PlpState::Eq(
                    *query_seq.get_unchecked(qpos.unwrap() as usize),
                ));
            } else {
                // diff
                locus_info.push_plp_state(PlpState::Diff(
                    *query_seq.get_unchecked(qpos.unwrap() as usize),
                ));
            }
        }
        if ref_pos_cur_or_pre == (ref_end as usize - 1) {
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
        .map(|locus_info| single_locus_bayes(locus_info, model))
        .collect::<Vec<_>>();
    qual
}

pub fn join_prob(cur_base: u8, plp_infos: &Vec<PlpState>, model: &Model) -> f32 {
    if plp_infos.len() == 0 {
        return 1e-10;
    }
    let value = plp_infos
        .iter()
        .map(|plp_state| match *plp_state {
            PlpState::Eq(called_base) => model
                .get_prob_with_ref_called_base(cur_base, called_base)
                .ln(),
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

pub fn single_locus_bayes(locus_info: &LocusInfo, model: &Model) -> u8 {
    let numerator = join_prob(locus_info.cur_base, &locus_info.plp_infos, model);
    if numerator < 1e-9 {
        return 0;
    }
    let denominator = ALL_BASES
        .into_iter()
        .map(|cur_base| join_prob(cur_base, &locus_info.plp_infos, model))
        .reduce(|acc, v| {
            // println!("{}", v);
            acc + v
        })
        .unwrap();

    let mut posterior = numerator / denominator;
    // println!("{}/{} = {}", numerator, denominator, posterior);

    posterior = if posterior > (1. - 1e-5) {
        1. - 1e-5
    } else {
        posterior
    };

    let phreq = (-10. * (1. - posterior).log10()) as u8;
    phreq
}

#[cfg(test)]
mod test {
    use crate::{join_prob, single_locus_bayes, LocusInfo, Model, PlpState};

    #[test]
    fn test_join_prob() {
        let cur_base = 'A' as u8;
        let model = Model::new("model/default.txt");
        let plp_infos = vec![PlpState::Eq('A' as u8), PlpState::Eq('A' as u8)];

        let prob = join_prob(cur_base, &plp_infos, &model);
        assert!((prob - 0.2025).abs() < 1e-4);
    }

    #[test]
    fn test_single_locus_bayes() {
        let cur_base = 'A' as u8;
        let model = Model::new("model/default.txt");

        let mut locus_info = LocusInfo::new(0, cur_base);
        locus_info.push_plp_state(PlpState::Eq('A' as u8));

        let q = single_locus_bayes(&locus_info, &model);
        assert_eq!(q, 12);

        locus_info.push_plp_state(PlpState::Diff('C' as u8));
        let q = single_locus_bayes(&locus_info, &model);
        assert_eq!(q, 2);



    }
}
