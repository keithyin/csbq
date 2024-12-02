use std::{collections::HashMap, fs};

use crossbeam::channel::Receiver;
use gskits::{
    gsbam::bam_record_ext::{BamRecord, BamRecordExt},
    pbar,
};
use mm2::AlignResult;
use models::TUnigramModel;
use rust_htslib::bam::ext::BamRecordExtensions;
// use num::Float;

pub mod hsc_csbq;
pub mod smc_csbq;
pub mod models;

pub static ALL_BASES: [u8; 4] = ['A' as u8, 'C' as u8, 'G' as u8, 'T' as u8];
pub static INS: u8 = '+' as u8;
pub static DEL: u8 = '-' as u8;


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
    model: &dyn TUnigramModel,
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

    let record_ext = BamRecordExt::new(record);
    let ref_start = record_ext.reference_start() as i64;
    let ref_end = record_ext.reference_end() as i64;

    let query_end = record_ext.query_alignment_end() as i64;

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

        if let Some(qpos_cursor_) = qpos_cursor {
            if qpos_cursor_ >= query_end {
                break;
            }
        }

        let ref_pos_cur_or_pre = rpos_cursor.unwrap() as usize;

        if let Some(rpos_) = rpos {
            let rpos_ = rpos_ as usize;
            let locus_info = unsafe { single_contig_locus_info.get_unchecked_mut(rpos_) };
            if let Some(qpos_) = qpos {
                let qpos_ = qpos_ as usize;
                unsafe {
                    if locus_info.cur_base == *query_seq.get_unchecked(qpos_) {
                        // eq
                        locus_info.push_plp_state(PlpState::Eq(*query_seq.get_unchecked(qpos_)));
                    } else {
                        // diff
                        locus_info.push_plp_state(PlpState::Diff(*query_seq.get_unchecked(qpos_)));
                    }
                }
            } else {
                locus_info.push_plp_state(PlpState::Del);
            }
        } else {
            let qpos_ = qpos.unwrap() as usize;
            let next_rpos = ref_pos_cur_or_pre + 1;
            if next_rpos < ref_end as usize {
                let locus_info = unsafe { single_contig_locus_info.get_unchecked_mut(next_rpos) };

                locus_info
                    .push_plp_state(PlpState::Ins(unsafe { *query_seq.get_unchecked(qpos_) }));
            }
        }
    }
}

pub fn calibrate_single_contig_use_bayes(
    single_contig_locus_info: &Vec<LocusInfo>,
    model: &dyn TUnigramModel,
) -> Vec<u8> {
    let qual = single_contig_locus_info
        .iter()
        .map(|locus_info| single_locus_bayes(locus_info, model))
        .collect::<Vec<_>>();
    qual
}

pub fn join_prob(cur_base: u8, plp_infos: &Vec<PlpState>, model: &dyn TUnigramModel,) -> f32 {
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

pub fn single_locus_bayes(locus_info: &LocusInfo, model: &dyn TUnigramModel,) -> u8 {
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
    use crate::{join_prob, models::{DefaultUnigramModel, UnigramModel}, single_locus_bayes, LocusInfo, PlpState};

    #[test]
    fn test_join_prob() {
        let cur_base = 'A' as u8;
        let model = DefaultUnigramModel{};
        let plp_infos = vec![PlpState::Eq('A' as u8), PlpState::Eq('A' as u8)];

        let prob = join_prob(cur_base, &plp_infos, &model);
        assert!((prob - 0.2025).abs() < 1e-4);
    }

    #[test]
    fn test_single_locus_bayes() {
        let cur_base = 'A' as u8;
        let model = UnigramModel::new("model/default.txt");

        let mut locus_info = LocusInfo::new(0, cur_base);
        locus_info.push_plp_state(PlpState::Eq('A' as u8));

        let q = single_locus_bayes(&locus_info, &model);
        assert_eq!(q, 12);

        locus_info.push_plp_state(PlpState::Diff('C' as u8));
        let q = single_locus_bayes(&locus_info, &model);
        assert_eq!(q, 2);
    }
}
