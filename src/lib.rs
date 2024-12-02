use std::{
    cmp,
    collections::HashMap,
    fs::File,
    io::{BufWriter, Write},
};

use crossbeam::channel::Receiver;
use gskits::{
    gsbam::bam_record_ext::{BamRecord, BamRecordExt},
    pbar,
};
use mm2::AlignResult;
use models::TModel;
use rust_htslib::bam::ext::BamRecordExtensions;
// use num::Float;

pub mod hsc_csbq;
pub mod models;
pub mod smc_csbq;

pub static ALL_BASES: [u8; 4] = ['A' as u8, 'C' as u8, 'G' as u8, 'T' as u8];
pub static INS: u8 = '+' as u8;
pub static DEL: u8 = '-' as u8;

pub const DNA_SEQ_ASCII2IDX: [u8; 128] = [
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
];

pub fn encode_base_ctx(bases: &[u8]) -> u8 {
    assert!(bases.len() <= 4);

    let mut base_ctx_enc = 0;
    bases.iter().for_each(|v| {
        let single = encode_single_base(*v);
        base_ctx_enc <<= 2;
        base_ctx_enc += single;
    });
    base_ctx_enc
}

pub fn encode_single_base(base: u8) -> u8 {
    unsafe { *DNA_SEQ_ASCII2IDX.get_unchecked(base as usize) }
}

pub fn decode_base_ctx(mut enc: u8, mut ctx_len: usize) -> String {
    let mut stack = vec![];
    while ctx_len > 0 {
        let v = enc & 0b11;
        stack.push(ALL_BASES[v as usize]);
        enc >>= 2;
        ctx_len -= 1;
    }
    unsafe { String::from_utf8_unchecked(stack.into_iter().rev().collect()) }
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
    ctx_len: u8,
    base_ctx_enc: u8,
    plp_infos: Vec<PlpState>,
}

impl LocusInfo {
    pub fn new(pos: usize, bases: &[u8]) -> Self {
        Self {
            pos,
            ctx_len: bases.len() as u8,
            base_ctx_enc: encode_base_ctx(bases),
            plp_infos: vec![],
        }
    }

    pub fn push_plp_state(&mut self, plp_state: PlpState) {
        self.plp_infos.push(plp_state);
    }

    pub fn get_other_choices(&self) -> Vec<u8> {
        let tot_choices = 2_u32.pow(self.ctx_len as u32);
        (0..tot_choices as u8)
            .into_iter()
            .filter(|v| *v != self.base_ctx_enc)
            .collect()
    }

    pub fn get_pos(&self) -> usize {
        self.pos
    }
}

pub fn collect_plp_worker(
    recv: Receiver<AlignResult>,
    all_contig_locus_info: &mut HashMap<i32, Vec<LocusInfo>>,
    model: &dyn TModel,
    use_pbar: bool,
) {
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
            collect_plp_info_from_record(&record, single_contig_locus_info, model);
        }
    }
    pb.as_ref().map(|pb_| pb_.finish());

}

pub fn cali_worker_for_hsc(
    recv: Receiver<AlignResult>,
    all_contig_locus_info: &mut HashMap<i32, Vec<LocusInfo>>,
    model: &dyn TModel,
    use_pbar: bool,
) -> HashMap<i32, Vec<u8>> {


    collect_plp_worker(recv, all_contig_locus_info, model, use_pbar);

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
    model: &dyn TModel,
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
                    if locus_info.base_ctx_enc == *query_seq.get_unchecked(qpos_) {
                        // eq
                        locus_info.push_plp_state(PlpState::Eq(encode_single_base(
                            *query_seq.get_unchecked(qpos_),
                        )));
                    } else {
                        // diff
                        locus_info.push_plp_state(PlpState::Diff(encode_single_base(
                            *query_seq.get_unchecked(qpos_),
                        )));
                    }
                }
            } else {
                locus_info.push_plp_state(PlpState::Del);
            }
        } else {
            // for insertion
            let qpos_ = qpos.unwrap() as usize;
            let next_rpos = ref_pos_cur_or_pre as i64 + model.shift_position_for_insertion();
            let next_rpos = cmp::max(0, next_rpos) as usize;
            if next_rpos < ref_end as usize {
                let locus_info = unsafe { single_contig_locus_info.get_unchecked_mut(next_rpos) };

                locus_info.push_plp_state(PlpState::Ins(unsafe {
                    encode_single_base(*query_seq.get_unchecked(qpos_))
                }));
            }
        }
    }
}

pub fn calibrate_single_contig_use_bayes(
    single_contig_locus_info: &Vec<LocusInfo>,
    model: &dyn TModel,
) -> Vec<u8> {
    let qual = single_contig_locus_info
        .iter()
        .map(|locus_info| model.locus_posterior_prob(locus_info))
        .map(|prob| (-10. * (1. - prob).log10()).round() as u8)
        .collect::<Vec<_>>();
    qual
}

/// used for train model
/// csv header is ctx_enc   base_enc    op
pub fn dump_locus_infos(locus_infos: &Vec<LocusInfo>, writer: &mut BufWriter<File>) {
    locus_infos.iter().for_each(|locus_info| {
        locus_info
            .plp_infos
            .iter()
            .for_each(|&plp_info| match plp_info {
                PlpState::Eq(base_enc) | PlpState::Diff(base_enc) => {
                    writeln!(writer, "{}\t{}\tM", locus_info.base_ctx_enc, base_enc).unwrap()
                }
                PlpState::Ins(base_enc) => {
                    writeln!(writer, "{}\t{}\tI", locus_info.base_ctx_enc, base_enc).unwrap()
                }
                PlpState::Del => writeln!(writer, "{}\t\tD", locus_info.base_ctx_enc).unwrap(),
            });
    });
}

#[cfg(test)]
mod test {
    // use crate::{single_locus_bayes, LocusInfo, PlpState};

    // #[test]
    // fn test_join_prob() {
    //     let cur_base = 'A' as u8;
    //     let model = DefaultUnigramModel {};
    //     let plp_infos = vec![PlpState::Eq('A' as u8), PlpState::Eq('A' as u8)];

    //     let prob = join_prob(cur_base, &plp_infos, &model);
    //     assert!((prob - 0.2025).abs() < 1e-4);
    // }

    // #[test]
    // fn test_single_locus_bayes() {
    //     let cur_base = 'A' as u8;
    //     let model = UnigramModel::new("model/default.txt");

    //     let mut locus_info = LocusInfo::new(0, &vec![cur_base]);
    //     locus_info.push_plp_state(PlpState::Eq('A' as u8));

    //     let q = single_locus_bayes(&locus_info, &model);
    //     assert_eq!(q, 12);

    //     locus_info.push_plp_state(PlpState::Diff('C' as u8));
    //     let q = single_locus_bayes(&locus_info, &model);
    //     assert_eq!(q, 2);
    // }

    use crate::{decode_base_ctx, encode_base_ctx};

    #[test]
    fn test_encode_base_ctx() {
        let enc = encode_base_ctx(b"ACG");
        assert_eq!(enc, 6);
        let decoded = decode_base_ctx(enc, 3);
        assert_eq!(&decoded, "ACG");

        let enc = encode_base_ctx(b"TTT");
        assert_eq!(enc, 63);
        let decoded = decode_base_ctx(enc, 3);
        assert_eq!(&decoded, "TTT");
    }
}
