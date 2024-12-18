use std::{
    collections::HashMap,
    env, fs,
    io::{BufWriter, Write},
    thread,
};

use clap::{self, Parser};
use csbq::{
    cali_worker_for_hsc, calibrate_single_contig_use_bayes, collect_plp_info_from_record,
    collect_plp_worker, dump_locus_infos,
    models::{TModel, TriGramModel},
};
use gskits::{
    ds::ReadInfo,
    fastx_reader::{fasta_reader::FastaFileReader, fastq_reader::FastqReader, read_fastx},
    gsbam::bam_record_ext::{BamReader, BamRecordExt, BamWriter},
    pbar,
    samtools::sort_by_tag,
    utils::{command_line_str, fastx_file_fmt, FastxFile},
};
use mm2::{build_aligner, targets_to_targetsidx};
use rust_htslib::bam::{header::HeaderRecord, Header, HeaderView, Read};
use time;

#[derive(Debug, Parser, Clone)]
#[command(version, about, long_about=None)]
pub struct Cli {
    #[arg(long = "threads")]
    pub threads: Option<usize>,

    #[arg(long = "mode", help = "smc/hsc/train")]
    pub mode: String,

    /// query file path
    #[arg(short = 'q')]
    pub query_file: String,

    /// target file path
    #[arg(short = 't')]
    pub target_file: String,

    /// model path
    #[arg(long = "model", help = "if not given, default model will be used")]
    pub model_path: Option<String>,

    /// output filepath, default None, the output will be ${target}.cali.fq or  ${target}.cali.bam  
    #[arg(short = 'o')]
    pub output_filepath: Option<String>,
}

impl Cli {
    pub fn get_output_filepath(&self) -> String {
        if let Some(oup_) = &self.output_filepath {
            oup_.to_string()
        } else {
            if self.mode.eq("train") {
                format!(
                    "{}.csbq.train.csv",
                    self.target_file.rsplit_once(".").unwrap().0
                )
            } else if self.target_file.ends_with(".bam") {
                format!("{}.cali.bam", self.target_file.rsplit_once(".").unwrap().0)
            } else {
                format!("{}.cali.fq", self.target_file.rsplit_once(".").unwrap().0)
            }
        }
    }
}

fn main() {
    let time_fmt = time::format_description::parse(
        "[year]-[month padding:zero]-[day padding:zero] [hour]:[minute]:[second]",
    )
    .unwrap();

    // let timer = tracing_subscriber::fmt::time::OffsetTime::new(time_offset, time_fmt);
    // let timer = tracing_subscriber::fmt::time::LocalTime::new(time_fmt);
    let time_offset =
        time::UtcOffset::current_local_offset().unwrap_or_else(|_| time::UtcOffset::UTC);
    let timer = tracing_subscriber::fmt::time::OffsetTime::new(time_offset, time_fmt);

    tracing_subscriber::fmt::fmt().with_timer(timer).init();

    let cli = Cli::parse();

    match cli.mode.as_str() {
        "hsc" | "train" => {
            hsc_csbq(
                &cli.query_file,
                &cli.target_file,
                cli.model_path.as_ref().map(|v| v.as_str()),
                cli.threads,
                &cli.get_output_filepath(),
                cli.mode.eq("train"),
            );
        }
        "smc" => {
            smc_csbq(
                &cli.query_file,
                &cli.target_file,
                cli.model_path.as_ref().map(|v| v.as_str()),
                cli.threads,
                &cli.get_output_filepath(),
            );
        }
        m => panic!("invalid mode. expected smc/hsc, but got {}", m),
    }
}

pub fn hsc_csbq(
    query_bam_file: &str,
    target_file: &str,
    model_path: Option<&str>,
    threads: Option<usize>,
    output_filepath: &str,
    dump_train_data: bool,
) {
    let targets = match fastx_file_fmt(target_file).unwrap() {
        FastxFile::Fasta => {
            let fa_iter = FastaFileReader::new(target_file.to_string());
            read_fastx(fa_iter)
        }

        FastxFile::Fastq => {
            let fq_iter = FastqReader::new(target_file.to_string());
            read_fastx(fq_iter)
        }
    };

    let target_name2idx = targets_to_targetsidx(&targets);
    let idx2target_seq = target_name2idx
        .iter()
        .map(|(_, (idx, _))| (*idx as i32, &targets[*idx].seq))
        .collect::<HashMap<_, _>>();

    let aligners = build_aligner(
        "map-ont",
        &mm2::params::IndexParams {
            kmer: None,
            wins: None,
        },
        &mm2::params::MapParams::default(),
        &mm2::params::AlignParams::default(),
        &mm2::params::OupParams::default(),
        &targets,
    );

    let query_files = vec![query_bam_file.to_string()];

    let model = if dump_train_data {
        Box::new(TriGramModel::new_empty()) as Box<dyn TModel>
    } else {
        Box::new(TriGramModel::new(
            model_path,
        )) as Box<dyn TModel>
    };

    let mut all_contig_locus_info = idx2target_seq
        .iter()
        .map(|(tid, seq)| (*tid, model.as_ref().init_locus_info(seq.as_bytes())))
        .collect::<HashMap<_, _>>();

    thread::scope(|s| {
        let aligners = &aligners;
        let target2idx = &target_name2idx;
        let query_files = &query_files;
        let (qs_sender, qs_recv) = crossbeam::channel::bounded(1000);
        s.spawn(move || {
            mm2::query_seq_sender(
                query_files,
                qs_sender,
                &&mm2::params::InputFilterParams::default(),
            );
        });

        let num_threads = threads.unwrap_or(num_cpus::get_physical());
        let (align_res_sender, align_res_recv) = crossbeam::channel::bounded(1000);
        for _ in 0..num_threads {
            let qs_recv_ = qs_recv.clone();
            let align_res_sender_ = align_res_sender.clone();
            s.spawn(move || {
                mm2::align_worker(
                    qs_recv_,
                    align_res_sender_,
                    aligners,
                    target2idx,
                    &mm2::params::OupParams::default(),
                )
            });
        }
        drop(qs_recv);
        drop(align_res_sender);

        // let model = if let Some(model_path_) = model_path {
        //     Box::new(UnigramModel::new(model_path_)) as Box<dyn TModel>
        // } else {
        //     Box::new(DefaultUnigramModel {}) as Box<dyn TModel>
        // };

        assert!(
            output_filepath.ends_with("fq") || output_filepath.ends_with("csv"),
            "{}",
            output_filepath
        );
        let file = fs::File::create(output_filepath).unwrap();
        let mut buf_writer = BufWriter::new(file);

        if dump_train_data {
            writeln!(&mut buf_writer, "ref_base_ctx_enc\tbase_enc\tcnt\top").unwrap();
            collect_plp_worker(
                align_res_recv,
                &mut all_contig_locus_info,
                model.as_ref(),
                true,
            );
            let pb = pbar::get_bar_pb(
                "dump locus info".to_string(),
                pbar::DEFAULT_INTERVAL,
                all_contig_locus_info
                    .values()
                    .into_iter()
                    .map(|infos| infos.len())
                    .sum::<usize>() as u64,
            );
            all_contig_locus_info.iter().for_each(|(_, locus_infos)| {
                pb.inc(1);
                dump_locus_infos(locus_infos, &mut buf_writer);
            });
            pb.finish();
        } else {
            let calibrated_qual = cali_worker_for_hsc(
                align_res_recv,
                &mut all_contig_locus_info,
                model.as_ref(),
                true,
            );

            let pb = pbar::get_spin_pb(
                format!("dump result to {}", output_filepath),
                pbar::DEFAULT_INTERVAL,
            );

            targets.iter().enumerate().for_each(|(tid, target)| {
                pb.inc(1);
                let tid = tid as i32;
                let name = &target.name;
                let seq = &target.seq;
                let qual = calibrated_qual
                    .get(&tid)
                    .unwrap()
                    .iter()
                    .map(|q| *q + 33)
                    .collect::<Vec<_>>();
                let qual_str = unsafe { String::from_utf8_unchecked(qual) };
                writeln!(&mut buf_writer, "@{}\n{}\n+\n{}", name, seq, qual_str).unwrap();
            });
            pb.finish();
        }
    });
}

fn smc_csbq(
    query_bam_file: &str,
    target_file: &str,
    model_path: Option<&str>,
    threads: Option<usize>,
    output_filepath: &str,
) {
    assert!(
        target_file.ends_with(".bam"),
        "only the bam target file supported in smc_csbq"
    );

    assert!(
        !output_filepath.eq(target_file),
        "output_filepath can't be target_file"
    );

    tracing::info!("sorting sbr.bam {}", query_bam_file);
    let sorted_sbr = sort_by_tag(query_bam_file, "ch", None);

    tracing::info!("sorting smc.bam {}", target_file);
    let sorted_smc = sort_by_tag(&target_file, "ch", None);

    let bam_records = read_bam(&sorted_smc);

    let target_name2idx = bam_records
        .iter()
        .enumerate()
        .map(|(tid, read_info)| (read_info.name.clone(), (tid, read_info.seq.len())))
        .collect::<HashMap<String, (usize, usize)>>();

    let idx2target_seq = bam_records
        .iter()
        .enumerate()
        .map(|(idx, read_info)| (idx, &read_info.seq))
        .collect::<HashMap<_, _>>();

    let model = Box::new(TriGramModel::new(
        model_path,
    )) as Box<dyn TModel>;

    let model = model.as_ref();

    thread::scope(|s| {
        let sorted_sbr = &sorted_sbr;
        let sorted_smc = &sorted_smc;

        let target_name2idx = &target_name2idx;
        let idx2target_seq = &idx2target_seq;

        let (sbr_and_smc_sender, sbr_and_smc_recv) = crossbeam::channel::bounded(1000);
        s.spawn(move || {
            asts::subreads_and_smc_generator(
                sorted_sbr,
                sorted_smc,
                &mm2::params::InputFilterParams::default(),
                sbr_and_smc_sender,
            );
        });

        let threads = threads.unwrap_or(num_cpus::get());
        let (align_res_sender, align_res_recv) = crossbeam::channel::bounded(1000);
        for _ in 0..(threads / 2) {
            let sbr_and_smc_recv_ = sbr_and_smc_recv.clone();
            let align_res_sender_ = align_res_sender.clone();
            s.spawn(move || {
                asts::align_worker(sbr_and_smc_recv_, align_res_sender_, target_name2idx);
            });
        }
        drop(sbr_and_smc_recv);
        drop(align_res_sender);

        let (cali_qual_sender, cali_qual_receiver) = crossbeam::channel::bounded(1000);
        for _ in 0..(threads / 2) {
            let align_res_recv_ = align_res_recv.clone();
            let cali_qual_sender_ = cali_qual_sender.clone();

            s.spawn(move || {
                for single_channel_align_res in align_res_recv_ {
                    if single_channel_align_res.records.len() == 0 {
                        continue;
                    }
                    let tid = unsafe { single_channel_align_res.records.get_unchecked(0).tid() };

                    let smc_read = idx2target_seq.get(&(tid as usize)).unwrap().as_bytes();

                    let mut smc_read_locus_info = model.init_locus_info(smc_read);

                    for record in &single_channel_align_res.records {
                        collect_plp_info_from_record(record, &mut smc_read_locus_info, model);
                    }

                    let qual = calibrate_single_contig_use_bayes(&smc_read_locus_info, model);
                    cali_qual_sender_.send((tid, qual)).unwrap();
                }
            });
        }
        drop(align_res_recv);
        drop(cali_qual_sender);

        let pb = pbar::get_spin_pb("do calibration".to_string(), pbar::DEFAULT_INTERVAL);
        let smc_name2cali_qual = cali_qual_receiver
            .into_iter()
            .map(|(tid, qual)| {
                pb.inc(1);
                (
                    unsafe { bam_records.get_unchecked(tid as usize).name.clone() },
                    qual,
                )
            })
            .collect::<HashMap<_, _>>();
        pb.finish();

        let mut target_bam_file = BamReader::from_path(target_file).unwrap();
        target_bam_file.set_threads(10).unwrap();
        let mut o_header = Header::from_template(target_bam_file.header());

        let mut hd = HeaderRecord::new(b"PG");
        hd.push_tag(b"ID", "csbq")
            .push_tag(b"PN", "csbq")
            .push_tag(b"CL", &command_line_str())
            .push_tag(b"VN", &env!("CARGO_PKG_VERSION"));

        if let Some(pp) = get_last_pg_from_bam_header(target_bam_file.header()) {
            hd.push_tag(b"PP", &pp);
        }
        o_header.push_record(&hd);

        let mut o_bam_file =
            BamWriter::from_path(output_filepath, &o_header, rust_htslib::bam::Format::Bam)
                .unwrap();
        o_bam_file.set_threads(10).unwrap();

        let pb = pbar::get_spin_pb(
            format!("dump result to {}", output_filepath),
            pbar::DEFAULT_INTERVAL,
        );

        for record in target_bam_file.records() {
            pb.inc(1);
            let mut record = record.unwrap();
            let record_ext = BamRecordExt::new(&record);
            let qname = record_ext.get_qname();
            let seq = record_ext.get_seq();
            if let Some(qual) = smc_name2cali_qual.get(&qname) {
                // let mut record_new = BamRecord::new();
                record.set(qname.as_bytes(), None, seq.as_bytes(), qual);
                record.remove_aux(b"rq").unwrap();
                record
                    .push_aux(
                        b"rq",
                        rust_htslib::bam::record::Aux::Float(baseq2channelq(&qual)),
                    )
                    .unwrap();
                // record_new
                //     .push_aux(b"RG", record.aux(b"RG").unwrap())
                //     .unwrap();
                // record_new
                //     .push_aux(b"ch", record.aux(b"ch").unwrap())
                //     .unwrap();
                // record_new
                //     .push_aux(b"np", record.aux(b"np").unwrap())
                //     .unwrap();

                // record_new
                //     .push_aux(
                //         b"rq",
                //         rust_htslib::bam::record::Aux::Float(baseq2channelq(&qual)),
                //     )
                //     .unwrap();

                // record_new
                //     .push_aux(b"op", record.aux(b"op").unwrap())
                //     .unwrap();
                // record_new
                //     .push_aux(b"sc", record.aux(b"sc").unwrap())
                //     .unwrap();
                // record = record_new;
            } else {
                tracing::warn!("no cli result for qname:{}", qname);
            }

            o_bam_file.write(&record).unwrap();
        }
        pb.finish();
    });
}

fn read_bam(bam_file: &str) -> Vec<ReadInfo> {
    let mut reader = BamReader::from_path(bam_file).unwrap();
    reader.set_threads(4).unwrap();
    reader
        .records()
        .into_iter()
        .map(|record| record.unwrap())
        .map(|record| ReadInfo::from_bam_record(&record, None))
        .collect()
}

fn get_last_pg_from_bam_header(header_view: &HeaderView) -> Option<String> {
    let header = Header::from_template(header_view);
    let header = header.to_hashmap();
    if let Some(pg_info) = header.get("PG") {
        if let Some(last) = pg_info.last() {
            return Some(
                last.get("ID")
                    .expect(&format!("No ID in PG header"))
                    .to_string(),
            );
        } else {
            return None;
        }
    } else {
        return None;
    }
}

fn baseq2channelq(baseq: &[u8]) -> f32 {
    let length = baseq.len() as f64;
    let err_rate = baseq
        .iter()
        .map(|v| *v as f64)
        .map(|v| 10.0_f64.powf(v / -10.0_f64))
        .reduce(|acc, v| acc + v)
        .unwrap()
        / length;
    (1.0_f64 - err_rate) as f32
}

#[cfg(test)]
mod test {
    use crate::baseq2channelq;

    #[test]
    fn test_baseq2channelq() {
        assert!((baseq2channelq(&[10, 10, 10]) - 0.9).abs() < 0.01);
    }
}
