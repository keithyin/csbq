use std::{collections::HashMap, fs, io::{BufWriter, Write}, thread};

use clap::{self, Parser};

use csbq::{cali_worker_for_hsc, LocusInfo, Model};
use gskits::{
    fastx_reader::{fasta_reader::FastaFileReader, fastq_reader::FastqReader, read_fastx},
    utils::{fastx_file_fmt, FastxFile},
};
use mm2::{build_aligner, targets_to_targetsidx};

#[derive(Debug, Parser, Clone)]
#[command(version, about, long_about=None)]
pub struct Cli {
    #[arg(long = "threads")]
    pub threads: Option<usize>,

    #[arg(long = "mode", help = "smc/hsc")]
    pub mode: String,

    /// query file path
    #[arg(short = 'q')]
    pub query_file: String,

    /// target file path
    #[arg(short = 't')]
    pub target_file: String,

    /// model path
    #[arg(long = "model")]
    pub model_path: String,

    /// output filepath, default None, the output will be ${target}.cali.fq or  ${target}.cali.bam  
    #[arg(short = 'o')]
    pub output_filepath: Option<String>,
}

impl Cli {
    pub fn get_output_filepath(&self) -> String {
        if let Some(oup_) = &self.output_filepath {
            oup_.to_string()
        } else {
            if self.target_file.ends_with(".bam") {
                format!("{}.cali.bam", self.target_file.rsplit_once(".").unwrap().0)
            } else {
                format!("{}.cali.fq", self.target_file.rsplit_once(".").unwrap().0)
            }
        }
    }

}

fn main() {
    let cli = Cli::parse();

    match cli.mode.as_str() {
        "hsc" => {
            hsc_csbq(&cli.query_file, &cli.target_file, &cli.model_path, cli.threads, &cli.get_output_filepath());
        },
        "smc" => {
            
        }
        m => panic!("invalid mode. expected smc/hsc, but got {}", m),
    }

    println!("hello world");
}

pub fn hsc_csbq(query_bam_file: &str, target_file: &str, model_path: &str, threads: Option<usize>, output_filepath: &str) {
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
        &mm2::cli::IndexArgs {
            kmer: Some(15),
            wins: Some(10),
        },
        &mm2::cli::MapArgs {},
        &mm2::cli::AlignArgs::default(),
        &mm2::cli::OupArgs::default(),
        &targets,
    );

    let query_files = vec![query_bam_file.to_string()];

    let mut all_contig_locus_info = idx2target_seq
        .iter()
        .map(|(tid, seq)| {
            (
                *tid,
                seq.as_bytes()
                    .iter()
                    .copied()
                    .enumerate()
                    .map(|(pos, base)| LocusInfo::new(pos, base))
                    .collect::<Vec<_>>(),
            )
        })
        .collect::<HashMap<_, _>>();

    thread::scope(|s| {
        let aligners = &aligners;
        let target2idx = &target_name2idx;
        let query_files = &query_files;
        let idx2target = &idx2target_seq;
        let (qs_sender, qs_recv) = crossbeam::channel::bounded(1000);
        s.spawn(move || {
            mm2::query_seq_sender(query_files, qs_sender);
        });

        let num_threads = threads.unwrap_or(num_cpus::get_physical());
        let (align_res_sender, align_res_recv) = crossbeam::channel::bounded(1000);
        for _ in 0..num_threads {
            let qs_recv_ = qs_recv.clone();
            let align_res_sender_ = align_res_sender.clone();
            s.spawn(move || mm2::align_worker(qs_recv_, align_res_sender_, aligners, target2idx));
        }
        drop(qs_recv);
        drop(align_res_sender);
        let model = Model::new(model_path);
        let calibrated_qual = cali_worker_for_hsc(
            align_res_recv,
            &mut all_contig_locus_info,
            idx2target,
            &model,
        );

        assert!(output_filepath.ends_with("fq"), "{}", output_filepath);
        let file = fs::File::open(output_filepath).unwrap();
        let mut buf_writer = BufWriter::new(file);
        targets.iter().enumerate().for_each(|(tid, target)| {
            let tid = tid as i32;
            let name = &target.name;
            let seq = &target.seq;
            let qual = calibrated_qual.get(&tid).unwrap().iter().map(|q| *q + 33).collect::<Vec<_>>();
            let qual_str = unsafe {
                String::from_utf8_unchecked(qual)
            };
            writeln!(&mut buf_writer, "@{}\n{}\n+{}", name, seq, qual_str).unwrap();
        });

    });
}
