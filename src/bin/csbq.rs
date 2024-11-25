
use std::thread;

use clap::{self, Parser};


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

    #[arg(long="mode", help="smc/hsc")]
    pub mode: String,

    /// query file path
    #[arg(short='q')]
    pub query_file: String,

    /// target file path
    #[arg(short='t')]
    pub target_file: String,

    /// model path
    #[arg(long="model")]
    pub model_path: String

}



fn main() {

    let cli = Cli::parse();
    

    println!("hello world");
}


pub fn hsc_csbq(query_bam_file: &str, target_file: &str, threads: Option<usize>) {
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

    let target2idx = targets_to_targetsidx(&targets);

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

    thread::scope(|s| {
        let aligners = &aligners;
        let target2idx = &target2idx;
        let query_files = &query_files;
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

        

    });
    
}
