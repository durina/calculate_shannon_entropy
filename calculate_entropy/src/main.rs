/*
Calculate shannon entropy at each position for a given sequence alignment
    + check if all the alignments are of equal length
    + count the occurrence of each character by position
        + handle characters
        + calculate % of ATCG in a column
            + report entropy only if the column > x% ATGC
            + default 80

Implement clap to parse cli

Libs
    get arguments
    check_format: enforce fasta only format
        len_check: read input file and check if the lengths are identical
    index: build a map of the seek location of the first nucleotide
    count: count the occurence of unique characters across a column
    report: calculate the shannon entropy of each column and export along with column position

Arguments
    get path to alignment file
    flag to include other alphabets
*/

mod utils;
use clap::Parser;
use utils::get_args::Cli;
use utils::check_ali_format::check_fasta;
use utils::calculate_entropy::report_entropy;
use log::{debug, error, info};
use env_logger;
use utils::struct_helper::FileBufferHelper;

fn main() {
    // Path to alignment file 
    // Mode of operation
    // Characters to ignore
    env_logger::init();
    let cli = Cli::parse();
    debug!("Parsing commandline arguments");
    for file in cli.input_alignment {
        debug!("Processing file: {:?}", file);
        if check_fasta(&file) {
            let mut alignment_file = FileBufferHelper::new(&file);
            info!("Alignment complies requirements {:?}", file);
            report_entropy(&mut alignment_file, cli.mode, cli.threshold, cli.nproc, &cli.output_suffix);
        } else {
            error!("Alignment failed");
        }
    }
}
