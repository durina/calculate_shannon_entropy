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
mod bin_utils;
use clap::Parser;
use bin_utils::get_args::Cli;
use check_fasta::check_fasta;
use bin_utils::calculate_entropy::report_entropy;
use log::{debug, info};
use env_logger;
fn main() {
    // Path to alignment file 
    // Mode of operation
    // Characters to ignore
    env_logger::init();
    let cli = Cli::parse();
    // debug!("Parsing commandline arguments");
    for file in &cli.input_alignment {
        debug!("Processing file: {:?}", file);
        match check_fasta(&file, true) {
            Ok(mut alignment_file) => {
                info!("Alignment complies requirements {:?}", file);
                report_entropy(&mut alignment_file, &cli);
            },
            Err(e) => eprintln!("{}", e)
        }
    }
}