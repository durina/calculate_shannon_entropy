use clap::{Parser, ValueEnum};
use std::ops::RangeInclusive;
use std::path::PathBuf;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    /// Path to Alignment file stored in fasta format
    #[arg(short='i', long="infile", required = true, action=clap::ArgAction::Append)]
    pub input_alignment: Vec<PathBuf>, 
    /// Keep tab of 'All' allowed DNA notations or only allow the 'Standard' ATGC. Recommended: 'Standard'
    #[arg(short='m', long="mode", value_enum, required = true)]
    pub mode: Mode,
    /// Set minimum percentage of 'Standard' ATGC notations to constitute the column. Default: 0.8.
    #[arg(short='t', long="threshold", value_parser=validate_percent, default_value_t=0.8)]
    pub threshold: f64,
    /// Suffix to be appended to the filename when storing the file. Default: "_output.csv"
    #[arg(short='s', long="output-suffix", default_value_t=String::from("_output"))]
    pub output_suffix: String,
    /// Specify delimiter to separate position and entropy. Defalt: ","
    #[arg(short='d', long="delimiter", default_value_t=',')]
    pub delimiter: char,
    /// Specify delimiter to separate position and entropy. Defalt: ","
    #[arg(short='n', long="threads", default_value_t=16)]
    pub nproc: usize
}

const PERCENTAGE: RangeInclusive<f64> = 0f64..=1f64;

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum Mode {
    Standard,
    All
}

fn validate_percent(input_str: &str) -> Result<f64, String> {
    let percent: f64 = input_str
        .parse()
        .unwrap_or(0.8f64);
    if PERCENTAGE.contains(&percent) {
        Ok(percent)
    } else {
        Err(
            format!("Threshold not in the range {} - {}",PERCENTAGE.start(), PERCENTAGE.end())
        )
    }
}