/*
    Open the file
    Check if file is in fasta format
        Look for ">"
            Count ">"
            Else print not in fasta format
        Check if number of characters in the lines after ">" and before ">" are similar throughout
            Else print alignment not correct
*/


use std::io::{Seek, BufRead, BufReader};
use std::fs::File;
use std::path::PathBuf;
use log::{debug, error, warn, info};

pub fn check_fasta(infile: &PathBuf) -> bool {
    // check if the first line is ">", except empty space
    // lines after empty lines start with ">"
    let alignment_file =  match File::open(infile) {
        Ok(file) => {  
                            info!("File opened successfully - {:?}", infile);
                            file
                        },
        Err(x) => {
                            error!("File {:?} could not be openned - {}", infile, x);
                            return false;
        }
    };
    let mut read_alignment_file = BufReader::new(&alignment_file);
    let mut prev_alignment_length = 0u64;
    let mut first_seq_line_pos = 0u64;
    let mut store_position = false;
    let mut line = String::new();
    let mut found_header = false;
    while read_alignment_file.read_line(&mut line).unwrap_or(0) > 1 {
        match &line[..1] {
            "\n" => {
                if store_position {
                    error!("Empty sequence enountered at {}", read_alignment_file.stream_position().unwrap());
                    return false;
                };
                if !check_alignment_length(&mut read_alignment_file, &found_header, &store_position, &first_seq_line_pos, &mut prev_alignment_length) {
                    return false
                };
                found_header = false;
            }
            ">" => {
                info!("Processing {}", line.trim());
                if !check_alignment_length(&mut read_alignment_file, &found_header, &store_position, &first_seq_line_pos, &mut prev_alignment_length) {
                    return false;
                };
                found_header = true;
                debug!("Accessing: {}", line.trim());
            }
            &_ if !found_header => {
                error!("Sequence interrupted by newline");
                return false;
            }
            &_ if found_header => {
                if store_position {
                    first_seq_line_pos = read_alignment_file.stream_position().unwrap();
                    store_position = false;
                }
            }
            &_ => warn!("Unplanned format: \n{}", line)
        }
    line.clear();
    }
    true
}

fn check_alignment_length(buffer: &mut BufReader<&File>, header_flag: &bool, pos_store_flag: &bool, seq_start_pos: &u64, prev_align_len: &mut u64) -> bool {
    if *header_flag && !*pos_store_flag && seq_start_pos > &0 {
        let current_alignment_length = buffer.stream_position().unwrap() - seq_start_pos;
        if *prev_align_len > 0u64 && &current_alignment_length != prev_align_len {
            error!("Sequence alignment length does not match previous alignment length!");
            false
        } else {
            *prev_align_len = current_alignment_length;
            debug!("Current_alignment_length == prev_align_len == {}", current_alignment_length);
            true
        }
    } else {
        warn!("Situation unmanaged");
        true
    }
}