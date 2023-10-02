/*
    Open the file
    Check if file is in fasta format
        Look for ">"
            Count ">"
            Else print not in fasta format
        Check if number of characters in the lines after ">" and before ">" are similar throughout
            Else print alignment not correct
*/


use std::io::{Seek, BufRead};
use std::fs::File;
use std::path::PathBuf;
use log::{debug, error, info, trace, warn};
pub mod lib_utils;
use lib_utils::struct_helper::FileBufferHelper;

const IUPAC_DNA_ALIGNMENT: &str = "ATGCUWSMKRYBDHVN-.";
const IUPAC_DNA_ALIGNMENT_LOWER: &str = "atgcuwsmkrybdhvn";
pub fn check_fasta(infile: &PathBuf, length_check: bool) -> Result<FileBufferHelper, &'static str> {
    // check if the first line is ">", except empty space
    // lines after empty lines start with ">"
    let mut alignment_file =  match File::open(infile) {
        Ok(_) => {
                            info!("File opened successfully - {:?}", infile);
                            FileBufferHelper::new(infile)
                        },
        Err(x) => {
                            error!("File {:?} could not be opened - {}", infile, x);
                            return Err("Alignment incorrect");
        }
    };
    let mut prev_alignment_length = 0u64;
    let mut first_seq_line_pos = 0u64;
    let mut store_position = false;
    let mut current_seq_line_pos = 0u64;
    let mut found_header = false;
    let mut interruption: bool = false;
    let mut header_info = String::new();

    while alignment_file.buffer_reader.read_line(&mut alignment_file.line).unwrap_or(0) >= 1 {
        match &alignment_file.line[..1] {
            // match header
            ">" => {
                if store_position && !found_header && length_check {
                    // if position of genome start is available
                    // and header flag is off
                    // and length check is on
                    // perform length check
                    info!("Checking length of sequence in alignment.");
                    let length = current_seq_line_pos - first_seq_line_pos;
                    if prev_alignment_length != 0 && prev_alignment_length == length {
                        debug!("{} matches alignment length of {}",
                            header_info, prev_alignment_length);

                    } else if prev_alignment_length == 0 {
                        info!("Alignment length set as {}", length);
                        prev_alignment_length = length;
                    } else if prev_alignment_length != length {
                        error!("{} does not match alignment length.", header_info);
                        return Err("Does not match alignment length., check log")
                    }
                } else if found_header {
                    // if header is encountered right after a header
                    // throw error
                    error!("No sequence encountered in between headers. Empty sequence encountered.\\
                    Remove headers without any sequence and try again.");
                    return Err("No sequence encountered in between headers. Empty sequence encountered.\\
                    Remove headers without any sequence and try again.")
                }
                header_info = alignment_file.line.trim().to_string();
                trace!("Processing {}", header_info);
                found_header = true;
                store_position = false;
                interruption = false;
                first_seq_line_pos = alignment_file.buffer_reader.stream_position()
                        .expect("Unable to retrieve buffer position");
            },
            x if IUPAC_DNA_ALIGNMENT.contains(x) || IUPAC_DNA_ALIGNMENT_LOWER.contains(x) => {
                // if the the line starts with any of the IUPAC DNA characters
                if !found_header && !store_position {
                    // if the sequences are found before the corresponding header
                    // throw error
                    error!("Encountered sequences before header");
                    return  Err("Encountered sequences before header")
                } else if found_header && !store_position {
                    // if the first line of sequence is encountered
                    // the location of the pointer stored when the header was encountered
                    // is valid
                    store_position = true;
                    found_header = false;
                } else if !found_header && store_position && !interruption {
                    // if subsequent lines are encountered without any interruptions
                    // store position of pointer
                    // not interrupted by empty lines
                    if length_check {
                        current_seq_line_pos = alignment_file.buffer_reader.stream_position()
                            .expect("Unable to retrieve position of stream.");
                    }
                } else if !found_header && store_position && interruption {
                    // if the sequences are interrupted by newline or
                    // any other non-IUPAC DNA character
                    // throw error
                    error!("Sequence interrupted by newline or non-IUPAC character.");
                    return Err("Sequence interrupted by newline or non-IUPAC character.")
                } else if found_header && store_position {
                    // found_header and store_position should be mutually exclusive
                    // if both are true, throw error
                    error!("Unexplained situation. Both store_position and header are true.\\
                    Possibly interrupted by newline or unidentified character.");
                    return Err("Unexplained situation. Both store_position and header are true.\\
                    Possibly interrupted by newline or unidentified character.")
                }
            },
            "\n" => {
                // if the headers or sequences are interrupted by newline
                // throw error if newlines are found in between
                //  consecutive headers
                //  consecutive lines of sequences
                if found_header {
                    error!("Header interrupted by newline before encountering sequence");
                    return Err("Header interrupted by newline before encountering sequence")
                } else if !found_header && store_position {
                    // two possibilities
                    //  interruption is after the end of sequence
                    //  interruption is in between the sequences
                    interruption = true;
                }
            },
            &_ => {
                // if the line starts with a character that is not ">",
                // IUPAC_DNA_characters or empty
                // throw error
                interruption = true;
                warn!("Unexpected scenarios encountered {}", alignment_file.line);
                // Err(format!("Unexpected scenarios encountered {}", alignment_file.line))
            }
        }
    alignment_file.line.clear();
    }
    Ok(alignment_file)
}