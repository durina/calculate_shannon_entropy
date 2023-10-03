/*
    Calculate the shannon entropy at every position
        Seek to the position
        Read one character to the buffer
        maintain counts of the character at every position
    
    If Mode::Standard
        Calculate the occurrence of A, T, G, C, -, and others
    Else if Mode::All
        Calculate the occurrence of all characters
*/

use std::{collections::HashMap, fs::File};
use std::{io::{BufRead, BufWriter, Write}, sync::{Mutex, Arc}};
use super::get_args::Mode;
use threadpool::ThreadPool;
use log::{debug, error, warn, info};
use check_fasta::lib_utils::struct_helper::FileBufferHelper;
use crate::bin_utils::get_args::Cli;

const STANDARD_DNA_NOTATIONS_UPPER: &str = "ATGC";
const ALL_DNA_NOTATIONS_UPPER: &str = "ATGCUWSMKRYBDHVN";
const ALL_DNA_NOTATIONS_LOWER: &str = "atgcuwsmkrybdhvn";
const DNA_ALIGNMENT_NOTATIONS: &str = "-.";

pub fn report_entropy(file: &mut FileBufferHelper, cli: &Cli) {
    // open validated alignment file from main()
    // initialise the HashMap of DNA notatations for every position of the alignment
    let count_vec: Vec<HashMap<char, f64>> = initialise_structs(file, &cli.mode);
    info!("Positions initialised");
    // reset file buffer to zero
    file.buffer_reset();
    debug!("Reset file buffer position to start");
    // count the occurrence of respective notations at every position
    // count_vec: Position wise count of DNA notations from the alignment
    // threshold: Fraction of positions needed to be filled across a
    //              position in the alignment for the position to be considered
    // mode: Consider only Standard, or All accepted IUPAC notations
    // nproc: # number of processors to be involved
    // read_alignment_file: file buffer mapped to the alignment file
    // suffix: suffix to be added while saving the final output file
    // file: location of the alignment file
    process_genomes(count_vec, file, cli);
}

// Initialise the each position in the alignment
fn initialise_structs(file: &mut FileBufferHelper, mode: &Mode) -> Vec<HashMap<char, f64>> {
    // reading into the first genome
    // estimate length of the alignment
    // initialise Vec of HashMaps
    debug!("Preparing the columns for analysis");
    let mut temp_genome = String::new();
    let mut temp_genome_length = 0usize;
    // get one sample genome sequence
    info!("Assessing length of alignments in {:?}", file.path);
    // read the sequence from file buffer
    // Assumptions:
    // The alignment starts with the header of a particular sequence denoted by ">"
    // The sequence follows the header on the following lines
    // The next sequence is either interrupted by a newline or the next header
    file.buffer_reset();
    while file.buffer_reader.read_line(&mut file.line).unwrap_or(0) >= 1 {
        if &file.line[..1] == ">" && temp_genome_length == 0 {
            continue
        } else if ALL_DNA_NOTATIONS_UPPER.contains(&file.line[..1]) ||
            DNA_ALIGNMENT_NOTATIONS.contains(&file.line[..1]) ||
            ALL_DNA_NOTATIONS_LOWER.contains(&file.line[..1]) {
            temp_genome_length += file.line.trim().len();
            temp_genome += file.line.trim();
        } else if &file.line[..1] == ">" && temp_genome_length != 0 {
            info!("While loop breaking");
            break
        }
        file.line.clear()
    }
    // crosscheck estimate of length of genome
    if temp_genome.len() != temp_genome_length {
        error!("Check the discrepancy in the length estimated from \
                genome and character count {} {} respectively", temp_genome.len(),
                                                                temp_genome_length);
    } else {
        info!("Length of every genome in this alignment \
                is {}, and the values agree", temp_genome_length);
    }
    // assign alphabets to be considered
    let alphabets = match mode {
        Mode::All => {
            format!("{}{}", ALL_DNA_NOTATIONS_UPPER, DNA_ALIGNMENT_NOTATIONS)
        },
        Mode::Standard => {
            format!("{}{}", STANDARD_DNA_NOTATIONS_UPPER, DNA_ALIGNMENT_NOTATIONS)
        }
    };
    // initialise the structs
    info!("Initialising hashmap with {} each with 0", alphabets);
    let hashmap_init_vals = HashMap::from_iter(alphabets
                                                            .chars().map(|n| (n, 0.0f64)));
    vec![hashmap_init_vals; temp_genome_length]
}

// Tabulate the frequency of each notation at the given position
fn process_genomes(count_vec: Vec<HashMap<char, f64>>, file: &mut FileBufferHelper,
                   cli: &Cli) {
    let mut genome = String::new();
    let mut header = String::new();
    let arc_count_vec: Arc<Mutex<Vec<HashMap<char, f64>>>> = Arc::new(
                                                        Mutex::new(
                                                            count_vec));
    let pool = ThreadPool::new(cli.nproc);
    let mut genome_count: f64 = 0.0f64;
    // analyse all genomes
    while file.buffer_reader.read_line(&mut file.line).unwrap_or(0) >= 1 {
        if file.line.starts_with('>') {
            debug!("Analysing {}", file.line);
            if header.is_empty() && genome.is_empty() {
                header =  file.line.trim().to_string();
            } else if !genome.is_empty() {
                debug!("Processing {}", header);
                let arc_clone = Arc::clone(&arc_count_vec);
                let temp_genome = genome.drain(..).collect();
                let temp_header = header.drain(..).collect();
                pool.execute(move || {
                    analyse_genomes(temp_genome, arc_clone, temp_header);
                });
                header =  file.line.trim().to_string();
            }
            genome_count += 1.0;
        } else if ALL_DNA_NOTATIONS_UPPER.contains(&file.line[..1]) ||
            DNA_ALIGNMENT_NOTATIONS.contains(&file.line[..1]) ||
            ALL_DNA_NOTATIONS_LOWER.contains(&file.line[..1]) {
            genome += file.line.trim()
        }
        file.line.clear();
    }
    genome_count += 1.0;
    let arc_clone = Arc::clone(&arc_count_vec);
    pool.execute(move || {analyse_genomes(genome, arc_clone, header)});
    pool.join();
    info!("Threadpool jobs complete");
    let final_vec = Arc::try_unwrap(arc_count_vec).unwrap()
                                                            .into_inner().unwrap();
    info!("Arc and mutex unwrapped successfully");
    finalise_counts(final_vec, genome_count, cli, file)
}
// Handle multiple thread requests, call update_counts

fn analyse_genomes(genome: String, arc_clone: Arc<Mutex<Vec<HashMap<char, f64>>>>,
                                                                        header: String) {
    if let Ok(mut char_map) = arc_clone.lock() {
        genome.chars().enumerate().for_each( |(idx, x)|
                update_counts(&mut char_map[idx], &x, &idx, &header)
        );
    }
}

// update counts at given location
fn update_counts(column: &mut HashMap<char, f64>, letter: &char, position: &usize,
                                                                        header: &String) {
    if let Some(val) = column.get_mut(&letter.to_ascii_uppercase()) {
        *val += 1.0
    } else {
        warn!("Position: {position} in {header} contains non-permissible character: {letter}");
        if let Some(val) = column.get_mut(&'.') {
            *val += 1.0
        }
    }
}

// arrive at Shannon entropy at each position
fn finalise_counts(map_vec: Vec<HashMap<char, f64>>, genome_count: f64,
                   cli: &Cli, file: &mut FileBufferHelper) {
    // calculate the shannon entropy at every position
    // shannon entropy = sum(-p log_2 p)
    let atgc: &str = if cli.mode == Mode::Standard {
        STANDARD_DNA_NOTATIONS_UPPER
    } else {
        ALL_DNA_NOTATIONS_UPPER
    };
    let genome_count_f64: f64 = f64::from(genome_count);
    let out_file_name: String = format!("{}_{}",file.path.to_str().unwrap(), cli.output_suffix);
    info!("Output file: {}", out_file_name);
    let out_file = File::create(out_file_name).expect("Unable to create file");
    let mut entropy_writer = BufWriter::new(out_file);

    // vector to store the counts of characters considered
    let mut atgc_count_vec: Vec<f64> = Vec::with_capacity(atgc.len());
    info!("Notations considered to calculate Shannon entropy: {}",  atgc);

    // headers of columns that contain the values of fraction of each
    // character present in a given position
    let count_headers = atgc.chars()
                                    .fold(String::new(), |final_str, x|
                                                        final_str + "\tCount_" + &x.to_string());
    // Headers of final output file
    let headers = format!("Position\
                        {count_headers}{delim}\
                        Genome_count{delim}\
                        Notation_share{delim}\
                        Fraction_notations{delim}\
                        Shannon_entropy{delim}\
                        Validity", delim=cli.delimiter);
    writeln!(entropy_writer, "{}", headers).expect("Unable to write to file");

    // for each position in the alignment, calculate the % of considered characters
    // and shannon entropy at the given position

    map_vec.iter().enumerate().for_each( | (idx, char_map) | {
            // idx: position
            // char_map: HashMap of characters considered and their counts
            atgc_count_vec = atgc.chars()
                            .map(|n| *char_map.get(&n).unwrap())
                            .collect();
            let atgc_share: f64 = atgc_count_vec.iter().sum();
            let atgc_fraction: f64 = atgc_share/genome_count_f64;
            let entropy: f64 = get_entropy(&atgc_count_vec, &genome_count_f64);
            let n_counts = atgc_count_vec.drain(..)
                                        .map(|x| format!("{x}\t")).collect::<String>();
            let validity = if atgc_fraction >= cli.threshold {
                format!("Valid. Threshold = {}", cli.threshold)
            } else {
                format!("Invalid. Threshold = {}", cli.threshold)
            };
            writeln!(entropy_writer, "{pos}\t{counts}{genome_count}\t{atgc_share}\t\
                                        {atgc_fraction}\t{shannon}\t{validity}",
                     pos=idx+1, counts=n_counts,shannon=entropy).unwrap();
        }
    );
}

fn get_entropy(notation_count: &Vec<f64>, genome_count: &f64) -> f64 {
    notation_count.iter().map(|&count| {
            let p = count/genome_count;
            -1.0*p*p.log(2.0)
        }
    ).sum() // sum (-plogp) where p = N/sum
}