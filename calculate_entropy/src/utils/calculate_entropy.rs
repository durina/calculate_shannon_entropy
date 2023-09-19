/*
    Calculate the shannon entropy at every position
        Seek to the position
        Read one character to the buffer
        maintain counts of the character at every position
    
    If Mode::Standard
        Calculate the occurence of A, T, G, C, -, and others
    Else if Mode::All
        Calculate the occurence of all characters
*/

use std::{collections::HashMap, path::PathBuf, fs::File, io::{self, BufRead, BufReader, Seek, SeekFrom, BufWriter, Write}, sync::{Mutex, Arc}};
use super::get_args::Mode;
use threadpool::ThreadPool;
use log::{debug, error, warn, info};

const STANDARD_DNA_NOTATIONS: &str = "ATGC";
const ALL_DNA_NOTATIONS: &str = "ATGCUWSMKRYBDHVN";
const DNA_ALIGNMENT_NOTATIONS: &str = "-.";

pub fn report_entropy(file: &PathBuf, mode: Mode, threshold: f64, nproc: usize, suffix: &String) {
    let alignment_file = File::open(&file).unwrap();
    let mut read_alignment_file = io::BufReader::new(&alignment_file);
    debug!("Opened the file buffer");
    let count_vec: Vec<HashMap<char, u32>> = initialise_structs(&file, &mode, &mut read_alignment_file);
    info!("Positions initialised");
    read_alignment_file.seek(SeekFrom::Start(0)).unwrap();
    debug!("Reset file buffer position to start");
    process_genomes(count_vec, threshold, mode, nproc, &mut read_alignment_file, suffix, &file);
}

// Initialise the each position in the aligment

fn initialise_structs(file: &PathBuf, mode: &Mode, buffer: &mut BufReader<&File>) -> Vec<HashMap<char, u32>> {
    debug!("Preparing the columns for analysis");
    let mut line = String::new();
    let mut temp_genome = String::new();
    let mut temp_genome_length = 0usize;
    let mut flag: bool = false;
    // get one sample genome sequence
    info!("Assessing length of alignments in {:?}", file);
    while buffer.read_line(&mut line).unwrap_or(0) > 1 {
        if line.starts_with('\n') {
            flag = true
        } else if line.starts_with('>') && flag {
            info!("While loop terminated");
            line.clear();
            break
        } else if line.starts_with('>') {
            debug!("{}", line.trim());
            flag = true;
        } else {
            temp_genome_length += line.trim().len();
            temp_genome = temp_genome + &line.trim();
            flag = true;
        }
        line.clear()
    }
    if temp_genome.len() != temp_genome_length {
        error!("Check the discrepancy in the length estimated from \
                genome and character count {} {} respectively", temp_genome.len(), temp_genome_length);
    } else {
        info!("Length of every genome in this alignment \
                is {}, and the values agree", temp_genome_length);
    }
    // assign alphabets to be considered
    let alphabets: String = if mode == &Mode::All {
        format!("{}{}" , ALL_DNA_NOTATIONS, DNA_ALIGNMENT_NOTATIONS)
    } else {
        format!("{}{}" , STANDARD_DNA_NOTATIONS, DNA_ALIGNMENT_NOTATIONS)
    };
    // initialise the structs
    info!("Initialising hasmap with {} each with 0", alphabets);
    let hashmap_init_vals = alphabets.chars().map(|n| (n, 0u32));
    let count_vec = (0..temp_genome.len())
                                            .map(|_: usize| HashMap::from_iter(hashmap_init_vals.clone()))
                                            .collect::<Vec<HashMap<char, u32>>>();
    count_vec
}

// Tabulate the frequency of each notation at the given position

fn process_genomes(count_vec: Vec<HashMap<char, u32>>, threshold: f64, mode: Mode, nproc: usize, buffer: &mut BufReader<&File>, suffix: &String, file: &PathBuf) {
    let mut line = String::new();
    let mut temp_genome = String::new();
    let mut header = String::new();
    let arc_count_vec: Arc<Mutex<Vec<HashMap<char, u32>>>> = Arc::new(Mutex::new(count_vec));
    let pool = ThreadPool::new(nproc);
    let mut genome_count: u32 = 0u32;
    // analyse all genomes
    while buffer.read_line(&mut line).unwrap_or(0) > 1 {
        if line.starts_with('>') {
            debug!("Analysing {}", line);
            if header.is_empty() && temp_genome.is_empty() {
                header =  line.trim().to_string();
            } else if !temp_genome.is_empty() {
                debug!("Processing {}", header);
                let arc_clone = Arc::clone(&arc_count_vec);
                let temp_string = temp_genome.drain(..).collect();
                let temp_header = header.drain(..).collect();
                pool.execute(move || {
                    analyse_genomes(temp_string, arc_clone, temp_header);
                });
                header =  line.trim().to_string();
            }
            genome_count += 1;
        } else if line.starts_with('\n') {
            debug!("Empty line");
            line.clear();
            continue;
        } else {
            temp_genome = temp_genome + &line.trim()
        }
        line.clear();
    }
    pool.join();
    info!("Threadpool jobs complete");
    let final_vec = Arc::try_unwrap(arc_count_vec).unwrap().into_inner().unwrap();
    info!("Arc and mutex unwrapped successfully");
    finalise_counts(final_vec, genome_count, threshold, mode, suffix, file)
}

// Handle multiple thread requests, call update_counts

fn analyse_genomes(genome: String, arc_clone: Arc<Mutex<Vec<HashMap<char, u32>>>>, header: String) {
    if let Ok(mut count_vec) = arc_clone.lock() {
        genome.chars().enumerate().for_each( |(idx, x)| {
                update_counts(&mut count_vec[idx], &x, &idx, &header)
            }
        );
    }
}

// update counts at given location
fn update_counts(column: &mut HashMap<char, u32>, letter: &char, position: &usize, header: &String) {
    if let Some(val) = column.get_mut(&letter.to_ascii_uppercase()) {
        *val += 1
    } else {
        warn!("Position: {position} in {header} contains non-permissible character: {letter}");
        if let Some(val) = column.get_mut(&'.') {
            *val += 1
        }
    }
}

// arrive at Shannon entropy at each position

fn finalise_counts(map_vec: Vec<HashMap<char, u32>>, genome_count: u32, threshold: f64, mode: Mode, suffix: &String, file: &PathBuf) {
    // calculate the shannon entropy at every position
    // shannon entropy = sum(-p log_2 p)
    let atgc: &str = if mode == Mode::Standard {
        STANDARD_DNA_NOTATIONS
    } else {
        ALL_DNA_NOTATIONS
    };
    info!("Output file: {}{suffix}", file.to_str().unwrap());
    let mut entropy_writer = BufWriter::new(File::create(format!("{}{}", file.to_str().unwrap(), suffix)).unwrap());
    let mut atgc_count_vec: Vec<&u32> = Vec::with_capacity(atgc.len());
    info!("Notations considered to calculate Shannon entropy: {}",  atgc);
    let count_headers = atgc.chars().fold(String::new(), |fin_str, x| fin_str + "\tCount_" + &x.to_string());
    let headers = format!("Position\
                        {count_headers}\t\
                        Genome_count\t\
                        Notation_share\t\
                        Fraction_notations\t\
                        Shannon_entropy\t\
                        Validity");
    writeln!(entropy_writer, "{}", headers).unwrap();
    map_vec.iter().enumerate().for_each( | (idx, map) | {
            atgc_count_vec = atgc.chars().map(|n| map.get(&n).unwrap())
                            .collect::<Vec<&u32>>();
            let atgc_share: f64 = atgc_count_vec.iter().map(|&&x| f64::from(x)).sum(); 
            let atgc_fraction = atgc_share/f64::from(genome_count);
            let entropy: f64 = get_entropy(&atgc_count_vec, &genome_count);
            let n_counts = atgc_count_vec.drain(..).map(|x| format!("{x}\t")).collect::<String>();
            let validity = if atgc_fraction >= threshold {
                format!("Valid. Threshold = {}", threshold)
            } else {
                format!("Invalid. Threshold = {}", threshold)
            };
            writeln!(entropy_writer, "{pos}\t{counts}{genome_count}\t{atgc_share}\t{atgc_fraction}\t{shannon}\t{validity}", pos=idx+1, counts=n_counts,shannon=entropy).unwrap();
        }
    );
}

fn get_entropy(notation_count: &Vec<&u32>, genome_count: &u32) -> f64 {
    notation_count.iter().map(|&&n| {
            let p = f64::from(n)/f64::from(*genome_count);
            -1.0*p*p.log(2.0)
        }
    ).sum() // sum (-plogp) where p = N/sum
}