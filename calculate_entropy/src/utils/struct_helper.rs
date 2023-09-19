
// struct to handle file buffers


use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

pub struct FileBufferHelper<'a>{
    pub(crate) path: &'a PathBuf,
    pub buffer_reader: BufReader<&'a File>,
    pub line: String
}

impl FileBufferHelper {
    pub fn new(file: &PathBuf) -> FileBufferHelper {
        let mut line = String::new();
        let file_open = File::open(&file).unwrap();
        let mut read_file = BufReader::new(&file_open);
        Self {
            path: &file,
            buffer_reader: read_file,
            line,
        }
    }
}