
// struct to handle file buffers


use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

pub struct FileBufferHelper<'a> {
    pub(crate) path: &'a PathBuf,
    pub buffer_reader: BufReader<File>,
    pub line: String
}

impl<'a> FileBufferHelper<'a> {
    pub fn new(file: &'a PathBuf) -> FileBufferHelper {
        let mut line = String::new();
        let file_open = File::open(file.clone()).unwrap();
        Self {
            path: &file,
            buffer_reader: BufReader::new(file_open),
            line,
        }
    }
}