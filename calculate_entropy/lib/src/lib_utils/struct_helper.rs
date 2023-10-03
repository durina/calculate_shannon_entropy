
// struct to handle file buffers


use std::fs::File;
use std::io::{BufReader, Seek, SeekFrom};
use std::path::PathBuf;
use log::debug;

pub struct FileBufferHelper<'a> {
    pub path: &'a PathBuf,
    pub buffer_reader: BufReader<File>,
    pub line: String
}

impl<'a> FileBufferHelper<'a> {
    pub fn new(file: &'a PathBuf) -> FileBufferHelper {
        // initialise instant of FileBufferHelper
        let line = String::new();
        debug!("FileHelper created for: {:?}", file);
        let file_open = File::open(file.clone()).unwrap();
        Self {
            path: &file,
            buffer_reader: BufReader::new(file_open),
            line,
        }
    }

    pub fn buffer_reset(&mut self) {
        // reset buffer to position 0
        // file.buffer_reader.seek(SeekFrom::Start(0)).unwrap();
        self.buffer_reader.seek(SeekFrom::Start(0)).expect("Unable to reset buffer");
    }
}