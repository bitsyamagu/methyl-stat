
use std::fs::File;
// use std::io::{self, BufReader};
use std::io::{self, Write};
use std::io::prelude::*;
use std::path::Path;
// use std::path::PathBuf;
use std::io::SeekFrom;

#[derive(Debug)]
pub struct FastaIndex {
   pub indices: Vec<Index>,
   pub fasta_path: String
}

impl FastaIndex {
    fn push(&mut self, idx: Index){
        self.indices.push(idx);
    }
    pub fn get_length(&self, chr: &String) -> usize {
        let mut id = 0;
        for idx in &self.indices {
            if idx.name == *chr {
                break;
            }
            id += 1;
        }
        self.indices[id].length
    }
    pub fn get(&self, chr: &String, start: usize, end: usize) -> Result<String, std::string::FromUtf8Error> {
        let mut id = 0;
        for idx in &self.indices {
            if idx.name == *chr {
                break;
            }
            id += 1;
        }
        io::stderr().write_all(&format!("processing {}:{}-{}\n", chr, start, end).into_bytes());
        self.get_by_id(id, start, end)
    }
    pub fn get_by_id(&self, id: usize, start: usize, end_unchecked: usize) -> Result<String, std::string::FromUtf8Error> {
        let chr_offset = self.indices[id].offset;
        let line60 = self.indices[id].line_length;
        let line61 = self.indices[id].line_bytes;
        let seq_start =  chr_offset + (start/line60)*line61 + (start%line60);
        let mut end = end_unchecked;
        if  end > self.indices[id].length {
            end = self.indices[id].length;
        }
        let seq_end =  chr_offset + (end/line60)*line61 + (end%line60);
        let seq_len = seq_end - seq_start + 1;


        // println!("seq_start = {}", seq_start);
        // println!("seq_len = {}", seq_len);

        let mut fasta_file = File::open(&self.fasta_path).unwrap();
        let mut buffer = vec![0; seq_len];
        match fasta_file.seek(SeekFrom::Start(seq_start as u64)) {
            Ok(_n) => {
                fasta_file.read_exact(&mut buffer);
            },
            Err(_e) => {
                panic!("bad region");
            }
        }
        String::from_utf8(buffer)
    }
}
#[derive(Debug)]
pub struct Index {
    name: String,
    length: usize,
    offset: usize,
    line_length: usize,
    line_bytes: usize
}

pub fn read_index(index_path: &String) -> FastaIndex {
    // fasta
    let ipath = index_path.replace(".fai", "");
    /*
    let mut fasta_path = PathBuf::from(index_path);
    fasta_path.set_extension("fa");
    if !fasta_path.exists() {
        fasta_path.set_extension("fasta");
        if !fasta_path.exists() {
            panic!("no fasta file");
        }
    }
    let mut fai: FastaIndex = FastaIndex{ indices: vec![], fasta_path: String::from(fasta_path.to_str().unwrap())};
    */
    let mut fai: FastaIndex = FastaIndex{ indices: vec![], fasta_path: ipath};
    if let Ok(lines) = read_lines(index_path) {
        for line in lines {
            if let Ok(ip) = line {
                // println!("{}", ip);
                let buf: Vec<&str> = ip.split("\t").collect();
                let idx = Index{ name: String::from(buf[0]), length: parse_usize(buf[1]), offset: parse_usize(buf[2]),
                    line_length: parse_usize(buf[3]), line_bytes:  parse_usize(buf[4])};
                // println!("{:?}", fai);
                fai.push(idx);
            }
        }
    }
    fai
}
// fn parse_u32 (buf: &str) -> u32{
//     buf.parse::<u32>().unwrap()
// }
fn parse_usize (buf: &str) -> usize{
    buf.parse::<usize>().unwrap() 
}
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>> where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}
