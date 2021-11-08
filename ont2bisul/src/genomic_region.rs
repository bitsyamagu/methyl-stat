use fasta;

pub struct GenomicRegion {
    pub chr: String,
    pub start: u32,
    pub end: u32,
    pub path: String,
    pub sequence: String,
    pub char_seq: Vec<char>
}

impl GenomicRegion {
    pub fn new(chr: &String, start: u32, end: u32, path: &String) -> Self {
        Self { chr: chr.to_string(), start: start, end: end, path: path.to_string(), sequence: String::from(""), char_seq: Vec::<char>::new()}
    }
    pub fn init(&mut self){
        let fai = fasta::read_index(&self.path);
        let genomic_seq = fai.get(&self.chr, self.start as usize, self.end as usize);
        match genomic_seq {
            Ok(mut seq) => {
                &seq.retain(|c| !c.is_whitespace());
                // println!("{}", seq);
                self.sequence = seq;
                self.char_seq = self.sequence.chars().collect();
            },
            Err(msg) => {
                println!("{}", msg);
            }
        }
    }
    pub fn get_base(&self, pos: u32) -> char{
        self.char_seq[(pos - self.start) as usize]
    }
    pub fn get_length(&self, chr: &String) -> usize {
        let fai = fasta::read_index(&self.path);
        fai.get_length(chr)
    }
}
