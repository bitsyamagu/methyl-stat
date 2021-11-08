#[allow(non_snake_case)]
pub struct FastQ {
    pub raw_name: String,
    pub name: String,
    pub seq: String,
    pub qual: String,
    pub mod6mA: Option<Vec<u8>>,
    pub mod5mC: Option<Vec<u8>>
}
pub fn parse_fastq(src: &String) -> FastQ {
    let buf: Vec<&str> = src.split('\n').collect();
    let buf2: Vec<&str> = buf[0].split(' ').collect();
    let fq = FastQ {raw_name: buf[0].to_string(), name: buf2[0].to_string(), seq: buf[1].to_string(), qual: buf[3].to_string(), mod6mA: None, mod5mC: None};
    fq
}
