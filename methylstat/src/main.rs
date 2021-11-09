#[warn(unused_imports)]
extern crate rust_htslib;
extern crate bio_types;
extern crate libc;
// extern crate hdf5;

use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::Record;
use std::io::{self, BufRead, Write};
use std::fs::File;
use std::collections::HashMap;
use std::slice;
use bio_types::strand::ReqStrand;
use libc::{ c_void, c_char, size_t };
use std::ffi::{ CStr, CString };

#[macro_use]
mod fastq;
mod fasta;
mod genomic_region;

use genomic_region::GenomicRegion;

#[repr(C)]
pub struct ModBase {
    pub bases: *const c_char,
    pub width: u32,
    pub height: u32
}

#[link(name = "fast5", kind = "static")]
extern "C" {
  fn fast5_open(s: *const c_char);
  fn read_fastq(s: *const c_char) -> *const c_char;
  fn read_modbase(s: *const c_char) -> ModBase;
  fn get_data_index(s: *const c_char) -> i32;
  fn free_fastq(s: *const c_char);
  // fn free_modbase(s: *const c_char);
  fn fast5_close();
}
// https://qiita.com/taichitk/items/d72bcc32345caa08ce8a

pub struct HaplotypeProfile {
    pub left: u32,
    pub right: u32,
    pub discrepancy: u32
}

pub struct GenomicInterval {
    pub chr: String,
    pub start: u32,
    pub end: u32
}

fn collect_reads(bam: &mut bam::IndexedReader, region: &String) -> Vec<Record> {
    // println!("{}", region);
    match bam.fetch(region.as_bytes()) {
        Ok(_n) => { /* noop */ },
        Err(e) => {
                println!("{:?}\n---------------------------------\n", e);
                panic!("{:?}", e);
        }
    }
    let mut result = Vec::<Record>::new();
    for r in bam.records() {
        let record = r.ok().expect("Error reading BAM file.");
        result.push(record);
    }
    result
}

fn parse_interval(interval: &mut GenomicInterval, region_str: &str) {
    let buf: Vec<&str> = region_str.split(':').collect();
    interval.chr = buf[0].to_string();
    // println!("{:?}", buf);
    let buf2: Vec<&str> = buf[1].split('-').collect();
    interval.start = buf2[0].to_string().parse::<u32>().unwrap();
    interval.end = buf2[1].parse::<u32>().unwrap();
}
fn main() {
  // let interval = 100000;
  let argv: Vec<String> = std::env::args().collect();
    // https://stackoverflow.com/questions/24145823/how-do-i-convert-a-c-string-into-a-rust-string-and-back-via-ffi
    // let mut bam = bam::IndexedReader::from_path(&argv[2]).unwrap();

    let mut fai_path = String::from("");
    let mut bam = bam::IndexedReader::from_path(
       "/disk1/methyl/merged2/merged2.bam");

    let mut f5map: HashMap<String, String> = HashMap::new();
    let mut fast5_dir = "/nas2/methyl/workspace";
    let mut interval = GenomicInterval{ chr: String::from(""), start: 0, end: 0};

    for i in  1..argv.len() {
        if argv[i].as_str() == "--region" {
            parse_interval(&mut interval, argv[i+1].as_str());
        } else if argv[i].as_str() == "--index" {
            let index_file_name = argv[i+1].as_str();
            let index_file = File::open(index_file_name).expect("Couldn't open index");
            for raw_line in io::BufReader::new(index_file).lines() {
                if let Ok(line) = raw_line {
                    let buf: Vec<&str> = line.split('\t').collect();
                    // println!("# {}", buf[0]);
                    f5map.insert(buf[0].to_string(), buf[1].to_string());
                }
            }
        } else if argv[i].as_str() == "--fast5-dir" {
            fast5_dir = argv[i+1].as_str();
        } else if argv[i].as_str() == "--fasta-index" {
            fai_path = String::from(argv[i+1].to_string());
        } else if argv[i].as_str() == "--bam" {
            bam = bam::IndexedReader::from_path(argv[i+1].to_string());
        }
    }
	if interval.end == 0 {
    	panic!("Use --region option");
	}
    let mut file_out = File::create(format!("{}_{}-{}.txt", interval.chr, interval.start, interval.end)).unwrap();
    // testing on KCNQ1OT1 ************
    // let start = 2629557;
    // let end = 2721228;
        
    let chr =  interval.chr;
    let start = interval.start;
    let end = interval.end;
    // println!("{}:{}-{}", chr, start, end);
    let reads = collect_reads(&mut bam.unwrap(), &format!("{}:{}-{}", chr, start, end));
    // let reads = collect_reads(&mut bam, &format!("{}:{}-{}", chr, start, end));
    // let vcf = phased_vcf::load_phased_vcf(&String::from(format!("../merged2/phased.chr{}.vcf", chr)), &chr, start, end);
    
    // println!("genomic region: {} {} {}", chr, start-1, end);
    let mut genomic_region = GenomicRegion::new(&chr, start-1, end, &fai_path.to_string());
    genomic_region.init();
    // let chr_length = genomic_region.get_length(&chr) as u32;
    
    // println!("loop for: {}", reads.len());
    
    for mut read in reads {
        // println!("{}", String::from_utf8(read.qname().to_vec()).unwrap());
        // let fname = CString::new(argv[1].as_str()).unwrap();
    
        let mut fastq;
		let mut grpidx = -1;
        if let Some(fq) = get_fastq(&read,  &f5map, fast5_dir.to_string(), &mut grpidx) {
            fastq = fastq::parse_fastq(&fq);
            // println!("name: {}", fastq.name);
            // println!("seq: {}", fastq.seq);
        }else {
            continue;
        }
    
        let qname = String::from_utf8(read.qname().to_vec()).unwrap();
        let mpath = format!("/read_{}/Analyses/Basecall_1D_00{}/BaseCalled_template/ModBaseProbs", &qname, grpidx);
        // println!("{}", mpath);
        let modbase_path = CString::new(mpath).unwrap();
        let mb = unsafe { read_modbase(modbase_path.as_ptr()) };
        unsafe {
            fast5_close();
        }
        // dump_seq(&mb, &qname);
        add_mod_info(&mut fastq, &mb);
    
        // println!("{:?}", fastq.mod5mC);
        // check_haplotype(&fastq, &mut read, &vcf);
        dump_bed(&fastq, &mut read, &genomic_region, &mut file_out);
    }
}
pub fn dump_bed(fastq: &fastq::FastQ, read: &mut bam::Record, genome: &GenomicRegion, file: &mut std::fs::File){
    let mut start = read.pos();
    let mut end = read.cigar().end_pos();
    let seq: Vec<char> = fastq.seq.chars().collect();
    // println!("{}..{}\t{}", start, end, end-start);
    let mut strand = '-';
    if read.strand() == ReqStrand::Forward {
        strand = '+';
    }
    if start < genome.start as i64 {
        start = genome.start as i64;
    }
    if end > genome.end as i64 {
        end = genome.end as i64;
    }
    let mut last_pos_in_read: i64 = 0;
	let mod5mC = &fastq.mod5mC.as_ref().unwrap();

	let mut inserted;
	let mut inserted_bases;
	let mut inserted_len: usize = 0;

    for genomic_pos in start..end {
        // println!("genomic_pos: {} {} {}", genomic_pos, genome.start, genome.end);
        let genomic_base = genome.get_base(genomic_pos as u32);
		let mut direction = 1;
	    if strand == '-' {
			direction = -1;
			if !genome.can_read_base(genomic_pos + direction) {
				continue;
			}
		}
        let genomic_next_base = genome.get_base((genomic_pos+direction) as u32);
        if (((genomic_base != 'C' && genomic_base != 'c') || (genomic_next_base != 'g' && genomic_next_base != 'G')) && strand == '+') 
				|| (((genomic_base != 'G' && genomic_base != 'g') || (genomic_next_base != 'c' && genomic_next_base != 'C')) &&  strand == '-')  {
            continue;
        }
        let read_pos: Option<u32> = read.cigar().read_pos(genomic_pos as u32, true, true).unwrap();
        match read_pos {
            Some(read_pos) => {
                let pos_in_read: i64 = ref_to_read(fastq.seq.len(), strand, read_pos);
				if !(pos_in_read+1 < fastq.seq.len() as i64) {
					continue;
				}
                let base: char = seq[pos_in_read as usize];
				let nextbase: char = seq[(pos_in_read + 1) as usize];
				if nextbase != 'g' && nextbase != 'G' {
					continue;
				}
				
				if last_pos_in_read != pos_in_read && last_pos_in_read != 0 {
					inserted_len = (i64::abs(pos_in_read - last_pos_in_read) - 1) as usize;
					inserted = vec![0 as u8; inserted_len];
					inserted_bases = vec!['N' as char; inserted_len];
					for p in 0..inserted_len {
						let mut c: i64 = 1;
						if strand == '-' {
							c = -1;
						}
						inserted[p] = mod5mC[((last_pos_in_read+c) as i64 + c*(p as i64)) as usize];
						inserted_bases[p] = seq[((last_pos_in_read+c) as i64 + c*(p as i64)) as usize];
					}
				}else {
					inserted = vec![0 as u8; 0];
					inserted_bases = vec!['n' as char; 0];
				}

                writeln!(file, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:?}\t{:?}", genomic_pos, &fastq.name, mod5mC[pos_in_read as usize], genomic_base, base, strand, pos_in_read, inserted_len, &inserted, &inserted_bases);
				last_pos_in_read = pos_in_read;
            },
            None => {
                // panic?
                println!("no read_pos for {}", genomic_pos);
            }
        }
    }
}
pub fn ref_to_read(read_len: usize, strand: char, read_pos: u32) -> i64 {
	let pos: i64;
	match strand {
		'-' ... '-' => {
			pos = (read_len - 1 - read_pos as usize) as i64;
			// println("read_pos={} read_len={}", read_pos, read_len);
		},
		'+' ... '+' => {
			pos = read_pos as i64;
		},
		_ => { pos = 0; }
	}
	pos
}
pub fn complementary(b: char) -> char{
    let val: char;
    match b {
        'C' ... 'C'=> { val = 'G'; },
        'G' ... 'G'=> { val = 'C'; },
        'A' ... 'A'=> { val = 'T'; },
        'T' ... 'T'=> { val = 'A'; },
        'c' ... 'c'=> { val = 'f'; },
        'g' ... 'g'=> { val = 'c'; },
        'a' ... 'a'=> { val = 't'; },
        't' ... 't'=> { val = 'a'; },
        _ => { val = 'N'}
    }
    val
}
        
/**
 *  /// For a given position in the reference, get corresponding position within read.
    /// If reference position is outside of the read alignment, return None.
    ///
    /// # Arguments
    ///
    /// * `ref_pos` - the reference position
    /// * `include_softclips` - if true, softclips will be considered as matches or mismatches
    /// * `include_dels` - if true, positions within deletions will be considered (first reference matching read position after deletion will be returned)
    ///
    pub fn read_pos(
        &self,
        ref_pos: u32,
        include_softclips: bool,
        include_dels: bool,
    ) -> Result<Option<u32>> {
*/
/*
pub fn check_haplotype(fastq: &fastq::FastQ, read: &mut bam::Record, vcf: &Vec<phased_vcf::SNP>) -> Option<HaplotypeProfile>{
    let mut profile = None;
    let start_index = vcf.iter().position(|x| x.pos > read.pos() as u32).unwrap();
    let slice = &vcf[start_index..];
    let mut left = 0;
    let mut right = 0;
    let mut discrepancy = 0;
    for v in slice {
        let rpos = read.cigar().read_pos(v.pos-1, true, true).unwrap();
        match rpos {
            Some(rpos) => {
                /* println!("REV: pos: {} ref:{} alt:{}, read:{}, genotype: {}, alleles: {},{}", 
                 *            v.pos, v.ref_allele, v.alt_alleles, read_base, v.genotype, left_allele, right_allele); */
                let read_base: char;
                if read.strand() == ReqStrand::Reverse {
                    read_base = complementary(fastq.seq.chars().nth(fastq.seq.len() - 1 - rpos as usize).unwrap());
                }else {
                    read_base = fastq.seq.chars().nth(rpos as usize).unwrap();
                }
                let left_allele = &v.all_alleles[v.genotype[0..1].parse::<usize>().unwrap() as usize];
                let right_allele = &v.all_alleles[v.genotype[2..3].parse::<usize>().unwrap() as usize];
                let not_right: bool = read_base != right_allele.chars().nth(0).unwrap();
                let not_left: bool = read_base != left_allele.chars().nth(0).unwrap();
                if not_left {
                    right += 1;
                }
                if not_right {
                    left += 1;
                }
                if not_left && not_right {
                    discrepancy += 1;
                    left -= 1;
                    right -= 1;
                    // print!("*[{}]", fastq.name);
                }
            },
            None => { // end of the read
                break;
            }
        }
        // println!("left: {}, right: {}, discrepancy: {}", left, right, discrepancy);
        profile = Some(HaplotypeProfile{ left: left, right: right, discrepancy: discrepancy });
    }
    profile
}
*/
pub fn get_fastq(read: &Record, f5map: &HashMap<String, String>, fast5_dir: String, grpidx: &mut i32) -> Option<String> {
    let qname: String = format!("read_{}", String::from_utf8(read.qname().to_vec()).unwrap());
    let fname: &String;
    match f5map.get(&qname) {
        Some(f) => { fname = f; },
        None => { /*println!("no for {}",  qname);*/ return None; }
    }
    // println!("{}", fname);
    let cfname = CString::new(format!("{}/{}", fast5_dir, fname)).unwrap();
    // let gpath = format!("/read_{}/Analyses/", String::from_utf8(read.qname().to_vec()).unwrap());
    unsafe{
        fast5_open(cfname.as_ptr());
		if(*grpidx < 0){
			*grpidx = get_data_index(CString::new(read.qname()).unwrap().as_ptr());
		}
    }
    let dpath = format!("/read_{}/Analyses/Basecall_1D_00{}/BaseCalled_template/Fastq", String::from_utf8(read.qname().to_vec()).unwrap(), grpidx);
    // println!("{}", dpath);
    let fastq_path = CString::new(dpath).unwrap();
    let c_buf: *const c_char = unsafe { read_fastq(fastq_path.as_ptr()) };
    let c_str: &CStr = unsafe { CStr::from_ptr(c_buf) };
    // print!("{}", c_str.to_str().unwrap());
    let result: String = String::from(c_str.to_str().unwrap());

    Some(result)
}
// [A, 6mA, C, 5mC, G, T]
#[allow(non_snake_case)]
pub fn add_mod_info(fastq: &mut fastq::FastQ, mb: &ModBase/*, qname: &String*/) {
    let mut idx: usize = 0;
    let array: &[c_char] = unsafe {slice::from_raw_parts(mb.bases, (mb.width*mb.height) as usize)};
    // println!("{}: (width, height) = ({}, {})", qname, mb.width, mb.height);
    let COL_6mA = 1; let COL_5mC = 3;
    let mut prob6mA = vec![0 as u8; mb.width as usize];
    let mut prob5mC = vec![0 as u8; mb.width as usize];
    for row in 0..mb.width {
        for col in 0..mb.height {
            // print!("{} ", array[idx] as u8);
            if col == COL_6mA {
                prob6mA[row as usize] = array[idx] as u8;
            }else if col == COL_5mC {
                prob5mC[row as usize] = array[idx] as u8;
            }
            idx += 1;
        }
        // println!("");
    }
    fastq.mod6mA = Some(prob6mA);
    fastq.mod5mC = Some(prob5mC);
}
pub fn dump_seq(mb: &ModBase, qname: &String) {
    let mut idx: usize = 0;
    let array: &[c_char] = unsafe {slice::from_raw_parts(mb.bases, (mb.width*mb.height) as usize)};
    // println!("{}: (width, height) = ({}, {})", qname, mb.width, mb.height);
    for _i in 0..mb.width {
        for _j in 0..mb.height {
            print!("{} ", array[idx] as u8);
            idx += 1;
        }
        println!("");
    }
}

