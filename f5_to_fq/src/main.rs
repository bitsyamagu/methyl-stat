extern crate libc;
extern crate hdf5;

use std::io::{self, BufRead};
use std::fs::File;


use libc::{ c_void, c_char, size_t };
use std::ffi::{ CStr, CString };

#[link(name = "fast5", kind = "static")]
extern "C" {
  // fn get_last_basecall_index();
  fn fast5_open(s: *const c_char);
  fn read_fastq(s: *const c_char) -> *const c_char;
  fn fast5_close();
}

fn main() {
    // https://stackoverflow.com/questions/24145823/how-do-i-convert-a-c-string-into-a-rust-string-and-back-via-ffi

    let argv: Vec<String> = std::env::args().collect();
    let path = &argv[1];
    let file;
    match  hdf5::File::open(path) {
      Err(e) => {
          println!("{:?}", e);
          panic!("Error");
      },
      Ok(f) => {
          file = f;
      },
    }
    let root = file.group("/").unwrap();
    let read_ids: Vec<String> = root.member_names().unwrap();
    let mut last_gname:String = String::new();
    for read_id in read_ids {
	  let basecalls: Vec<String> = file.group(format!("/{}/Analyses", read_id).as_str()).unwrap().member_names().unwrap();
	  for gname in basecalls {
		  if(gname.starts_with("Basecall_1D_")){
		  	last_gname = gname.clone();
	  		// println!("last: /{}/Analyses/{}", read_id, last_gname);
		  }
	  }
	  // println!("/{}/Analyses/{}/BaseCalled_template/Fastq", read_id, last_gname);
	  let dpath = format!("/{}/Analyses/{}/BaseCalled_template/Fastq", read_id, last_gname);
	  let fastq_path = CString::new(dpath).unwrap();
	  let fname = CString::new(argv[1].as_str()).unwrap();
	  unsafe{
    	  fast5_open(fname.as_ptr());
	  }
 
	  let c_buf: *const c_char = unsafe { read_fastq(fastq_path.as_ptr()) };
	  let c_str: &CStr = unsafe { CStr::from_ptr(c_buf) };

	  print!("{}", c_str.to_str().unwrap());
	}
  unsafe {
    fast5_close();
  }
}
