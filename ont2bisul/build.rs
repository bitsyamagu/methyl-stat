extern crate cc;

fn main(){
    // println!("cargo:rustc-flags=-L/usr/local/hdf5/lib/libhdf5.a");
    println!("cargo:rustc-link-search=native=/usr/local/lib");
    println!("cargo:rustc-link-search=native=/usr/local/hdf5/lib");
    // println!("cargo:rustc-link-search=/ldisk1/yamagu/rust/fast5/rust_cpp/target/debug/build/fast5-rs-f676715db272de7c/out");
    cc::Build::new()
        .cpp(true)
        .warnings(true)
        .flag("-std=c++11")
        .flag("-Wall")
        .flag("-Wextra")
        .flag("-v")
        .flag("-g")
        .file("src/cpp/src/fast5.cpp")
        .include("src/cpp/include")
        .include("/usr/include/boost169")
        .include("/usr/local/hdf5/include")
        .compile("libfast5.a");
    println!("cargo:rustc-flags=-lhdf5");
}
