This system is developped and tested in CentOS7 and CentOS8 environment.

# methyl-stat

'methyl-stat' is a collection of command-line tools for analyzing methylation aware ONT sequencing data.

- methylstat
  - methylstat mapps modification scores to reference bases. 
- methylcall
  - methylcall calls methylation for each base from results of methylstat.
- methylblock
  - methylblock detects highly methylated regions from results of methylcall.
- f5_to_fq
  - fastq extraction tool for reads with ModbaseProbs data.
- ont2bisul
  - ont2bisul transforms an ONT's BAM file to a virtually Bisul-converted BAM file to utilize IGV for detailed visualization.

# Installation
## Dependencies
methylstat reuqires following softwares and libraries:
- hdflib
   - https://www.hdfgroup.org/downloads/hdf5/source-code/ 
- hdf5 package (for h5ls commmand)
   - dnf install hdf5
- gcc (We tested with GNU C++-8.5.0)
- rust 1.39 or more
- jdk11 
- perl
- nextflow (for nextflow pipeline script)
- gatk-package-4.1.4.1-spark.jar or later
- htsjdk-2.17.0.jar or later
- commons-math3-3.6.1.jar (not 3.7)
## Building methylstat
```
cd methylstat
cargo build
```
## Building utilities
```
cd methylstat-util
make
```

## Building ont2bisul
- Build and install hdf5 libs
```
dnf install boost boost-level boost-filesystem

wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.1/src/hdf5-1.12.1.tar

cd hdf5-1.12.1
./configure --prefix=/usr/local/hdf5
make
make install
```
- Build ont2bisul
```
# dependency
dnf install perl-IPC-Cmd

git clone git@github.com:bitsyamagu/methyl-stat.git
cd methyl-stat/ont2bisul
cargo build
```
- Install
```
cp target/debug/ont2bisul /usr/local/bin
```

# Analysis

## Input
We are using long-read sequences produced from PromethION and base-called by guppy.

1. Map reads to reference genome by minimap2.
1. Merge all reads into a bam file, and sort and make index of it.
1. Create FAST5-to-readId index
  - perl fast5index.pl workspace > fast5_index.txt
1. (to be continued)

## ont2bisul usage
Example:
```
ont2bisul --region chr1:12000000-125000000 --bam my_test.bam \
     --fast5-dir workspace/ --index fast5index.txt \
     --fasta-index /path/to/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fai \
     --out chr_1_12000000-125000000.BS.bam
```
Options:
- --region:
  - For target region to convert
- --bam:
  - Your source BAM file. 
- --fast5-dir
  - The directory containing fast5 files
- --index:
  - The index file of read names in fast5 files. You can generate this by using fast5index.pl 
- --fasta-index  
  - Path to fai file of reference sequence
- --out:
  -  Bisul converted BAM file
