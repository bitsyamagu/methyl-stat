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
- gcc
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
# Analysis

## Input
We are using long-read sequences produced from PromethION and base-called by guppy.

1. Map reads to reference genome by minimap2.
1. Merge all reads into a bam file, and sort and make index of it.
1. Create FAST5-to-readId index
  - perl fast5index.pl workspace > fast5_index.txt
1. (to be continued)
