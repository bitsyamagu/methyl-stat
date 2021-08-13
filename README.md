# methyl-stat

'methyl-stat' is a collection of command-line tools for analyzing methylation aware ONT sequencing data.

- methylstat
  - methylstat mapps modification scores to reference bases. 
- methylcall
  - methylcall calls methylation for each base from results of methylstat.
- methylblock
  - methylblock detects highly methylated regions from results of methylcall.
- ont2bisul
  - ont2bisul transforms an ONT's BAM file to a virtually Bisul-converted BAM file to utilize IGV for detailed visualization.
