#pragma once

#ifdef __cplusplus
extern "C" {
#endif

  void fast5_open(const char*);
  // int get_last_basecall_index();
  const char* read_fastq(const char*);
  void fast5_close();

#ifdef __cplusplus
}
#endif
