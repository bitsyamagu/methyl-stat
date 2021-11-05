#pragma once

#ifdef __cplusplus
extern "C" {
#endif

struct modbase {
    char* bases;
    uint32_t  width;
    uint32_t  height;
};

  void fast5_open(const char*);
  const char* read_fastq(const char*);
  struct modbase  read_modbase(const char*);
  void free_fastq(char const*);
  int get_data_index(char const*);
  void free_modbase( char const*);
  void fast5_close();

#ifdef __cplusplus
}
#endif
