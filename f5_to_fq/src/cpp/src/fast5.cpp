#include <iostream>
#include <string>
#include <boost/filesystem.hpp>
#include "fast5.hpp"
#include "hdf5.h"

using namespace std;

hid_t file = 0L;

#define DIM0 1

// Example:
//    https://support.hdfgroup.org/ftp/HDF5/examples/examples-by-api/hdf5-examples/1_10/C/H5T/h5ex_t_string.c

void fast5_open(const char* path){
  // std::string spath = std::string(path);
  // std::cerr << path << std::endl;
 
  file = H5Fopen (path, H5F_ACC_RDONLY, H5P_DEFAULT);
}
/*
int get_last_basecall_index(){
	char* path_temp = "/Analyses/Basecall_1D_%03d/BaseCalled_template/Fastq";
	int last_index = 0;
	for(int i = 0; i<10; i++){
		char path[256];
		sprintf(path, path_temp, i);
		// hid_t h = H5Dopen(file, path, H5P_DEFAULT);
		//fprintf(stdout, "hid_t h: %ld\n", (long)h);
	}
	return 0;
}*/
const char* read_fastq(const char* path){
    // boost::filesystem::create_directory("dirname");
    // std::cerr << path << std::endl;
    hid_t dset = H5Dopen (file, path, H5P_DEFAULT);
    hid_t filetype = H5Dget_type (dset);
    hsize_t     dims[1] = {DIM0};

    size_t sdim = H5Tget_size (filetype);
    sdim++; // for NULL terminator

    hid_t space = H5Dget_space (dset);
    int ndims = ndims = H5Sget_simple_extent_dims (space, dims, NULL);

    char** rdata = (char **) malloc (dims[0] * sizeof (char *));
    rdata[0] = (char *) malloc (dims[0] * sdim * sizeof (char));

    for (int i=1; i<dims[0]; i++)
         rdata[i] = rdata[0] + i * sdim;

    hid_t memtype = H5Tcopy (H5T_C_S1);
    herr_t status = H5Tset_size (memtype, sdim);

    status = H5Dread (dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata[0]);

    char* result = rdata[0];

    /*
    for (int i=0; i<dims[0]; i++)
        printf ("[%d]: %s\n", i, rdata[i]);
        */
    // free (rdata[0]);
    free (rdata);
    status = H5Dclose (dset);
    status = H5Sclose (space);
    status = H5Tclose (filetype);
    status = H5Tclose (memtype);

    //cerr << "size = " << ( dims[0] * sdim * sizeof (char)) << endl;
    return result;
}

void fast5_close(){
    std::cerr << "close()" << std::endl;
    H5Fclose (file);
}
