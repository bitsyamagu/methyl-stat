#include <iostream>
#include <string>
#include <boost/filesystem.hpp>
#include <cstring>
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
int get_data_index(const char* readid){
    hid_t grpid;
    hsize_t num_obj;
    hsize_t num_u_ana;

	const char* ptn = "/read_%s/Analyses";
	char fqpath[256];
	sprintf(fqpath, ptn, readid);
    grpid = H5Gopen(file, fqpath, H5P_DEFAULT);
    H5Gget_num_objs(grpid, &num_u_ana);
    char name[256];
    char dpath[256];

    const char* pattern = "Basecall_1D_";
    int check_len = strlen(pattern);
    int i = 0;
    for(i = 0; i<num_u_ana; i++){
        H5Gget_objname_by_idx(grpid, i, name, 255);
        if(strncmp(name, pattern, check_len) != 0){
            break;
        }
    }
	H5Gclose(grpid);
    return i-1;
}
const char* read_fastq(const char* path){
    // boost::filesystem::create_directory("dirname");
    // std::cerr << path << std::endl;
    hid_t dset = H5Dopen (file, path, H5P_DEFAULT);
    hid_t filetype = H5Dget_type (dset);
    hsize_t     dims[1] = {DIM0};

    size_t sdim = H5Tget_size (filetype);
    sdim++; // for NULL terminator

    hid_t space = H5Dget_space (dset);
    // int ndims = H5Sget_simple_extent_dims (space, dims, NULL);

    char** rdata = (char **) malloc (dims[0] * sizeof (char *));
    rdata[0] = (char *) malloc (dims[0] * sdim * sizeof (char));

    for (hsize_t i=1; i<dims[0]; i++)
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
    // free (rdata);
    status = H5Dclose (dset);
    status = H5Sclose (space);
    status = H5Tclose (filetype);
    status = H5Tclose (memtype);

    //cerr << "size = " << ( dims[0] * sdim * sizeof (char)) << endl;
    return result;
}

struct modbase read_modbase(const char* path){
    // boost::filesystem::create_directory("dirname");
    // std::cerr << path << std::endl;
    hid_t dset = H5Dopen (file, path, H5P_DEFAULT);
    hid_t filetype = H5Dget_type (dset);
    hsize_t     dims[1] = {DIM0};

    size_t sdim = H5Tget_size (filetype);
    // sdim++; // for NULL terminator

    hid_t space = H5Dget_space (dset);
    int ndims = H5Sget_simple_extent_dims (space, dims, NULL);

    char** rdata = (char **) malloc (dims[0] * sizeof (char *));
    rdata[0] = (char *) malloc (dims[0] * dims[1] *sdim * sizeof (char));
    sprintf(rdata[0], "sdim %d, ndims %d dims[0] %d dims[1] %d\n", sdim, ndims, dims[0], dims[1]);
    // sdim 2, ndims 2 dims[0] 7154 dims[1] 6
    // printf("%s", rdata[0]);

    for (int i=1; i<dims[0]; i++)
         rdata[i] = rdata[0] + i * dims[1];

    hid_t memtype = H5Tcopy (H5T_STD_U8LE);
    herr_t status = H5Tset_size (memtype, sdim);

    status = H5Dread (dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata[0]);

    if(false){
    // if(strstr(path, "read_77ab8eeb-fa02-4af0-abe8-1e8016e3da96") != NULL){
        for(int j = 0; j<dims[0]; j++){
        // for(int j = 0; j<4; j++){
            unsigned char line[6];
            for(int idx = 0; idx<6; idx++){
                line[idx] = (unsigned char)rdata[j][idx];
            }
            // printf("# %s %hu %hu %hu %hu %hu %hu\n", path, line[0], line[1], line[2], line[3], line[4], line[5]);
        }
    }
    struct modbase result = {rdata[0], dims[0], dims[1]};

    // for (int i=0; i<dims[0]; i++)
    //     printf ("[%d]: %d\n", i, rdata[i]);
    // free (rdata[0]);
    // free (rdata);
    status = H5Dclose (dset);
    status = H5Sclose (space);
    status = H5Tclose (filetype);
    status = H5Tclose (memtype);

    //cerr << "size = " << ( dims[0] * sdim * sizeof (char)) << endl;
    return result;
}
void free_fastq(char const* mem){
    free(const_cast<char*>(mem));
}
void free_modbase(char const* mem){
    free(const_cast<char*>(mem));
}

void fast5_close(){
    // std::cerr << "close()" << std::endl;
    H5Fclose (file);
}
