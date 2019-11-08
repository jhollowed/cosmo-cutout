#ifndef UTIL_H_INCLUDED
#define UTIL_H_INCLUDED

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <iomanip>
#include <string>
#include <string.h>
#include <stdexcept>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <numeric>
#include <algorithm>
#include <sys/stat.h>
#include <dirent.h>
#include <errno.h>
#include <vector>
#include <omp.h>
#include <mpi.h>

using namespace std;

// Hacc Types
#define POSVEL_T float
#define ID_T int64_t

//////////////////////////////////////////////////////

struct Buffers_read {

    // Buffers to fill with data read from input LC, in processLC.cpp
    vector<POSVEL_T> x;
    vector<POSVEL_T> y;
    vector<POSVEL_T> z;
    vector<POSVEL_T> vx;
    vector<POSVEL_T> vy;
    vector<POSVEL_T> vz;
    vector<POSVEL_T> a;
    vector<ID_T> id;
    vector<int> rotation;
    vector<int32_t> replication;
    vector<float> theta;
    vector<float> phi;
};

struct Buffers_write {

    // Buffers to fill with data to write to cut out in processLC.cpp
    vector<POSVEL_T> x;
    vector<POSVEL_T> y;
    vector<POSVEL_T> z;
    vector<POSVEL_T> vx;
    vector<POSVEL_T> vy;
    vector<POSVEL_T> vz;
    vector<POSVEL_T> redshift;
    vector<ID_T> id;
    vector<int> rotation;
    vector<int32_t> replication;
    vector<POSVEL_T> theta;
    vector<POSVEL_T> phi;
    
    // Buffers to fill with MPI file writing offset values
    vector<int> np_rough_count; // length of particle vectors after rough cut (pre-rotation) for each rank
    vector<int> np_count; // length of output data vecotrs after refined cut (post-rotation) for each rank
    vector<int> np_offset; // cumulative sum of np_count
};

struct particle_pos {

    // struct for containing individual "primary" particle quantities
    POSVEL_T x;
    POSVEL_T y;
    POSVEL_T z;
    POSVEL_T theta;
    POSVEL_T phi;
    POSVEL_T a;
    ID_T id;
    int myrank;
};

struct particle_vel {

    // struct for containing individual "secondary" particle quantities
    POSVEL_T vx;
    POSVEL_T vy;
    POSVEL_T vz;
    int rotation;
    int32_t replication;
    int myrank;
};


//======================================================================================


//////////////////////////////////////////////////////
//
//               reading functions
//
//////////////////////////////////////////////////////

MPI_Datatype createParticles_pos();
MPI_Datatype createParticles_vel();

template<typename T>
bool comp_rank(const T &a, const T &b){
    // compare the myrank fields of two structs. Type T should be either a
    // particle_pos or particle_vel
    //
    // Params:
    // :param a: a particle struct, as defined in util.h
    // :param b: a particle struct, as defined in util.h
    // :return: a bool indicating whether or not the identifier of the rank
    //          possessing a is smaller in value than the identifier of the
    //          rank posessing b

    return a.myrank < b.myrank;
}

void resize_read_buffers(Buffers_read &r, int size, bool positionOnly, int extraSpace=0);

void comp_rank_scatter(size_t Np, vector<int> &idxRemap, int numranks);

bool does_file_exist(string filename);

void readHaloFile(string haloFileName, vector<float> &haloPos,
                  vector<string> &haloTags, vector<float> &haloProps,
                  string massDef = "sod");

int getLCSubdirs(string dir, vector<string> &subdirs);

int getLCFile(string dir, string &file);

int getLCSteps(int maxStep, int minStep, string dir, vector<string> &step_strings);

int prepStepSubdir(string subdir_name, bool overwrite, bool print, bool verbose);


//////////////////////////////////////////////////////
//
//               cosmo Functions
//
//////////////////////////////////////////////////////

float aToZ(float a);

int zToStep(float z, int totSteps=499, float maxZ=200.0);

float stepToZ(float step, int totSteps=499, float maxZ=200.0);


//////////////////////////////////////////////////////
//
//           matrix/vector operations
//
//////////////////////////////////////////////////////

void sizeMismatch();

vector<vector<float> > scalarMultiply(const vector<vector<float> > &matrix, 
                                      float scalar);

vector<vector<float> > squareMat(const vector<vector<float> > &matrix);

vector<float> matVecMul(const vector<vector<float> > &matrix, 
                        const vector<float> &vec);

float vecPairAngle(const vector<float> &v1,
                   const vector<float> &v2);

float dot(const vector<float> &v1, 
          const vector<float> &v2);

void cross(const vector<float> &v1, 
           const vector<float> &v2,
           vector<float> &v1xv2);

void normCross(const vector<float> &a,
           const vector<float> &b,
           vector<float> &k);

double determinant_3x3(const vector<vector<float> > &m);

vector<vector<float> > scale_adjoint_3x3(const vector<vector<float> > &m, float s = 1.0);
    
vector<vector<float> > invert_3x3(const vector<vector<float> > &m);


//////////////////////////////////////////////////////
//
//            coord rotation functions
//
//////////////////////////////////////////////////////

void cross_prod_matrix(const vector<float> &k, 
                       vector<vector<float> > &K);

void rotation_matrix(const vector<vector<float> > &K, const float B, 
                     vector<vector<float> > &R);

#endif
