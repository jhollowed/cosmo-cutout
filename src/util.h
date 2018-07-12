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

// Hacc Types
#define POSVEL_T float
#define ID_T int64_t

using namespace std;

//////////////////////////////////////////////////////
//
//               reading functions
//
//////////////////////////////////////////////////////

void readHaloFile(string haloFileName, vector<float> &haloPos, vector<string> &haloTags);

int getLCSubdirs(string dir, vector<string> &subdirs);

int getLCFile(string dir, string &file);

int getLCSteps(int maxStep, int minStep, string dir, vector<string> &step_strings);


//////////////////////////////////////////////////////
//
//               cosmo Functions
//
//////////////////////////////////////////////////////

float aToZ(float a);

float zToStep(float z, int totSteps=500, float maxZ=200.0);


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

void rotation_matrix(int rank, const vector<vector<float> > &K, 
                     const float B, 
                     vector<vector<float> > &R);

#endif
