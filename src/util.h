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

void readHaloFile(string haloFileName, vector<float> &haloPos, vector<ID_T> &haloTags);

int getLCSubdirs(string dir, vector<string> &subdirs);

int getLCFile(string dir, string &file);

int getLCSteps(int maxStep, int minStep, string dir, vector<string> &step_strings);

//////////////////////////////////////////////////////
//
//               helper Functions
//
//////////////////////////////////////////////////////

float aToZ(float a);

float zToStep(float z, int totSteps=499, float maxZ=200.0);

//////////////////////////////////////////////////////
//
//            coord rotation functions
//
//////////////////////////////////////////////////////

void sizeMismatch();

float vecPairAngle(const vector<float> &v1,
                   const vector<float> &v2);

void cross(const vector<float> &v1, 
           const vector<float> &v2,
           vector<float> &v1xv2);

void normCross(const vector<float> &a,
           const vector<float> &b,
           vector<float> &k);

void rotate(const vector<float> &k_vec,
            const float B, 
            const vector<float> &v_vec, 
            vector<float> &v_rot); 
#endif
