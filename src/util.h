#ifndef UTIL_H_INCLUDED
#define UTIL_H_INCLUDED

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <string.h>
#include <stdexcept>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <sstream>
#include <omp.h>
#include <sys/stat.h>
#include <dirent.h>
#include <errno.h>
#include <vector>

// Generic IO
//#include "GenericIO.h"

using namespace std;
//using namespace gio;


//////////////////////////////////////////////////////
//
//         Helper Functions
//
//////////////////////////////////////////////////////

float redshift(float a);

float zToStep(float z, int totSteps=499, float maxZ=200.0);

//////////////////////////////////////////////////////
//
//         Coord Rotation functions
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

//////////////////////////////////////////////////////
//
//             Reading functions
//
//////////////////////////////////////////////////////

int getLCSubdirs(string dir, vector<string> &subdirs);

int getLCFile(string dir, string &file);

int getLCSteps(int minStep, string dir, vector<string> &step_strings);

#endif
