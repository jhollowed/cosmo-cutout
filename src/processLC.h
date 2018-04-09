#ifndef PROCESSLC_H_INCLUDED
#define PROCESSLC_H_INCLUDED

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

#include "util.h"

#define PI 3.14159265
#define ARCSEC 3600.0

// Generic IO
#include "GenericIO.h"

// Cosmotools
#define REAL double

using namespace std;
using namespace gio;

void processLC(string dir_name, string out_dir, vector<string> step_strings, 
               vector<float> theta_bounds, vector<float> phi_bounds);

void processLC(string dir_name, string out_dir, vector<string> step_strings, 
               vector<float> halo_pos, float boxLength);

#endif
