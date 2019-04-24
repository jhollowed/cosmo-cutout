#include "util.h"

using namespace std;


//////////////////////////////////////////////////////
//
//             reading functions
//
//////////////////////////////////////////////////////


MPI_Datatype createParticles_pos(){
    // This function creates and returns an MPI struct type which has a field per 
    // "primary" particle quantity (position, redshift, id and rank). 
    // One particle shall be represented by one "particles_mpi" object, which is 
    // based upon the "particle_pos" struct in util.h
    //
    // Params:
    // :return: a struct of custom MPI type "particles_mpi"

    MPI_Datatype particles_mpi;
    MPI_Datatype type[9] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT,
                            MPI_FLOAT, MPI_FLOAT, MPI_INT64_T, MPI_INT};
    int blocklen[9] = {1,1,1,1,1,1,1,1,1};
    MPI_Aint disp[9] = {
                         offsetof(particle_pos, x),
                         offsetof(particle_pos, y),
                         offsetof(particle_pos, z),
                         offsetof(particle_pos, d),
                         offsetof(particle_pos, theta),
                         offsetof(particle_pos, phi),
                         offsetof(particle_pos, a),
                         offsetof(particle_pos, id),
                         offsetof(particle_pos, myrank)
                        };
    MPI_Type_struct(9, blocklen, disp, type, &particles_mpi);
    MPI_Type_commit(&particles_mpi);
    return particles_mpi;
}


MPI_Datatype createParticles_vel(){
    // This function creates and returns an MPI struct type which has a field per 
    // "secondary" particle quantity (velocity, rotation and replication info). 
    // One particle shall be represented by one "particles_mpi" object, which is 
    // based upon the "particle_vel" struct in util.h
    //
    // Params:
    // :return: a struct of custom MPI type "particles_mpi"

    MPI_Datatype particles_mpi;
    MPI_Datatype type[6] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, 
                             MPI_INT, MPI_INT32_T, MPI_INT};
    int blocklen[6] = {1,1,1,1,1,1};
    MPI_Aint disp[6] = {
                         offsetof(particle_vel, vx),
                         offsetof(particle_vel, vy),
                         offsetof(particle_vel, vz),
                         offsetof(particle_vel, rotation),
                         offsetof(particle_vel, replication),
                         offsetof(particle_vel, myrank)
                        };
    MPI_Type_struct(6, blocklen, disp, type, &particles_mpi);
    MPI_Type_commit(&particles_mpi);
    return particles_mpi;
}


//======================================================================================


void comp_rank_scatter(size_t Np, vector<int> &idxRemap, int numranks){
    // Constructs a vector divided into numranks chunks sharing a common identifier
    //
    // Params:
    // :param Np: desired length of the constructed vector
    // :param idxRemap: vector within which to store result
    // :param numranks: number of unique populations resultant in idxRemap
    // :return: None

    float np_per_rank = float(Np) / numranks;
    for(int j = 0; j < Np; ++j){
        idxRemap.push_back(int(j / np_per_rank));
    }
}


//======================================================================================


bool comp_by_theta(const particle_pos &a, const particle_pos &b){
    // Compares two particle_pos structs by their 'theta' field
    //
    // Params:
    // :param a: a particle_pos struct
    // :param b: a particle_pos struct
    // :return: 1 if a.theta < b.theta, 0 otherwise
    
    return a.theta < b.theta;
}


//======================================================================================


bool does_file_exist(string filename){
    // Checks if a file exists     //
    // Params:
    // :param filename: a filename as a string
    // :return: a boolean value (true if file exists, false if not, 
    //          or false if file is locked/inaccessible)

    std::ifstream f(filename);
    return f.good();
}


//======================================================================================


void readHaloFile(string haloFileName, vector<float> &haloPos, 
                                       vector<string> &haloTags, 
                                       vector<float> &haloProps,
                                       string massDef){
    // This function reads halo identifiers and positions from an input text file, where 
    // that file is expected to have one halo per row, each space-delimited row appearing as 
    //
    // tag1 redshift1 shell1 mass1 radius1 concentration1 concentration_err1 x1 y1 z1
    // tag2 redshift2 shell2 mass2 radius2 concentration2 concentration_err2 x2 y2 z2
    // tag3 redshift3 shell3 mass3 radius3 concentration3 concentration_err3 x3 y3 z3
    // ...
    //
    // and parses them into three vectors to contain the ids and positions, as such:
    //
    // vector  haloTags: {tag1, tag2,...}
    // vector haloPos: {x1, y1, z1, x2, y2, z2,...}
    // vector haloProps: {redshift1, shell1, mass1, radius1, conc1, conc_err1, redshift2, shell2, mass2, radius2,...}
    //
    // - The positions are expected to be able to be cast to floats
    // - The tags are expected to be able to be cast as strings
    // - The redshifts, masses, radii, and concentrations  are expected to be able to be cast as floats.
    // - The shell is expected to be a snapshot timestep number that can be cast as an int
    // - The radius and concentration can be omitted if the argument massDef=='fof'. 
    //   If those quantities are not present in the input file, and the massDef arg is not set properly, 
    //   or vice-verse, everything will break!
    //  
    // A code module to generate such a text file for a given HACC simulation can be cloned at  
    // https://github.com/jhollowed/cosmology/blob/master/lightcone/list_lc_halos.py
    //
    // Params:
    // :param haloFileName: the ascii file containing the halo data, in the format 
    //                      described above
    // :param haloPos: float vector to hold halo positions
    // :param haloTags: string vector to hold halo tags
    // :param haloProps: float vector to hold halo info including mass, redshift, and tag if massDef="fof" 
    //                   and also radius, concentration, and concentration error if massDef="sod"
    // :param massDef: the mass definition given in the ids of the haloFile, either 
    //                 'sod', or 'fof'. Defaults to 'sod'
    
    ifstream haloFile(haloFileName.c_str());
    copy(istream_iterator<string>(haloFile), 
         istream_iterator<string>(), 
         back_inserter(haloTags));

    int rowLen;
    if(massDef == "fof")
        rowLen = 7;
    else if(massDef == "sod")
        rowLen = 10;
    else{
        cout << "Unknown mass definition " << massDef << 
                ". Valid options are \'fof\' or \'sod\'" << endl;
        MPI_Abort(MPI_COMM_WORLD, 0);
    }
    
    // ensure input file is as expected, more or less 
    if(haloTags.size() % rowLen != 0){ 
        if(massDef == "sod")
            cout << "\nWhen massDef==\'sod\', each halo position given in input file must " <<
                    "have 10 quantities in the space-delimited form: " <<
                    "id redshift m200c r200c cdelta cdelta_err x y z " << endl;
        else
            cout << "\nWhen massDef==\'fof\', each halo position given in input file must " <<
                    "have 7 quantities in the space-delimited form: " <<
                    "id redshift m_fof x y z " << endl;
        MPI_Abort(MPI_COMM_WORLD, 0);
    } 
    
    for(int i=0; i<haloTags.size()/rowLen; ++i){
        
        // tag is every 7th or 10th element
        haloTags[i] = haloTags[i*rowLen];
        // redshift, shell, mass, and radius
        haloProps.push_back( strtof(haloTags[i*rowLen+1].c_str(), 0) );
        haloProps.push_back( strtof(haloTags[i*rowLen+2].c_str(), 0) );
        haloProps.push_back( strtof(haloTags[i*rowLen+3].c_str(), 0) );
        if(massDef == "sod")
            haloProps.push_back( strtof(haloTags[i*rowLen+4].c_str(), 0) );
            haloProps.push_back( strtof(haloTags[i*rowLen+5].c_str(), 0) );
            haloProps.push_back( strtof(haloTags[i*rowLen+6].c_str(), 0) );
        // x, y, z
        haloPos.push_back( strtof(haloTags[i*rowLen+(rowLen-3)].c_str(), 0) );
        haloPos.push_back( strtof(haloTags[i*rowLen+(rowLen-2)].c_str(), 0) );
        haloPos.push_back( strtof(haloTags[i*rowLen+(rowLen-1)].c_str(), 0) );
    }
    haloTags.resize(haloTags.size()/rowLen);
}


//======================================================================================


int getLCSubdirs(string dir, vector<string> &subdirs) {
    // This function writes all of the subdirectory names present in a lightcone
    // output directory to the string vector 'subdirs'. The assumptions are that 
    // each subdirectory name somewhere contains either the character couple "lc" 
    // or the string "STEP", and that no non-directory items lie under dir/. Any
    // subdirs not conforming to this naming scheme are ignored.
    //
    // Params:
    // :param dir: path to a lightcone output directory
    // :param subdirs: a vector to contain the subdirectory names under dir
    // :return: none
    
    // open dir
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL) {
        cout << "Error(" << errno << ") opening lightcone data subdirs at " << dir << endl;
        return errno;
    }

    // find all items within dir/
    while ((dirp = readdir(dp)) != NULL) {
        if (string(dirp->d_name).find("lc")!=string::npos or 
            string(dirp->d_name).find("STEP")!=string::npos){ 
            subdirs.push_back(string(dirp->d_name));
        }   
    }
    closedir(dp);
    return 0;
}


//======================================================================================


int getLCFile(string dir, string &file) {
    // This functions returns the GIO header file present in a lightcone output 
    // step subdirectory (header files are those that are unhashed (#)). 
    // This function enforces that only one file header is found, implying
    // that the output of only one single lightcone step is contained in the
    // directory dir/. In short, a step-wise directory structure is assumed, 
    // as described in the documentation comments under the function header 
    // for getLCSteps(). 
    // Assumptions are that the character couple "lc" appear somewhere in the 
    // file name, and that there are no subdirectories or otherwise unhashed
    // file names present in directory dir/.
    // As a exception, .SubInput GIO files are allowed to be present, and if so, 
    // are simply ignored.
    //
    // Params:
    // :param dir: the path to the directory containing the output gio files
    // :param file: string object at which to store the found header file
    // :return: the header file found in directory dir/ as a string

    // open dir/
    DIR *dp;
    struct dirent *dirp;

    if((dp  = opendir(dir.c_str())) == NULL) {
        cout << "Error(" << errno << ") opening lightcone data files at " << dir << endl;
        return errno;
    }

    // find all header files in dir/
    vector<string> files;
    while ((dirp = readdir(dp)) != NULL) {
        if (string(dirp->d_name).find("lc") != string::npos & 
            string(dirp->d_name).find("#") == string::npos  & 
            string(dirp->d_name).find("SubInput") == string::npos   ){ 
            files.push_back(string(dirp->d_name));
        }   
    }

    // enforce exactly one header file found
    if(files.size() == 0){
        cout << "\nNo valid header files found in dir" << dir << endl;
        MPI_Abort(MPI_COMM_WORLD, 0);
    }
    if(files.size() > 1){     
        cout << "Too many header files in directory " << dir << 
                ". LC Output files should be separated by step-respective subdirectories" << endl;
        MPI_Abort(MPI_COMM_WORLD, 0);
    }

    // done
    closedir(dp);
    file = files[0];
    return 0;
}


//======================================================================================


int getLCSteps(int maxStep, int minStep, string dir, vector<string> &step_strings){
    // This function finds all simulation steps that are present in a lightcone
    // output directory that are within some specified timestep range. The 
    // expected directory structure is that given in Figure 8 of Creating 
    // Lightcones in HACC; there should be a top directory (string dir) under 
    // which is a subdirectory for each step that ran through the lightcone code. 
    // The name of these subdirectories is expected to take the form:
    //
    // {N non-digit characters, somewhere containing "lc"}{N digits composing the step number}
    // 
    // For example, lc487, and lcGals487 are valid. lightcone487_output, is not.
    // The assumptions stated in the documentation comments under the 
    // getLCSubdirs() function header are of course made here as well.
    // 
    // Params:
    // :param maxStep: the maximum step of interest (corresponding to the minimum
    //                 redshift desired to appear in the cutout)
    // :param minStep: the minimum step of interest (corresponding to the maximum
    //                 redshift desired to appear in the cutout)
    // :param dir: the path to a lightcone output directory. It is assumed that 
    //             the output data for each lightcone step are organized into 
    //             subdirectories. The expected directory structure is described 
    //             in the documentation comments under the function header for
    //             getLCSubdirs()
    // :param step_strings: a vector to contain steps found in the lightcone
    //                      output, as strings, which satisfy maxStep >= step >= minStep. 
    // :return: none

    // find all lc step subdirs
    vector<string> subdirs;
    getLCSubdirs(dir, subdirs);
    
    // extract step numbers from each subdirs
    vector<int> stepsAvail;
    for(int i=0; i<subdirs.size(); ++i){
        for(string::size_type j = 0; j < subdirs[i].size(); ++j){
            if( isdigit(subdirs[i][j]) ){
                stepsAvail.push_back( atoi(subdirs[i].substr(j).c_str()) );
                break;
            }
        }   
    }

    // fill step_strings vector (see remarks at function header)
    sort(stepsAvail.begin(), stepsAvail.end());
     
    for(int k=0; k<stepsAvail.size(); ++k){
        if(stepsAvail[k] <= maxStep and stepsAvail[k] >= minStep){
            ostringstream stepStringStream;
            stepStringStream << stepsAvail[k];
            const string& stepString = stepStringStream.str(); 
            step_strings.push_back( stepString );
        }
    }

    return 0;
}


//======================================================================================


int prepStepSubdir(string subdir_name, bool overwrite, bool print, bool verbose){
    // This function opens/creates a directory given by subdir_name and handles its
    // contents depending on the values of the following arguments. This is used
    // in processLC.cpp (use case 2) to prepare step-wise cutout subdirectories 
    // for writing. This function either
    //
    // 1) creates a non-existent directory
    // 2) opens an existing empty directory
    // 3) throws an error if the directory already exists and is non-empty
    // 4) opens and clears the contents of the directory if already exists
    //
    // Params:
    // :param subdir_name: the name of the directory to open
    // :param overwrite: if subdir_name is non-empty, rm all of its contents. If any
    //                   files ending in any extension other than .bin are found, an
    //                   error is thrown for safety
    // :param print: whether or not to print function progress to stdout
    // :param verbose: whether or not to print extra verbose output
    // :return: an error code. 0 is no error. 1 means overwrite=true, but non-binary 
    //          files were unexpectedly found. 2 means the directory was non-empty 
    //          and overwrite==false. 
    
    

    DIR *dir = opendir(subdir_name.c_str());
    struct dirent *d;
    char full_path[256];
    char *ext;
    int nf = 0;
    vector<string> files_to_remove;
    int error = 0;
    
    // doesn't exist; create the subdir
    if(dir == NULL){
        mkdir(subdir_name.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IXOTH);
        if(print){
            cout << "Created subdir: " << subdir_name << endl;
        }
    }
    
    // does exist; handle
    else{
        while((d = readdir(dir)) != NULL){ if(++nf>2){ break;} }
        
        // empty but exists; do nothing
        if(nf<=2){
            closedir(dir);
            if(print){
                cout << "Entered subdir: " << subdir_name << endl;
            }
        }

        // not empty; abort
        else if(nf > 2 and overwrite == false){
            closedir(dir);
            if(print){
                cout << "\nDirectory " << subdir_name << " is non-empty; skipping halo" << endl;
            }
            error = 2;
        }

        // not empty; overwrite
        else if(nf > 2 and overwrite == true){
            if(print){
                cout << "\nDirectory " << subdir_name << " is non-empty";
                cout << " and --overwrite flag was passed; removing binary files here" << endl;
            }
            while( (d = readdir(dir)) != NULL){
                
                // just to be safe
                if (0==strcmp(d->d_name, ".") || 
                    0==strcmp(d->d_name, "..")) { continue; }

                sprintf(full_path, "%s%s%s", subdir_name.c_str(), "/", d->d_name);
                ext = strrchr(d->d_name, '.');

                // make sure directory only contains binary files
                if( strcmp(ext, ".bin") != 0){
                    if(print){
                        cout << "\nDirectory contains non-cutout files! Not deleting anything. " << 
                                "Maybe passed wrong path?" << endl;
                    }
                    error = 1;
                    break;
                }
                files_to_remove.push_back(full_path);
            }
            
            if(error == 0){
                for(int l = 0; l < files_to_remove.size(); ++l){
                    if(verbose and print){
                        cout << "\nDELETING " << files_to_remove[l] << endl;
                    }
                    remove(files_to_remove[l].c_str());
                }
            }
        }
    }    
    closedir(dir);
    return error;
}


//======================================================================================


//////////////////////////////////////////////////////
//
//                cosmo functions
//
//////////////////////////////////////////////////////


float aToZ(float a) {
    // Converts scale factor to redshift.
    //
    // Params:
    // :param a: the scale factor
    // :return: the redshift corresponding to input a
    
    return (1.0f/a)-1.0f;
}


//======================================================================================


int zToStep(float z, int totSteps, float maxZ){
    // Function to convert a redshift to a step number, rounding 
    // toward a = 0.
    //
    // Params:
    // :param z: the input redshift
    // :totSteps: the total number of steps of the simulation of 
    //            interest. Note-- the initial conditions are not 
    //            a step! totSteps should be the maximum snapshot 
    //            number.
    // :maxZ: the initial redshift of the simulation
    // :return: the simulation step corresponding to the input redshift, 
    //          rounded toward a = 0

    float amin = 1/(maxZ + 1);
    float amax = 1.0;
    float adiff = (amax-amin)/(totSteps);
    
    float a = 1/(1+z);
    int step = floor((a-amin) / adiff);
    return step;
}


//======================================================================================


float stepToZ(float step, int totSteps, float maxZ){
    // Function to convert a step number to a redshift
    //
    // Params:
    // :param step: the input step
    // :totSteps: the total number of steps of the simulation of 
    //            interest. Note-- the initial conditions are not 
    //            a step! totSteps should be the maximum snapshot 
    //            number.
    // :maxZ: the initial redshift of the simulation
    // :return: the redshift corresponding to the input simulation step

    float amin = 1/(maxZ + 1);
    float amax = 1.0;
    float adiff = (amax-amin)/(totSteps);
    
    float a = amin + (adiff*step);
    float z = 1/a - 1;
    return z;
}


//======================================================================================


//////////////////////////////////////////////////////
//
//           matrix/vector operations
//
//////////////////////////////////////////////////////


void sizeMismatch(){ 
    cout << "\ninput vectors must have the same length" << endl;
    MPI_Abort(MPI_COMM_WORLD, 0);
}


//======================================================================================


vector<vector<float> > scalarMultiply(const vector<vector<float> > &matrix, float scalar){
    
    // Multiply a matrix by a scalar value
    //
    // Params:
    // :param matrix: a vector<vector<float>> object 
    // :param scalar: scalar float
    // :return: a vector<vector<float>> object, which is matrix*scalar
    
    int rows = matrix.size();
    int cols = matrix[0].size();
    vector<vector<float> > ans(rows, vector<float>(cols));

    for( int n = 0; n < rows; ++n )
        for( int m = 0; m < cols; ++m )
            ans[n][m] = matrix[n][m] * scalar;
    return ans;
}


//======================================================================================


vector<vector<float> > squareMat(const vector<vector<float> > &matrix){
    
    // Square a matrix
    //
    // Params:
    // :param matrix: a vector<vector<float>> object 
    // :return: a vector<vector<float>> object,which is matrix^2
    
    int rows = matrix.size();
    int cols = matrix[0].size();
    vector<vector<float> > ans(rows, vector<float>(cols));
    
    for( int n = 0; n < rows; ++n )
        for( int m = 0; m < cols; ++m )
            for( int y=0; y < rows; ++y)
                ans[n][m] += matrix[n][y] * matrix[y][m];
    return ans;
}


//======================================================================================


vector<float> matVecMul(const vector<vector<float> > &matrix, const vector<float> &vec){
    
    // Multiply a matrix by a vector
    //
    // Params:
    // :param matrix: a vector<vector<float>> object
    // :param vec: a vector<float> object
    // :return: a vector<float> object,which is matrix*vec
    
    int vecSize = vec.size();
    int rows = matrix.size();
    int cols = matrix[0].size();
    if(vecSize != rows){
        cout << "\nmatrix and vector dimensions do not match" << endl;
        MPI_Abort(MPI_COMM_WORLD, 0);
    }
    
    vector<float> ans(vecSize);
    
    for( int n = 0; n < rows; ++n )
        for( int m = 0; m < cols; ++m )
            ans[n] += (matrix[n][m] * vec[m]);
    return ans;
}


//======================================================================================


float vecPairAngle(const vector<float> &v1,
                   const vector<float> &v2){
    // Return the angle between two vectors, in radians
    //
    // Params:
    // :param v1: some three-dimensional vector
    // :param v2: some other three-dimensional vector
    // :return: the angle between v1 and v2, in radians
   
    // find (v1·v2), |v1|, and ||
    float v1dv2 = std::inner_product(v1.begin(), v1.end(), v2.begin(), 0.0);
    float mag_v1 = sqrt(std::inner_product(v1.begin(), v1.end(), v1.begin(), 0.0));
    float mag_v2 = sqrt(std::inner_product(v2.begin(), v2.end(), v2.begin(), 0.0));

    float theta = acos( v1dv2 / (mag_v1 * mag_v2) );
    return theta; 
} 


//======================================================================================


float dot(const vector<float> &v1, 
          const vector<float> &v2){
    // This function calculates the dot product of two N-dimensional vectors
    //
    // Params:
    // :param v1: some N-dimensional vector
    // :param v2: some other N-dimensional vector
    // :return: the dot product as a float
    
    int n1 = v1.size();
    int n2 = v2.size();
    if(n1 != n2){ sizeMismatch(); }
    
    float dot = 0;
    for(int i = 0; i < n1; ++i){
        dot += v1[i] * v2[i];
    }
    return dot;
}


//======================================================================================


void cross(const vector<float> &v1, 
           const vector<float> &v2,
           vector<float> &v1xv2){
    // This function calculates the cross product of two three 
    // dimensional vectors
    //
    // Params:
    // :param v1: some three-dimensional vector
    // :param v2: some other three-dimensional vector
    // :param v1xv2: vector of size 0to hold the resultant 
    //               cross-product of vectors v1 and v2
    // :return: none
    
    int n1 = v1.size();
    int n2 = v2.size();
    if(n1 != n2){ sizeMismatch(); }

    v1xv2.push_back( v1[1]*v2[2] - v1[2]*v2[1] );
    v1xv2.push_back( -(v1[0]*v2[2] - v1[2]*v2[0]) );
    v1xv2.push_back( v1[0]*v2[1] - v1[1]*v2[0] );
}


//======================================================================================


void normCross(const vector<float> &a,
           const vector<float> &b,
           vector<float> &k){
    // This function returns the normalized cross product of two three-dimensional
    // vectors. The notation, here, is chosen to match that of the Rodrigues rotation 
    // formula for the rotation vector k, rather than matching the notation of cross() 
    // above. Feel free to contact me with urgency if you find this issue troublesome.
    //
    // Parms:
    // :param a: some three-dimensional vector
    // :param b: some other three-dimensional vector
    // :param k: vector of size 0 to hold the resultant normalized cross-product of 
    //           vectors a and b
    // :return: none

    int na = a.size();
    int nb = b.size();
    if(na != nb){ sizeMismatch(); }

    vector<float> axb;
    cross(a, b, axb);
    float mag_axb = sqrt(std::inner_product(axb.begin(), axb.end(), axb.begin(), 0.0));

    for(int i=0; i<na; ++i){
        if(mag_axb == 0){ 
            k.push_back(0);
        } else {
            k.push_back( axb[i] / mag_axb );
        }
    }
}


//======================================================================================


//////////////////////////////////////////////////////////////////////
// 
// The following functions come from vvector.h in the OpenGL Utility 
// Toolkit (glut). I've rewritten them here as functions rather than 
// macros for ease of use. See 
// https://github.com/markkilgard/glut/blob/master/lib/gle/vvector.h
//
//  vvector.h
// 
//  FUNCTION:
//  This file contains a number of utilities useful for handling
//  3D vectors
//  
//  HISTORY:
//  Written by Linas Vepstas, August 1991
//  Added 2D code, March 1993
//  Added Outer products, C++ proofed, Linas Vepstas October 1993
//
// ///////////////////////////////////////////////////////////////////


double determinant_3x3(const vector<vector<float> > &m){
   
   // Computes the determinant of a 3x3 matrix
   //
   // Params:
   // :param m: a vector<vector<float>> object (matrix)
   // :return: the determinant of m as a double

    double d; 
    d = m[0][0] * (m[1][1]*m[2][2] - m[1][2] * m[2][1]);		
    d -= m[0][1] * (m[1][0]*m[2][2] - m[1][2] * m[2][0]);	
    d += m[0][2] * (m[1][0]*m[2][1] - m[1][1] * m[2][0]);	
    return d;
}


//======================================================================================


vector<vector<float> > scale_adjoint_3x3(const vector<vector<float> > &m, float s){
    
    // Computes the adjoint of a 3x3 matrix, and scales it
    //
    // Params:
    // :param m: a vector<vector<float>> object (matrix)
    // :param s: the scaling factor

    vector<vector<float> > a(3, vector<float>(3));

    a[0][0] = (s) * (m[1][1] * m[2][2] - m[1][2] * m[2][1]);	
    a[1][0] = (s) * (m[1][2] * m[2][0] - m[1][0] * m[2][2]);	
    a[2][0] = (s) * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);	
								
    a[0][1] = (s) * (m[0][2] * m[2][1] - m[0][1] * m[2][2]);	
    a[1][1] = (s) * (m[0][0] * m[2][2] - m[0][2] * m[2][0]);	
    a[2][1] = (s) * (m[0][1] * m[2][0] - m[0][0] * m[2][1]);	
								
    a[0][2] = (s) * (m[0][1] * m[1][2] - m[0][2] * m[1][1]);	
    a[1][2] = (s) * (m[0][2] * m[1][0] - m[0][0] * m[1][2]);	
    a[2][2] = (s) * (m[0][0] * m[1][1] - m[0][1] * m[1][0]);
    return a;
}


//======================================================================================


vector<vector<float> > invert_3x3(const vector<vector<float> > &m){

    // Inverts a 3x3 matrix
    //
    // Params:
    // :param m: the matrix to invert as a vector<vector<float>> object
    // :return: the 3x3 matrix inversion of m

    double det_inv;					    
    double  det = determinant_3x3(m);	    
    det_inv = 1.0 / (det);				
    vector<vector<float> > m_inv;
    m_inv = scale_adjoint_3x3(m, det_inv);
    return m_inv;
}


//======================================================================================


//////////////////////////////////////////////////////
//
//         coord rotation functions
//
//////////////////////////////////////////////////////


void cross_prod_matrix(const vector<float> &k, vector<vector<float> > &K){

    // Computes the cross-product matrix, K, for the unit vector k, 
    // as used in the Rodrigues rotation formula
    //
    // Parms:
    // :param k: the axis of rotation
    // :param K: a vector<vector(float>> object in which to store the result
    // :return: None

    float Karr[3][3] = { 
        { 0.0, -k[2], k[1]},
        { k[2], 0.0, -k[0]},
        {-k[1], k[0], 0.0 }
    };
    
    K.push_back( vector<float>(Karr[0], Karr[0]+3) );
    K.push_back( vector<float>(Karr[1], Karr[1]+3) );
    K.push_back( vector<float>(Karr[2], Karr[2]+3) );
}


//======================================================================================


void rotation_matrix(const vector<vector<float> > &K, const float B, vector<vector<float> > &R){
    
    // Computes the rotation matrix, R, as found in the Rodrigues 
    // rotation formula
    //
    // Parms:
    // :param K: the cross-product matrix for the unit vector of rotation k
    // :param B: the angle of rotation
    // :param R: a vector<vector<float>> object in which to store the result
    // :return: None
    
    vector<vector<float> > K2 = squareMat(K);
    vector<vector<float> > Ksin = scalarMultiply(K, sin(B));
    vector<vector<float> > K2cos = scalarMultiply(K2, 1-cos(B));
    
    // sum components above to find rotation matrix
    float Rarr[3][3] = {
        { 1.0+Ksin[0][0]+K2cos[0][0],     Ksin[0][1]+K2cos[0][1],     Ksin[0][2]+K2cos[0][2]},
        {     Ksin[1][0]+K2cos[1][0], 1.0+Ksin[1][1]+K2cos[1][1],     Ksin[1][2]+K2cos[1][2]},
        {     Ksin[2][0]+K2cos[2][0],     Ksin[2][1]+K2cos[2][1], 1.0+Ksin[2][2]+K2cos[2][2]}
    };
    
    R.push_back( vector<float>(Rarr[0], Rarr[0]+3) );
    R.push_back( vector<float>(Rarr[1], Rarr[1]+3) );
    R.push_back( vector<float>(Rarr[2], Rarr[2]+3) );  
}
