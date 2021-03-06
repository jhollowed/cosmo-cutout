// Joe Hollowed
// CPAC 2018

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
#include <sys/time.h>

// source files
#include "util.h"
#include "processLC.h"

using namespace std;

//////////////////////////////////////////////////////
//
//            driver function
//
//////////////////////////////////////////////////////

int main( int argc, char** argv ) {

    // This code generates a cutout from a larger lightcone run by finding all
    // particles/objects residing within a volume defined by ɵ and ϕ bounds 
    // in spherical coordaintes, centered on the observer (assumed to be the 
    // spatial origin)
    //
    // Two use cases are supported: 
    //
    /////////////////////////////////////////
    //
    // Define the ɵ and ϕ bounds explicitly:
    // 
    // lc_cutout <input lightcone dir> <output dir> <min redshift> <max redshift> 
    // --theta <ɵ_center> <dɵ> --phi <ϕ_center> <dϕ>
    // 
    // where the ϕ_center argument is the azimuthal coordinate of the center of the 
    // field of view that one wishes to cut out of the lightcone, and dϕ is the 
    // angualar distance from this center to the edge of the cutout and likewise 
    // for the similar ɵ args. That is, the result will be a sky area that spans 
    // 
    // (ɵ_center - dɵ)  ≤  ɵ    ≤  (ɵ_center + dɵ)
    // (ϕ_center - dϕ)  ≤  ɵ    ≤  (ϕ_center + dϕ)
    //
    // The expected units are DEGREES. The --theta and --phi flags can be replaced 
    // with -t and -p
    //
    /////////////////////////////////////////
    //
    // Allow the ɵ and ϕ bounds to be computed interanally to obtain a cutout of 
    // a certain width, in Mpc/h, centered on a certain cartesian positon, 
    // (x, y, z) Mpc/h (intended to be used for cutting out cluster-sized halos):
    //
    // lc_cutout <input lightcone dir> <output dir> <min redshift> <max redshift> 
    // --halo <x> <y> <z> --boxLength <box length in arcmin>
    //
    // The --halo and --boxLength flags can be replaced with -h and -b
    //
    // We want to express the positions of  all of our lightcone objects in 
    // spherical coordinates, to perform the cutout, and we want that coordinate 
    // system to be rotated such that the halo of intererst lies on the equator at
    // 
    // (r_halo, 90°, 0°)
    // where r_halo = (x^2 + y^2 + z^2)^(1/2)
    // 
    // Let's call the position vector of the halo before this rotation 
    // a = [x, y, z], and after, b = [x_rot, y_rot, z_rot] = [r_halo, 0, 0]
    //
    // We perform this rotation for each lightcone object via the Rodrigues rotation
    // formula, which answers the following question: given a position vector v, a 
    // normalized axis of rotation k, and an angle of rotation β, what is an 
    // analytical form for a new vector v_rot which is v rotated by an anlge β 
    // about k?
    //
    // First, we find k by taking the cross product of two vectors defining the 
    // plane of rotation. The obvious choice of these two vectors are a and b, as 
    // defined above;
    //
    // k = (a ⨯ b) / ||a ⨯ b||
    //
    // then, for any other position vector v, v_rot is given by
    //
    // v_rot = vcosβ + (k ⨯ v)sinβ + k(k·v)(1-cosβ)
    //
    // This coordinate rotation is required because the bounding functions which 
    // define the field of view of the observer, while constant in theta-phi space, 
    // are nonlinear in cartesian space. The distortion is maximized near the poles 
    // of the spherical coordinate system, and minmized at the equator. Areas 
    // defined by constant theta-phi bounds then appear trapezoidal to the observer 
    // when far from the equator. It is important that our cutout areas are 
    // maintained as square for at least two reasons:
    //
    // - At the moment, fft restrictions of flat-sky lensing code require that the 
    //   cutout is square
    // - The cutouts returned will not actually have all side lengths of boxLength 
    //   if we don't do this rotation, which the user explicitly requested
    //
    ////////////////////////////////////////
    //
    // Note that the first use case describe does not perform the coordinate 
    // rotation which is described in the second. So, cutouts returned will
    // not necessarily be square, or symmetrical.
    //
    // Additional flags that can be passed are:
    //
    // --verbose: lots more output, useful for debugging
    // --timeit: reading, redistribution, computation, and write out will all be
    //           timed with MPI wall time, and reported in output
    // --overwrite: if write-out directories are not empty (there are binary files
    //              present there from a previous cutout run), then overwrite the
    //              contents rather than stopping execution with an error. If anything
    //              other than binary files are found in the directory, then an error
    //              is still raised and nothing is deleted.
    // --posOnly: only output "first-order" particle quantities to the resultant cutout, 
    //            including x, y, z, a, and id. vx, vy, vz, replication, and rotation 
    //            will be ommitted. Should speed up redistribution step.
    // 
    // All four of these additional options default to false

    // start MPI 
    MPI_Init(&argc, &argv);
    int myrank, numranks;
    MPI_Comm_size(MPI_COMM_WORLD, &numranks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if(myrank == 0){ cout << "\n\n---------- Starting on " << numranks << 
                           " MPI ranks ----------" << endl; }
    char cart[3] = {'x', 'y', 'z'};

    string input_lc_dir, out_dir;
    input_lc_dir = string(argv[1]);
    out_dir = string(argv[2]);

    // check format of in/out dirs
    if(input_lc_dir[input_lc_dir.size()-1] != '/'){
        ostringstream modified_in_dir;
        modified_in_dir << input_lc_dir << "/";
        input_lc_dir = modified_in_dir.str();
    }
    if(out_dir[out_dir.size()-1] != '/'){
        ostringstream modified_out_dir;
        modified_out_dir << out_dir << "/";
        out_dir = modified_out_dir.str();
    }
    struct stat info;
    if( stat( out_dir.c_str(), &info ) != 0 ){
        cout << "\nCannot access specified output directory " << out_dir;
        MPI_Abort(MPI_COMM_WORLD, 0);
    }

    if(myrank == 0){ cout << "Using lightcone at ";  
                   cout << input_lc_dir << endl;
                   cout << "Writing at ";
                   cout << out_dir << endl;}

    // build step_strings vector by locating the step present in the lightcone
    // data directory that is nearest the redshift (or snapshot) requested by the user
    // If args 2 and 3 are both above 20, assume the bounds are being given by snapshot 
    // number, otherwise assume redshift
    float minBound = atof(argv[3]);
    float maxBound = atof(argv[4]);
    int minStep;
    int maxStep;
    float minZ;
    float maxZ;
    
    if(minBound > 20 and maxBound > 20){
        maxStep = int(minBound);
        minStep = int(maxBound);
        maxZ = stepToZ(minStep);
        minZ = stepToZ(maxStep);
    }
    else{
        minZ = minBound;
        maxZ = maxBound;
        maxStep = zToStep(minZ);    
        minStep = zToStep(maxZ);
    }
    
    vector<string> step_strings;
    getLCSteps(maxStep, minStep, input_lc_dir, step_strings);
    if(myrank == 0){ 
        cout << "MAX STEP: " << maxStep << endl;
        cout << "MIN STEP: " << minStep << endl;
        cout << "steps to include from z= " << minZ << " to z=" << maxZ << ": " << endl;
        for(int i=0; i<step_strings.size(); ++i){ cout << step_strings[i] << " ";}
        cout << endl;
    }

    // might not use all of these but whatever
    vector<float> theta_cut(2);
    vector<float> phi_cut(2);
    vector<float> haloPos;
    vector<float> haloProps;
    vector<string> haloTags;
    float boxLength;
    bool verbose = false;
    bool timeit = false;
    bool overwrite = false;
    bool positionOnly = false;
    bool forceWriteProps = false;
    bool propsOnly = false;
    string massDef="sod";

    // check that supplied arguments are valid
    vector<string> args(argv+1, argv + argc);
    bool customThetaBounds = int((find(args.begin(), args.end(), "-t") != args.end()) ||
            (find(args.begin(), args.end(), "--theta") != args.end()));
    bool customPhiBounds = int((find(args.begin(), args.end(), "-p") != args.end()) ||
            (find(args.begin(), args.end(), "--phi") != args.end()));
    bool customHalo = int((find(args.begin(), args.end(), "-h") != args.end()) ||
            (find(args.begin(), args.end(), "--halo") != args.end()));
    bool customHaloFile = int((find(args.begin(), args.end(), "-f") != args.end()) ||
            (find(args.begin(), args.end(), "--haloFile") != args.end()));
    bool customBox = int((find(args.begin(), args.end(), "-b") != args.end()) ||
            (find(args.begin(), args.end(), "--boxLength") != args.end()));
    bool customMassDef = int((find(args.begin(), args.end(), "-m") != args.end()) ||
            (find(args.begin(), args.end(), "--massDef") != args.end()));

    // there are two general use cases of this cutout code, as described in the 
    // docstring below the declaration of this main function. Here, the program aborts
    // to prevent confused input arguments which mix those two use cases.
    if( (customHalo || customHaloFile) ^ customBox ){ 
        cout << "\n-h (or -f) and -b options must accompany eachother" << endl;
        MPI_Abort(MPI_COMM_WORLD, 0);
    }
    if( customMassDef && !customHaloFile){
        cout << "\n-m does nothing if not used along with -f";
        MPI_Abort(MPI_COMM_WORLD, 0);
    }
    if( customThetaBounds ^ customPhiBounds ){
        cout << "\n-t and -p options must accompany eachother";
        MPI_Abort(MPI_COMM_WORLD, 0);
    }
    if( (customHalo & customThetaBounds) || (customHaloFile & customThetaBounds) ){
        cout << "\n-t and -p options must not be used in the case " <<
                "that -h (or -f) and -b arguments are passed";
        MPI_Abort(MPI_COMM_WORLD, 0);
    }
    if( !customThetaBounds && !customPhiBounds && !customHalo && !customHaloFile && !customBox ){
        cout << "\nValid options are -h, -f, -b, -t, -p, -v, -m, and --timeit. See github readme for help";
        MPI_Abort(MPI_COMM_WORLD, 0);
    }

    // search argument vector for options, update default parameters if found
    // Note to future self: strcmp() returns 0 if the input strings are equal. 
    for(int i=5; i<argc; ++i){

        if(strcmp(argv[i],"-t")==0 || strcmp(argv[i],"--theta")==0){
            float theta_center = strtof(argv[++i], NULL) * ARCSEC;
            float dtheta = strtof(argv[++i], NULL) * ARCSEC;
            theta_cut[0] = theta_center - dtheta;
            theta_cut[1] = theta_center + dtheta;
        }
        else if(strcmp(argv[i],"-p")==0 || strcmp(argv[i],"--phi")==0){
            float phi_center = strtof(argv[++i], NULL) * ARCSEC;
            float dphi  = strtof(argv[++i], NULL) * ARCSEC;
            phi_cut[0] = phi_center - dphi;
            phi_cut[1] = phi_center + dphi;
        }
        else if (strcmp(argv[i],"-b")==0 || strcmp(argv[i],"--boxLength")==0){
            boxLength = strtof(argv[++i], NULL);
        }
        else if (strcmp(argv[i],"-m")==0 || strcmp(argv[i],"--massDef")==0){
            massDef = argv[++i];
        }
        else if(strcmp(argv[i],"-h")==0 || strcmp(argv[i],"--halo")==0){
            haloPos.push_back(strtof(argv[++i], NULL));
            haloPos.push_back(strtof(argv[++i], NULL));
            haloPos.push_back(strtof(argv[++i], NULL));
            // set halo properties to -1 (unknown)
            haloProps.push_back(-1.0); 
            haloProps.push_back(-1.0); 
            haloProps.push_back(-1.0); 
            if(massDef == "sod")
                haloProps.push_back(-1.0); 
                haloProps.push_back(-1.0); 
                haloProps.push_back(-1.0); 
        }
        else if(strcmp(argv[i],"-f")==0 || strcmp(argv[i],"--haloFile")==0){
            string haloFileName(argv[++i]);
            readHaloFile(haloFileName, haloPos, haloTags, haloProps, massDef);
        }
        else if (strcmp(argv[i],"-v")==0 || strcmp(argv[i],"--verbose")==0){
            verbose = true;
        }
        else if (strcmp(argv[i],"--timeit") == 0){
            timeit = true;
        }
        else if (strcmp(argv[i],"--overwrite") == 0){
            overwrite = true;
        }
        else if (strcmp(argv[i],"--posOnly") == 0){
            positionOnly = true;
        }
        else if (strcmp(argv[i],"--forceWriteProps") == 0){
            forceWriteProps = true;
        }
        else if (strcmp(argv[i],"--propsOnly") == 0){
            propsOnly = true;
        }
    }

    // if customHaloFile == true, then create an output subdirectory per halo in out_dir
    vector<string> halo_out_dirs;
    if(customHaloFile){
        for(int h=0; h<haloTags.size(); ++h){

            ostringstream halo_subdir;
            halo_subdir << out_dir << "halo_" << haloTags[h] << "/";
            
            DIR *dir = opendir(halo_subdir.str().c_str());
            struct dirent *d;
            int nf = 0; 

            // if subdir already exists, make sure it's valid
            if(dir != NULL){
                while((d = readdir(dir)) != NULL){ if(++nf>2){ break;} }
                closedir(dir);
            }
            // Otherwise, create the subdir
            else{
                mkdir(halo_subdir.str().c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IXOTH);
                if(myrank == 0 and h%100 == 0){ 
                    cout << "Created output directory: " << halo_subdir.str() << endl; 
                }
            }
            halo_out_dirs.push_back(halo_subdir.str());
        }
    }else if(customHalo){
        // the processLC function which takes halo positions as an argument
        // expects the output directorie(s) to be given as a vector
        halo_out_dirs.push_back(out_dir);
    }
    
    // print summary 
    if(myrank == 0){
        if(customHalo || customHaloFile){
            cout << endl << haloPos.size()/3 << " target halo(s) will be cut out of the lightcone";
            
            if(verbose == true){
                cout << ":" << endl;
                for(int k=0; k<haloTags.size(); ++k){
                    cout << "Halo " << haloTags[k] << ": " << endl;
                    for(int i=0;i<3;++i){ cout << cart[i] << "=" << haloPos[k+i] << " ";}
                }
            }
            cout << endl << "box length: " << boxLength << " arcmin" << endl;
        
        }else{
            cout << "theta bounds: ";
            cout << theta_cut[0]/ARCSEC << " -> " << theta_cut[1]/ARCSEC <<" deg"<< endl;
            cout << "phi bounds: ";
            cout << phi_cut[0]/ARCSEC << " -> " << phi_cut[1]/ARCSEC <<" deg" << endl;
        }


        cout << "\nverbose is set to " << verbose << endl;
        cout << "timeit is set to " << timeit << endl;
        cout << "overwrite is set to " << overwrite << endl;
        cout << "posOnly is set to " << positionOnly << endl;
    }

    // call overloaded processing function

    double start;
    double stop;
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    if(customHalo || customHaloFile){
        processLC(input_lc_dir, halo_out_dirs, step_strings, haloPos, haloProps, 
                  boxLength, myrank, numranks, verbose, timeit, overwrite, positionOnly, 
                  forceWriteProps, propsOnly);
    }else{
        processLC(input_lc_dir, out_dir, step_strings, theta_cut, phi_cut, 
                  myrank, numranks, verbose, timeit, overwrite, positionOnly);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    stop = MPI_Wtime();
    
    double duration = stop - start;
    if(myrank == 0 and timeit == true){ cout << "\nExecution time: " << duration << " s" << endl; }

    MPI_Finalize();
    return 0;
}
