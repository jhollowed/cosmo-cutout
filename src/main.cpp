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

// source files
#include "util.h"
#include "processLC.h"
 
// Generic IO
//#include "GenericIO.h"

using namespace std;
//using namespace gio;


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
    // lc_cutout <input lightcone dir> <output dir> <max redshift> --theta <ɵ_center> <dɵ> --phi <ϕ_center> <dϕ>
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
    // lc_cutout <input lightcone dir> <output dir> <max redshift> --halo <x> <y> <z> --boxLength <box length>
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
    
    MPI_Init(&argc, &argv);

    cout << "\n\n---------- Starting ----------" << endl;
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
	ostringstream badOut;
	badOut << "Cannot access specified output directory " << out_dir;
        throw runtime_error( badOut.str() );
    }
    
    cout << "using lightcone at ";
    cout << input_lc_dir << endl;
     
    // build step_strings vector by locating the step present in the lightcone
    // data directory that is nearest the redshift requested by the user
    float maxZ = atof(argv[3]);
    int minStep = zToStep(maxZ);    
    vector<string> step_strings;
    cout << "MINSTEP: " << minStep << endl;
    getLCSteps(minStep, input_lc_dir, step_strings);
    cout << "steps to include to z=" << maxZ << ": ";
    for(int i=0; i<step_strings.size(); ++i){ cout << step_strings[i] << " ";}
    cout << endl;

    // might not use all of these but whatever
    vector<float> theta_cut(2);
    vector<float> phi_cut(2);
    vector<float> haloPos(3);
    float boxLength;

    // check that supplied arguments are valid
    vector<string> args(argv+1, argv + argc);
    bool customThetaBounds = int((find(args.begin(), args.end(), "-t") != args.end()) ||
                             (find(args.begin(), args.end(), "--theta") != args.end()));
    bool customPhiBounds = int((find(args.begin(), args.end(), "-p") != args.end()) ||
                           (find(args.begin(), args.end(), "--phi") != args.end()));
    bool customHalo = int((find(args.begin(), args.end(), "-h") != args.end()) ||
                      (find(args.begin(), args.end(), "--halo") != args.end()));
    bool customBox = int((find(args.begin(), args.end(), "-b") != args.end()) ||
                     (find(args.begin(), args.end(), "--boxLength") != args.end()));
    
    // there are two general use cases of this cutout code, as described in the 
    // docstring below the declaration of this main function. Here, exceptons are 
    // thrown to prevent confused input arguments which mix those two use cases.
    if(customHalo^customBox){ 
        throw invalid_argument("-h and -b options must accompany eachother");
    }
    if(customThetaBounds^customPhiBounds){
        throw invalid_argument("-t and -p options must accompany eachother");
    }
    if(customHalo & customThetaBounds){
        throw invalid_argument("-t and -p options must not be used in the case " \
                                    "that -h and -b arguments are passed");
    }
    if(!customThetaBounds && !customPhiBounds && !customHalo && !customBox){
        throw invalid_argument("Valid options are -h, -b, -t, and -p");
    }

    // search argument vector for options, update default parameters if found
    // Note to future self: strcmp() returns 0 if the input strings are equal. 
    for(int i=4; i<argc; ++i){

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
        else if(strcmp(argv[i],"-h")==0 || strcmp(argv[i],"--halo")==0){
            haloPos[0] = strtof(argv[++i], NULL);
            haloPos[1] = strtof(argv[++i], NULL);
            haloPos[2] = strtof(argv[++i], NULL);  
        }
        else if (strcmp(argv[i],"-b")==0 || strcmp(argv[i],"--boxLength")==0){
            boxLength = strtof(argv[++i], NULL);
        }
    }

    if(customHalo){
        cout << "target halo: ";
        for(int i=0;i<3;++i){ cout << cart[i] << "=" << haloPos[i] << " ";}
        cout << endl;
        cout << "box length: " << boxLength << " Mpc";
    }else{
        cout << "theta bounds: ";
        cout << theta_cut[0]/ARCSEC << " -> " << theta_cut[1]/ARCSEC <<" deg"<<endl;
        cout << "phi bounds: ";
        cout << phi_cut[0]/ARCSEC << " -> " << phi_cut[1]/ARCSEC <<" deg";
    }
        
    // call overloaded processing function
    if(customHalo){
        processLC(input_lc_dir, out_dir, step_strings, haloPos, boxLength);
    }else{
        processLC(input_lc_dir, out_dir, step_strings, theta_cut, phi_cut);
    }

    MPI_Finalize();
    return 0;
}
