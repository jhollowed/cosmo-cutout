#include "processLC.h"

using namespace std;
using namespace gio;

//////////////////////////////////////////////////////

struct Buffers {
    
    // From LC output
    vector<float> x;
    vector<float> y;
    vector<float> z;
    vector<float> vx;
    vector<float> vy;
    vector<float> vz;
    vector<float> a;
    vector<int> step;
    vector<int64_t> id;
    vector<int> rotation;
    vector<int32_t> replication;
    // New data columns
    vector<float> theta;
    vector<float> phi;
};


//////////////////////////////////////////////////////
//
//                Cutout function
//             	    Use Case 1
//           Custom theta - phi bounds
//
//////////////////////////////////////////////////////

void processLC(string dir_name, string out_dir, vector<string> step_strings, 
               vector<float> theta_cut, vector<float> phi_cut){

    ///////////////////////////////////////////////////////////////
    //
    //                          Setup
    //
    ///////////////////////////////////////////////////////////////
    
    Buffers b;

    // find all lc sub directories for each step in step_strings
    cout << endl << endl;
    cout << "Reading directory: " << dir_name << endl;
    vector<string> subdirs;
    getLCSubdirs(dir_name, subdirs);
    cout << "Found subdirs:" << endl;
    for (vector<string>::const_iterator i = subdirs.begin(); i != subdirs.end(); ++i)
         cout << *i << ' ';

    // find the prefix (chars before the step number) in the subdirectory names.
    // It is assumed that all subdirs have the same prefix.
    string subdirPrefix;
    for(string::size_type j = 0; j < subdirs[0].size(); ++j){
        if( isdigit(subdirs[0][j]) == true){
            subdirPrefix = subdirs[0].substr(0, j);
            break;
        }
    }
    
    ///////////////////////////////////////////////////////////////
    //
    //                 Loop over step subdirs
    //
    ///////////////////////////////////////////////////////////////
    
    
    // perform cutout on data from each lc output step
    size_t max_size = 0;
    int step;
    for (int i=0; i<step_strings.size();++i){
        
        // continue if this is the step at z=0 (lightcone volume zero)
        step =atoi(step_strings[i].c_str());
        if(step == 499){ continue;}
        
        // find header file
        cout<< "\n\n---------- Working on step " << step_strings[i] << "----------" << endl;
        string file_name;
        ostringstream file_name_stream;
        file_name_stream << dir_name << subdirPrefix << step_strings[i];
        
        getLCFile(file_name_stream.str(), file_name);
        file_name_stream << "/" << file_name; 
       
        // create gio reader and open file header 
        cout << "Opening file: " << file_name_stream.str() << endl;
        GenericIO reader(MPI_COMM_SELF, file_name_stream.str());
        reader.openAndReadHeader(GenericIO::MismatchRedistribute);
        
        // set size of buffers to be the size required by the largest data rank
        int nRanks = reader.readNRanks();
        size_t current_size;
        for (int j=0; j<nRanks; ++j) {
            current_size = reader.readNumElems(j);
            max_size = current_size > max_size ? current_size : max_size;
        }
        max_size +=10;
        cout<< "max size: " << max_size << endl; 
        b.x.resize(max_size);
        b.y.resize(max_size);
        b.z.resize(max_size);
        b.vx.resize(max_size);
        b.vy.resize(max_size);
        b.vz.resize(max_size);
        b.a.resize(max_size);
        b.id.resize(max_size);
        b.step.resize(max_size);
        b.rotation.resize(max_size);
        b.replication.resize(max_size);
        b.theta.resize(max_size);
        b.phi.resize(max_size);
        cout<<"done resizing"<<endl;
        
        // create cutout subdirectory for this step
        ostringstream step_subdir;
        step_subdir << out_dir << subdirPrefix << "Cutout" << step_strings[i] << "/";
        mkdir(step_subdir.str().c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IXOTH);
        cout << "Created subdir: " << step_subdir.str() << endl;

        ///////////////////////////////////////////////////////////////
        //
        //           Create output files + start reading
        //
        ///////////////////////////////////////////////////////////////
        
        // create binary files for cutout output
        ofstream id_file;
        ofstream theta_file;
        ofstream phi_file;
        ofstream a_file;
        ofstream x_file;
        ofstream y_file;
        ofstream z_file;
        ofstream vx_file;
        ofstream vy_file;
        ofstream vz_file;
        ofstream rotation_file;
        ofstream replication_file;

        ostringstream id_file_name;
        ostringstream theta_file_name;
        ostringstream phi_file_name;
        ostringstream a_file_name;
        ostringstream x_file_name;
        ostringstream y_file_name;
        ostringstream z_file_name;
        ostringstream vx_file_name;
        ostringstream vy_file_name;
        ostringstream vz_file_name;
        ostringstream rotation_file_name;
        ostringstream replication_file_name;

        id_file_name << step_subdir << "/id." << step << ".bin";
        theta_file_name << step_subdir << "/theta." << step << ".bin";
        phi_file_name << step_subdir << "/phi." << step << ".bin";
        a_file_name << step_subdir << "/a." << step << ".bin";
        x_file_name << step_subdir << "/x."<< step <<".bin";
        y_file_name << step_subdir << "/y."<< step <<".bin";
        z_file_name << step_subdir << "/z."<< step <<".bin";
        vx_file_name << step_subdir << "/vx."<< step <<".bin";
        vy_file_name << step_subdir << "/vy."<< step <<".bin";
        vz_file_name << step_subdir << "/vz."<< step <<".bin";
        vz_file_name << step_subdir << "/rotation."<< step <<".bin";
        vz_file_name << step_subdir << "/replication."<< step <<".bin";
        
        cout<<"starting to open files"<<endl;
        id_file.open(id_file_name.str().c_str(), ios::out | ios::binary);
        theta_file.open(theta_file_name.str().c_str(), ios::out | ios::binary);
        phi_file.open(phi_file_name.str().c_str(), ios::out | ios::binary);
        a_file.open(a_file_name.str().c_str(), ios::out | ios::binary);
        x_file.open(x_file_name.str().c_str(), ios::out | ios::binary);
        y_file.open(y_file_name.str().c_str(), ios::out | ios::binary);
        z_file.open(z_file_name.str().c_str(), ios::out | ios::binary);
        vx_file.open(vx_file_name.str().c_str(), ios::out | ios::binary);
        vy_file.open(vy_file_name.str().c_str(), ios::out | ios::binary);
        vz_file.open(vz_file_name.str().c_str(), ios::out | ios::binary);
        rotation_file.open(rotation_file_name.str().c_str(), ios::out | ios::binary);
        replication_file.open(replication_file_name.str().c_str(), ios::out | ios::binary);
        cout<<"done opening files"<<endl;
      
        // start reading 
        reader.addVariable("x", &b.x[0]); 
        reader.addVariable("y", &b.y[0]); 
        reader.addVariable("z", &b.z[0]); 
        reader.addVariable("vx", &b.vx[0]); 
        reader.addVariable("vy", &b.vy[0]); 
        reader.addVariable("vz", &b.vz[0]); 
        reader.addVariable("a", &b.a[0]); 
        reader.addVariable("step", &b.step[0]); 
        reader.addVariable("id", &b.id[0]); 
        reader.addVariable("rotation", &b.rotation[0]); 
        reader.addVariable("replication", &b.replication[0]); 

        ///////////////////////////////////////////////////////////////
        //
        //                 Do cutting and write
        //
        ///////////////////////////////////////////////////////////////

        for (int j=0; j<nRanks; ++j) {
            
            size_t current_size = reader.readNumElems(j);
            if(j%20==0){ cout << "Reading: " << current_size << endl; }
            reader.readData(j, false);
    
            if(j%20==0){ cout << "Converting positions..." << j+1 << "/" << nRanks << endl; }
            for (int k=0; k<current_size; ++k) {
                
                // limit cutout to the first octant
                if (b.x[k] > 0.0 && b.y[k] > 0.0 && b.z[k] > 0.0) {
                    
                    // spherical coordinate transformation
                    float r = (float)sqrt(b.x[k]*b.x[k] + b.y[k]*b.y[k] + b.z[k]*b.z[k]);
                    b.theta[k] = acos(b.z[k]/r) * 180.0 / PI * ARCSEC;
                    b.phi[k] = atan(b.y[k]/b.x[k]) * 180.0 / PI * ARCSEC;
                    
                    // do cut and write
                    if (b.theta[k] > theta_cut[0] && b.theta[k] < theta_cut[1] && 
                        b.phi[k] > phi_cut[0] && b.phi[k] < phi_cut[1] ) {
                        
                        id_file.write( (char*)&b.id[k], sizeof(int64_t));
                        theta_file.write( (char*)&b.theta[k], sizeof(float));
                        phi_file.write( (char*)&b.phi[k], sizeof(float));
                        x_file.write((char*)&b.x[k],sizeof(float));
                        y_file.write((char*)&b.y[k],sizeof(float));
                        z_file.write((char*)&b.z[k],sizeof(float));
                        vx_file.write((char*)&b.vx[k],sizeof(float));
                        vy_file.write((char*)&b.vy[k],sizeof(float));
                        vz_file.write((char*)&b.vz[k],sizeof(float));
                        a_file.write( (char*)&b.a[k], sizeof(float));
                        rotation_file.write( (char*)&b.rotation[k], sizeof(int));
                        replication_file.write( (char*)&b.replication[k], sizeof(int32_t));
                    }
                }
            }
        }
        
        reader.close();
        id_file.close();
        theta_file.close();
        phi_file.close();
        x_file.close();
        y_file.close();
        z_file.close();
        vx_file.close();
        vy_file.close();
        vz_file.close();
        a_file.close();
        rotation_file.close();
        replication_file.close();
    }
}

//////////////////////////////////////////////////////
//
//                Cutout function
//             	    Use Case 2
//               Custom halo cutout
//
//////////////////////////////////////////////////////

void processLC(string dir_name, string out_dir, vector<string> step_strings, 
               vector<float> halo_pos, float boxLength){

    ///////////////////////////////////////////////////////////////
    //
    //                          Setup
    //
    ///////////////////////////////////////////////////////////////
    
    Buffers b;

    // find all lc sub directories for each step in step_strings
    cout << endl << endl;
    cout << "Reading directory: " << dir_name << endl;
    vector<string> subdirs;
    getLCSubdirs(dir_name, subdirs);
    cout << "Found subdirs:" << endl;
    for (vector<string>::const_iterator i = subdirs.begin(); i != subdirs.end(); ++i)
         cout << *i << ' ';

    // find the prefix (chars before the step number) in the subdirectory names.
    // It is assumed that all subdirs have the same prefix.
    string subdirPrefix;
    for(string::size_type j = 0; j < subdirs[0].size(); ++j){
        if( isdigit(subdirs[0][j]) == true){
            subdirPrefix = subdirs[0].substr(0, j);
            break;
        }
    }
    
    ///////////////////////////////////////////////////////////////
    //
    //                  Start coordinate rotation
    //
    ///////////////////////////////////////////////////////////////

    cout<< "\n\n---------- Setting up for coordinate rotation ----------" << endl;
    // do coordinate rotation to center halo at (r, 90, 0) in spherical coords
    float halo_r = (float)sqrt(halo_pos[0]*halo_pos[0] + 
                               halo_pos[1]*halo_pos[1] + 
                               halo_pos[2]*halo_pos[2]);
    float tmp[] = {halo_r, 0, 0};
    vector<float> rotated_pos(tmp, tmp+3);
    cout << "Finding axis of rotation to move (" << 
            halo_pos[0]<< ", " << halo_pos[1]<< ", " << halo_pos[2]<< ") to (" <<
            rotated_pos[0] << ", " << rotated_pos[1] << ", " << rotated_pos[2]<<")" << endl;

    // get angle and axis of rotation -- this only needs to be calculated once for all
    // steps, and it will be used to rotate all other position vectors in the 
    // loops below
    vector<float> k;
    normCross(halo_pos, rotated_pos, k);
    float B = vecPairAngle(halo_pos, rotated_pos);
    cout << "Rotation is " << B*(180/PI) << "° about axis k = (" << 
            k[0]<< ", " << k[1]<< ", " << k[2]<< ")" << endl;

    // calculate theta_cut and phi_cut, in arcsec, given the specified boxLength
    float halfBoxLength = boxLength / 2.0;
    float dtheta = atan(halfBoxLength / halo_r);
    float dphi = dtheta;
     
    // calculate theta-phi bounds of cutout under coordinate rotation
    vector<float> theta_cut(2);
    theta_cut[0] = (PI/2 - dtheta) * 180.0/PI * ARCSEC;
    theta_cut[1] = (PI/2 + dtheta) * 180.0/PI * ARCSEC;
    vector<float> phi_cut(2);
    phi_cut[0] = (0 - dphi) * 180.0/PI * ARCSEC;
    phi_cut[1] = (0 + dphi) * 180.0/PI * ARCSEC;
    cout << "theta bounds set to: ";
    cout << theta_cut[0]/ARCSEC << "° -> " << theta_cut[1]/ARCSEC <<"°"<< endl;
    cout << "phi bounds set to: ";
    cout << phi_cut[0]/ARCSEC << "° -> " << phi_cut[1]/ARCSEC <<"°" << endl;
    cout << "theta-phi bounds result in box width of " << 
            tan(dtheta) * halo_r * 2 << " Mpc at distance to halo of " << halo_r << endl 
            << "        " << "= " << dtheta*2*180.0/PI << "°x" << dphi*2*180.0/PI << "° field of view";
    
    ///////////////////////////////////////////////////////////////
    //
    //                 Loop over step subdirs
    //
    ///////////////////////////////////////////////////////////////
    
    // perform cutout on data from each lc output step
    size_t max_size = 0;
    int step;
    for (int i=0; i<step_strings.size();++i){
        
        // continue if this is the step at z=0 (lightcone volume zero)
        step =atoi(step_strings[i].c_str());
        if(step == 499){ continue;}
        
        // find header file
        cout<< "\n\n---------- Working on step " << step_strings[i] << "----------" << endl;
        string file_name;
        ostringstream file_name_stream;
        file_name_stream << dir_name << subdirPrefix << step_strings[i];
        
        getLCFile(file_name_stream.str(), file_name);
        file_name_stream << "/" << file_name; 
       
        // create gio reader and open file header 
        cout << "Opening file: " << file_name_stream.str() << endl;
        GenericIO reader(MPI_COMM_SELF, file_name_stream.str());
        reader.openAndReadHeader(GenericIO::MismatchRedistribute);
        
        // set size of buffers to be the size required by the largest data rank
        int nRanks = reader.readNRanks();
        size_t current_size;
        for (int j=0; j<nRanks; ++j) {
            current_size = reader.readNumElems(j);
            max_size = current_size > max_size ? current_size : max_size;
        }
        max_size +=10;
        cout<< "max size: " << max_size << endl; 
        b.x.resize(max_size);
        b.y.resize(max_size);
        b.z.resize(max_size);
        b.vx.resize(max_size);
        b.vy.resize(max_size);
        b.vz.resize(max_size);
        b.a.resize(max_size);
        b.id.resize(max_size);
        b.step.resize(max_size);
        b.rotation.resize(max_size);
        b.replication.resize(max_size);
        b.theta.resize(max_size);
        b.phi.resize(max_size);
        cout<<"done resizing"<<endl;
        
        // create cutout subdirectory for this step
        ostringstream step_subdir;
        step_subdir << out_dir << subdirPrefix << "Cutout" << step_strings[i];
        mkdir(step_subdir.str().c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IXOTH);
        cout << "Created subdir: " << step_subdir.str() << endl;
        
        ///////////////////////////////////////////////////////////////
        //
        //           Create output files + start reading
        //
        ///////////////////////////////////////////////////////////////
        
        // create binary files for cutout output
        ofstream id_file;
        ofstream theta_file;
        ofstream phi_file;
        ofstream a_file;
        ofstream x_file;
        ofstream y_file;
        ofstream z_file;
        ofstream vx_file;
        ofstream vy_file;
        ofstream vz_file;
        ofstream rotation_file;
        ofstream replication_file;
        ofstream thetaRot_file;
        ofstream phiRot_file;

        ostringstream id_file_name;
        ostringstream theta_file_name;
        ostringstream phi_file_name;
        ostringstream a_file_name;
        ostringstream x_file_name;
        ostringstream y_file_name;
        ostringstream z_file_name;
        ostringstream vx_file_name;
        ostringstream vy_file_name;
        ostringstream vz_file_name;
        ostringstream rotation_file_name;
        ostringstream replication_file_name;
        ostringstream thetaRot_file_name;
        ostringstream phiRot_file_name;

        id_file_name << step_subdir.str() << "/id." << step << ".bin";
        theta_file_name << step_subdir.str() << "/theta." << step << ".bin";
        phi_file_name << step_subdir.str() << "/phi." << step << ".bin";
        a_file_name << step_subdir.str() << "/a." << step << ".bin";
        x_file_name << step_subdir.str() << "/x."<< step <<".bin";
        y_file_name << step_subdir.str() << "/y."<< step <<".bin";
        z_file_name << step_subdir.str() << "/z."<< step <<".bin";
        vx_file_name << step_subdir.str() << "/vx."<< step <<".bin";
        vy_file_name << step_subdir.str() << "/vy."<< step <<".bin";
        vz_file_name << step_subdir.str() << "/vz."<< step <<".bin";
        rotation_file_name << step_subdir.str() << "/rotation."<< step <<".bin";
        replication_file_name << step_subdir.str() << "/replication."<< step <<".bin";
        thetaRot_file_name << step_subdir.str() << "/thetaRot." << step << ".bin";
        phiRot_file_name << step_subdir.str() << "/phiRot." << step << ".bin";
        
        cout<<"starting to open files"<<endl;
        id_file.open(id_file_name.str().c_str(), ios::out | ios::binary);
        theta_file.open(theta_file_name.str().c_str(), ios::out | ios::binary);
        phi_file.open(phi_file_name.str().c_str(), ios::out | ios::binary);
        a_file.open(a_file_name.str().c_str(), ios::out | ios::binary);
        x_file.open(x_file_name.str().c_str(), ios::out | ios::binary);
        y_file.open(y_file_name.str().c_str(), ios::out | ios::binary);
        z_file.open(z_file_name.str().c_str(), ios::out | ios::binary);
        vx_file.open(vx_file_name.str().c_str(), ios::out | ios::binary);
        vy_file.open(vy_file_name.str().c_str(), ios::out | ios::binary);
        vz_file.open(vz_file_name.str().c_str(), ios::out | ios::binary);
        rotation_file.open(rotation_file_name.str().c_str(), ios::out | ios::binary);
        replication_file.open(replication_file_name.str().c_str(), ios::out | ios::binary);
        thetaRot_file.open(thetaRot_file_name.str().c_str(), ios::out | ios::binary);
        phiRot_file.open(phiRot_file_name.str().c_str(), ios::out | ios::binary);
        cout<<"done opening files"<<endl;
     
        // start reading 
        reader.addVariable("x", &b.x[0]); 
        reader.addVariable("y", &b.y[0]); 
        reader.addVariable("z", &b.z[0]); 
        reader.addVariable("vx", &b.vx[0]); 
        reader.addVariable("vy", &b.vy[0]); 
        reader.addVariable("vz", &b.vz[0]); 
        reader.addVariable("a", &b.a[0]); 
        reader.addVariable("step", &b.step[0]); 
        reader.addVariable("id", &b.id[0]); 
        reader.addVariable("rotation", &b.rotation[0]); 
        reader.addVariable("replication", &b.replication[0]); 

        ///////////////////////////////////////////////////////////////
        //
        //                 Do cutting and write
        //
        ///////////////////////////////////////////////////////////////
        
        
        for (int j=0; j<nRanks; ++j) {
            
            size_t current_size = reader.readNumElems(j);
            if(j%20==0){ cout << "Reading: " << current_size << endl; }
            reader.readData(j, false);
    
            if(j%20==0){ cout << "Converting positions..." << j+1 << "/" << nRanks << endl; }
            for (int n=0; n<current_size; ++n) {
                
                // limit cutout to first octant for speed
                if (b.x[n] > 0.0 && b.y[n] > 0.0 && b.z[n] > 0.0){
                    
                    // do coordinate rotation center halo at (r, 90, 0)
                    // B and k are the angle and axis of rotation, respectively,
                    // calculated near the beginning of this function
                    float tmp[] = {b.x[n], b.y[n], b.z[n]};
                    vector<float> v(tmp, tmp+3);
                    vector<float> v_rot;
                    rotate(k, B, v, v_rot);

                    // spherical coordinate transformation. Don't write these theta
                    // and phi values to their respective buffers, since we want to save
                    // the true theta and phi to file, not these rotated versions
                    float r = (float)sqrt(v_rot[0]*v_rot[0] + v_rot[1]*v_rot[1] + v_rot[2]*v_rot[2]);
                    float v_theta = acos(v_rot[2]/r) * 180.0 / PI * ARCSEC;
                    float v_phi = atan(v_rot[1]/v_rot[0]) * 180.0 / PI * ARCSEC;
                            
                    // do cut and write
                    if (v_theta > theta_cut[0] && v_theta < theta_cut[1] && 
                        v_phi > phi_cut[0] && v_phi < phi_cut[1] ) {
                        
                        // spherical coordiante transform of original (pre-rotate) position  
                        b.theta[n] = acos(b.z[n]/r) * 180.0 / PI * ARCSEC;
                        b.phi[n] = atan(b.y[n]/b.x[n]) * 180.0 / PI * ARCSEC;
                        
                        id_file.write( (char*)&b.id[n], sizeof(int64_t));
                        x_file.write((char*)&b.x[n],sizeof(float));
                        y_file.write((char*)&b.y[n],sizeof(float));
                        z_file.write((char*)&b.z[n],sizeof(float));
                        vx_file.write((char*)&b.vx[n],sizeof(float));
                        vy_file.write((char*)&b.vy[n],sizeof(float));
                        vz_file.write((char*)&b.vz[n],sizeof(float));
                        a_file.write( (char*)&b.a[n], sizeof(float));
                        rotation_file.write( (char*)&b.rotation[n], sizeof(int));
                        replication_file.write( (char*)&b.replication[n], sizeof(int32_t));
                        theta_file.write( (char*)&b.theta[n], sizeof(float));
                        phi_file.write( (char*)&b.phi[n], sizeof(float));
                        thetaRot_file.write( (char*)&v_theta, sizeof(float));
                        phiRot_file.write( (char*)&v_phi, sizeof(float));
                    }
                }
            }
        }
        
        reader.close();
        id_file.close();
        theta_file.close();
        phi_file.close();
        x_file.close();
        y_file.close();
        z_file.close();
        vx_file.close();
        vy_file.close();
        vz_file.close();
        a_file.close();
        rotation_file.close();
        replication_file.close();
        thetaRot_file.close();
        phiRot_file.close();
    }
}
