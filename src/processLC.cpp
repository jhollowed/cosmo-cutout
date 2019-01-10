#include "processLC.h"

using namespace std;
using namespace gio;


//////////////////////////////////////////////////////
//
//                Cutout function
//                  Use Case 1
//           Custom theta - phi bounds
//
//////////////////////////////////////////////////////

void processLC(string dir_name, string out_dir, vector<string> step_strings, 
               vector<float> theta_cut, vector<float> phi_cut, int myrank, int numranks,
               bool verbose, bool timeit, bool overwrite, bool positionOnly){

    ///////////////////////////////////////////////////////////////
    //
    //                          Setup
    //
    ///////////////////////////////////////////////////////////////

    // find all lc sub directories for each step in step_strings
    if(myrank == 0){ cout << "\nReading directory: " << dir_name << endl; }
    vector<string> subdirs;
    getLCSubdirs(dir_name, subdirs);
    if(myrank==0){ 
        cout << "Found subdirs:" << endl;
        for (vector<string>::const_iterator i = subdirs.begin(); i != subdirs.end(); ++i){
            cout << *i << ' ';
        }
        cout << endl;
    }

    // find the prefix (chars before the step number) in the subdirectory names.
    // It is assumed that all subdirs have the same prefix.
    string subdirPrefix;
    for(string::size_type j = 0; j < subdirs[0].size(); ++j){
        if( isdigit(subdirs[0][j]) > 0){
            subdirPrefix = subdirs[0].substr(0, j);
            break;
        }
    }
    if(myrank == 0){ cout << "Subdir prefix is: " << subdirPrefix << endl; }

    ///////////////////////////////////////////////////////////////
    //
    //                 Loop over step subdirs
    //
    ///////////////////////////////////////////////////////////////


    // perform cutout on data from each lc output step
    size_t max_size = 0;
    int step;
    for (int i=0; i<step_strings.size();++i){
        
        // instances of buffer structs at file header
        Buffers_read r;
        Buffers_write w;

        // continue if this is the step at z=0 (lightcone volume zero)
        step =atoi(step_strings[i].c_str());
        if(step == 499){ continue;}

        // find header file
        if(myrank == 0){
            cout<< "\n---------- Working on step "<< step_strings[i] <<"----------" << endl; 
        }
        string file_name;
        ostringstream file_name_stream;
        file_name_stream << dir_name << subdirPrefix << step_strings[i];

        getLCFile(file_name_stream.str(), file_name);
        file_name_stream << "/" << file_name; 

        // setup gio
        size_t Np = 0;
        unsigned Method = GenericIO::FileIOPOSIX;
        const char *EnvStr = getenv("GENERICIO_USE_MPIIO");
        if(EnvStr && string(EnvStr) == "1"){
            Method = GenericIO::FileIOMPI;  
        }

        // create gio reader, open lightcone file header in new scope
        {
            if(myrank == 0){ cout << "Opening file: " << file_name_stream.str() << endl; }
            GenericIO GIO(MPI_COMM_WORLD, file_name_stream.str(), Method);
            GIO.openAndReadHeader(GenericIO::MismatchRedistribute);

            MPI_Barrier(MPI_COMM_WORLD);
            Np = GIO.readNumElems();
            if(myrank == 0){
                cout << "Number of elements in lc step at rank " << myrank << ": " << 
                Np << endl; 
            }

            // resize buffers   
            r.x.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
            r.y.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
            r.z.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
            r.vx.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
            r.vy.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
            r.vz.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
            r.a.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
            r.id.resize(Np + GIO.requestedExtraSpace()/sizeof(ID_T));
            r.rotation.resize(Np + GIO.requestedExtraSpace()/sizeof(int));
            r.replication.resize(Np + GIO.requestedExtraSpace()/sizeof(int32_t));

            // do reading
            GIO.addVariable("x", r.x, true); 
            GIO.addVariable("y", r.y, true); 
            GIO.addVariable("z", r.z, true); 
            GIO.addVariable("vx", r.vx, true); 
            GIO.addVariable("vy", r.vy, true); 
            GIO.addVariable("vz", r.vz, true); 
            GIO.addVariable("a", r.a, true); 
            GIO.addVariable("id", r.id, true); 
            GIO.addVariable("rotation", r.rotation, true); 
            GIO.addVariable("replication", r.replication, true);

            GIO.readData(); 
        }

        // resize again to remove reader extra space
        r.x.resize(Np);
        r.y.resize(Np);
        r.z.resize(Np);
        r.vx.resize(Np);
        r.vy.resize(Np);
        r.vz.resize(Np);
        r.a.resize(Np);
        r.id.resize(Np);
        r.rotation.resize(Np);
        r.replication.resize(Np);
        if(myrank == 0){ cout<<"done resizing"<<endl; }

        ///////////////////////////////////////////////////////////////
        //
        //           Create output files + start reading
        //
        ///////////////////////////////////////////////////////////////
        
        // open cutout subdirectory for this step...
        ostringstream step_subdir;
        step_subdir << out_dir << subdirPrefix << "Cutout" << step_strings[i];
        
        DIR *dir = opendir(step_subdir.str().c_str());
        struct dirent *d;
        int nf = 0;
        
        // if subdir already exists, make sure it's empty, because overwriting
        // binary files isn't always clean
        if(dir != NULL){
            while((d = readdir(dir)) != NULL){ if(++nf>2){ break;} }
            closedir(dir);

            if(nf > 2){
                cout << "\nDirectory " << step_subdir.str() << " is non-empty" << endl;
                MPI_Abort(MPI_COMM_WORLD, 0);
            }
            if(myrank == 0){ cout << "Entered subdir: " << step_subdir.str() << endl; }
        }
        // Otherwise, create the subdir
        else{
            mkdir(step_subdir.str().c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IXOTH);
            if(myrank == 0){ cout << "Created subdir: " << step_subdir.str() << endl; }
        }

        // create binary files for cutout output
        MPI_File id_file;
        MPI_File redshift_file;
        MPI_File x_file;
        MPI_File y_file;
        MPI_File z_file;
        MPI_File vx_file;
        MPI_File vy_file;
        MPI_File vz_file;
        MPI_File rotation_file;
        MPI_File replication_file;
        MPI_File theta_file;
        MPI_File phi_file;

        MPI_Request id_req;
        MPI_Request redshift_req;
        MPI_Request x_req;
        MPI_Request y_req;
        MPI_Request z_req;
        MPI_Request vx_req;
        MPI_Request vy_req;
        MPI_Request vz_req;
        MPI_Request rotation_req;
        MPI_Request replication_req;
        MPI_Request theta_req;
        MPI_Request phi_req;
        
        ostringstream id_file_name;
        ostringstream redshift_file_name;
        ostringstream x_file_name;
        ostringstream y_file_name;
        ostringstream z_file_name;
        ostringstream vx_file_name;
        ostringstream vy_file_name;
        ostringstream vz_file_name;
        ostringstream rotation_file_name;
        ostringstream replication_file_name;
        ostringstream theta_file_name;
        ostringstream phi_file_name; 

        id_file_name << step_subdir.str() << "/id." << step << ".bin";
        redshift_file_name << step_subdir.str() << "/redshift." << step << ".bin";
        x_file_name << step_subdir.str() << "/x."<< step <<".bin";
        y_file_name << step_subdir.str() << "/y."<< step <<".bin";
        z_file_name << step_subdir.str() << "/z."<< step <<".bin";
        vx_file_name << step_subdir.str() << "/vx."<< step <<".bin";
        vy_file_name << step_subdir.str() << "/vy."<< step <<".bin";
        vz_file_name << step_subdir.str() << "/vz."<< step <<".bin";
        rotation_file_name << step_subdir.str() << "/rotation."<< step <<".bin";
        replication_file_name << step_subdir.str() << "/replication."<< step <<".bin";
        theta_file_name << step_subdir.str() << "/theta." << step << ".bin";
        phi_file_name << step_subdir.str() << "/phi." << step << ".bin";

        if(myrank == 0){ cout<<"starting to open files"<<endl; }

        MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(id_file_name.str().c_str()), 
                MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &id_file);
        MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(x_file_name.str().c_str()), 
                MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &x_file);
        MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(y_file_name.str().c_str()),
                MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &y_file);
        MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(z_file_name.str().c_str()),
                MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &z_file);
        MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(vx_file_name.str().c_str()),
                MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &vx_file);
        MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(vy_file_name.str().c_str()),
                MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &vy_file);
        MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(vz_file_name.str().c_str()),
                MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &vz_file);
        MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(redshift_file_name.str().c_str()),
                MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &redshift_file);
        MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(rotation_file_name.str().c_str()),
                MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &rotation_file);
        MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(replication_file_name.str().c_str()),
                MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &replication_file);
        MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(theta_file_name.str().c_str()),
                MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &theta_file);
        MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(phi_file_name.str().c_str()),
                MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &phi_file);
        
        if(myrank == 0){ cout<<"done opening files"<<endl; }

        ///////////////////////////////////////////////////////////////
        //
        //                         Do cutting
        //
        ///////////////////////////////////////////////////////////////

        if(myrank == 0){ cout << "Converting positions..." << endl; }

        for (int n=0; n<Np; ++n) {

            // limit cutout to first octant for speed
            if (r.x[n] > 0.0 && r.y[n] > 0.0 && r.z[n] > 0.0){

                // spherical coordinate transformation
                float d = (float)sqrt(r.x[n]*r.x[n] + r.y[n]*r.y[n] + r.z[n]*r.z[n]);
                float theta = acos(r.z[n]/d) * 180.0 / PI * ARCSEC;
                float phi = atan(r.y[n]/r.x[n]) * 180.0 / PI * ARCSEC;

                // do cut and write
                if (theta > theta_cut[0] && theta < theta_cut[1] && 
                        phi > phi_cut[0] && phi < phi_cut[1] ) {

                    // get redshift from scale factor
                    float zz = aToZ(r.a[n]);  

                    // spherical corrdinate transform of positions
                    w.theta.push_back(theta);
                    w.phi.push_back(phi);

                    // other columns
                    w.x.push_back(r.x[n]);
                    w.y.push_back(r.y[n]);
                    w.z.push_back(r.z[n]);
                    w.vx.push_back(r.vx[n]);
                    w.vy.push_back(r.vy[n]);
                    w.vz.push_back(r.vz[n]);
                    w.redshift.push_back(zz);
                    w.id.push_back(r.id[n]);
                    w.rotation.push_back(r.rotation[n]);
                    w.replication.push_back(r.replication[n]);
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);

        ///////////////////////////////////////////////////////////////
        //
        //                   write out
        //
        ///////////////////////////////////////////////////////////////

        // define MPI file writing offset for the current rank --
        // This offset will be the sum of elements in all lesser ranks,
        // multiplied by the type size for each file    
        w.np_count.clear();
        w.np_count.resize(numranks);
        w.np_offset.clear();
        w.np_offset[0] = 0;
        

        int cutout_size = int(w.redshift.size());
        
        // get number of elements in each ranks portion of cutout
        MPI_Allgather(&cutout_size, 1, MPI_INT, &w.np_count[0], 1, MPI_INT, 
                      MPI_COMM_WORLD);
        
        // compute each ranks writing offset
        for(int j=1; j < numranks; ++j){
            w.np_offset.push_back(w.np_offset[j-1] + w.np_count[j-1]);
        }
        
        // print out offset vector for verification
        if(myrank == 0){
            if(numranks < 20){
                cout << "rank object counts: [";
                for(int m=0; m < numranks; ++m){ cout << w.np_count[m] << ","; }
                cout << "]" << endl;
                cout << "rank offsets: [";
                for(int m=0; m < numranks; ++m){ cout << w.np_offset[m] << ","; }
                cout << "]" << endl;
            } else {
               int numEmpty = count(w.np_count.begin(), w.np_count.end(), 0);
               cout << numranks - numEmpty << " of " << numranks << 
               " ranks found members within cutout field of view" << endl;
            }
        }

        MPI_Offset offset_posvel = sizeof(POSVEL_T) * w.np_offset[myrank];
        MPI_Offset offset_id = sizeof(ID_T) * w.np_offset[myrank];
        MPI_Offset offset_float = sizeof(float) * w.np_offset[myrank];
        MPI_Offset offset_int = sizeof(int) * w.np_offset[myrank];
        MPI_Offset offset_int32 = sizeof(int32_t) * w.np_offset[myrank];

        // write
        MPI_File_seek(id_file, offset_id, MPI_SEEK_SET);
        MPI_File_iwrite(id_file, &w.id[0], w.id.size(), MPI_INT64_T, &id_req);
        MPI_Wait(&id_req, MPI_STATUS_IGNORE);

        MPI_File_seek(x_file, offset_posvel, MPI_SEEK_SET);
        MPI_File_iwrite(x_file, &w.x[0], w.x.size(), MPI_FLOAT, &x_req);
        MPI_Wait(&x_req, MPI_STATUS_IGNORE);

        MPI_File_seek(y_file, offset_posvel, MPI_SEEK_SET);
        MPI_File_iwrite(y_file, &w.y[0], w.y.size(), MPI_FLOAT, &y_req);
        MPI_Wait(&y_req, MPI_STATUS_IGNORE);
        
        MPI_File_seek(z_file, offset_posvel, MPI_SEEK_SET);
        MPI_File_iwrite(z_file, &w.z[0], w.z.size(), MPI_FLOAT, &z_req);
        MPI_Wait(&z_req, MPI_STATUS_IGNORE);

        MPI_File_seek(vx_file, offset_posvel, MPI_SEEK_SET);
        MPI_File_iwrite(vx_file, &w.vx[0], w.vx.size(), MPI_FLOAT, &vx_req);
        MPI_Wait(&vx_req, MPI_STATUS_IGNORE);

        MPI_File_seek(vy_file, offset_posvel, MPI_SEEK_SET);
        MPI_File_iwrite(vy_file, &w.vy[0], w.vy.size(), MPI_FLOAT, &vy_req);
        MPI_Wait(&vy_req, MPI_STATUS_IGNORE);
        
        MPI_File_seek(vz_file, offset_posvel, MPI_SEEK_SET);
        MPI_File_iwrite(vz_file, &w.vz[0], w.vz.size(), MPI_FLOAT, &vz_req);
        MPI_Wait(&vz_req, MPI_STATUS_IGNORE);
        
        MPI_File_seek(theta_file, offset_posvel, MPI_SEEK_SET);
        MPI_File_iwrite(theta_file, &w.theta[0], w.theta.size(), MPI_FLOAT, &theta_req);
        MPI_Wait(&theta_req, MPI_STATUS_IGNORE);
        
        MPI_File_seek(phi_file, offset_posvel, MPI_SEEK_SET);
        MPI_File_iwrite(phi_file, &w.phi[0], w.phi.size(), MPI_FLOAT, &phi_req);
        MPI_Wait(&phi_req, MPI_STATUS_IGNORE);
        
        MPI_File_seek(redshift_file, offset_posvel, MPI_SEEK_SET);
        MPI_File_iwrite(redshift_file, &w.redshift[0], w.redshift.size(), MPI_FLOAT, &redshift_req);
        MPI_Wait(&redshift_req, MPI_STATUS_IGNORE);
        
        MPI_File_seek(id_file, offset_posvel, MPI_SEEK_SET);
        MPI_File_iwrite(id_file, &w.id[0], w.id.size(), MPI_FLOAT, &id_req);
        MPI_Wait(&id_req, MPI_STATUS_IGNORE);
        
        MPI_File_seek(rotation_file, offset_posvel, MPI_SEEK_SET);
        MPI_File_iwrite(rotation_file, &w.rotation[0], w.rotation.size(), 
                        MPI_FLOAT, &rotation_req);
        MPI_Wait(&rotation_req, MPI_STATUS_IGNORE);
        
        MPI_File_seek(replication_file, offset_posvel, MPI_SEEK_SET);
        MPI_File_iwrite(replication_file, &w.replication[0], w.replication.size(), 
                        MPI_FLOAT, &replication_req);
        MPI_Wait(&replication_req, MPI_STATUS_IGNORE);

        MPI_File_close(&id_file);
        MPI_File_close(&x_file);
        MPI_File_close(&y_file);
        MPI_File_close(&z_file);
        MPI_File_close(&vx_file);
        MPI_File_close(&vy_file);
        MPI_File_close(&vz_file);
        MPI_File_close(&redshift_file);
        MPI_File_close(&theta_file);
        MPI_File_close(&phi_file);
        MPI_File_close(&rotation_file);
        MPI_File_close(&replication_file);
    }
}


//////////////////////////////////////////////////////
//
//                Cutout function
//                  Use Case 2
//               Custom halo cutout
//
//////////////////////////////////////////////////////

void processLC(string dir_name, vector<string> out_dirs, vector<string> step_strings, 
               vector<float> halo_pos, vector<float> halo_props, float boxLength, int myrank, 
               int numranks, bool verbose, bool timeit, bool overwrite, bool positionOnly){

    ///////////////////////////////////////////////////////////////
    //
    //                          Setup
    //
    ///////////////////////////////////////////////////////////////

    // find all lc sub directories for each step in step_strings
    if(myrank == 0){ cout << "\nReading directory: " << dir_name << endl; }
    vector<string> subdirs;
    getLCSubdirs(dir_name, subdirs);
    if(myrank==0){ 
        cout << "Found subdirs:" << endl;
        for (vector<string>::const_iterator i = subdirs.begin(); i != subdirs.end(); ++i){
            cout << *i << ' ';
        }
        cout << endl;
    }

    // find the prefix (chars before the step number) in the subdirectory names.
    // It is assumed that all subdirs have the same prefix.
    string subdirPrefix;
    for(string::size_type j = 0; j < subdirs[0].size(); ++j){
        if( isdigit(subdirs[0][j]) > 0){
            subdirPrefix = subdirs[0].substr(0, j);
            break;
        }
    }
    if(myrank == 0){ cout << "Subdir prefix is: " << subdirPrefix << endl; }

    ///////////////////////////////////////////////////////////////
    //
    //                  Start coordinate rotation
    //
    ///////////////////////////////////////////////////////////////

    if(myrank == 0){
        cout<< "\n\n---------- Setting up for coordinate rotation ----------" << endl; 
    }

    // do coordinate rotation to center each input halo at (r, 90, 0) in spherical coords...
    // halo_pos is a vector of halo positions, each having three components
    // k is a vector of vectors (the rotation axis per target halo)
    // B is a vector (the rotation angle per target halo)
    // theta is a vector of vectors (the altitude angular bounds per target halo)
    // phi is a vector of vectors (the azimuthal angular bounds per target halo)
    
    int numHalos = halo_pos.size()/3;
    bool printHalo;
    
    // rotation axis k and angle B per halo
    vector<vector<float> > k(numHalos);
    vector<float> B(numHalos);

    // cross-product matrix K and rotation matrix R per halo
    vector<vector<vector<float> > > K(numHalos);
    vector<vector<vector<float> > > R(numHalos);
    vector<vector<vector<float> > > R_inv(numHalos);
    
    // constant (equitorial) angular bounds per halo and
    // rough constant (rotated) angular bounds per halo
    vector<vector<float> > theta_cut(numHalos);
    vector<vector<float> > phi_cut(numHalos);
    vector<vector<float> > theta_cut_rough(numHalos);
    vector<vector<float> > phi_cut_rough(numHalos);

    for(int h=0; h<halo_pos.size(); h+=3){ 
        
        int haloIdx = h/3;
        printHalo = (numHalos < 20) | (haloIdx%20==0) ? 1:0;
        
        // get next three values in halo_pos
        float tmp_pos[] = {halo_pos[h], halo_pos[h+1], halo_pos[h+2]};
        vector<float> this_halo_pos(tmp_pos, tmp_pos+3);

        // get next three or four values in halo_props (redshift, step, sod_mass, sod_radius) or 
        // (redshift, step, fof_mass) 
        // last element of tmp_props will be garbage if numProps == 3 (massDef == 'fof'),
        // but that's okay. 
        int numProps = halo_props.size() / numHalos;
        float tmp_props[] = {halo_props[numProps*haloIdx], halo_props[numProps*haloIdx+1], 
                             halo_props[numProps*haloIdx+2], halo_props[numProps*haloIdx+3]};
        if(numProps != 3 and numProps != 4){
            cout << "Something went wrong... " << numProps << " properties found per halo" << 
                    " in input vectors. Is massDef correct? Please contact developer." << endl;
            MPI_Finalize();
            exit(EXIT_FAILURE);
        } 
        vector<float> this_halo_props(tmp_props, tmp_props+numProps);
        
        // find distance magnitude (new rotated halo position)
        float halo_r = (float)sqrt(this_halo_pos[0]*this_halo_pos[0] + 
                                   this_halo_pos[1]*this_halo_pos[1] + 
                                   this_halo_pos[2]*this_halo_pos[2]);
        
        // get the halo-centric angular bounds of the cutout...
        // calculate dtheta and dphi in radians gven boxlength in arcmin
        float halfBoxLength = ((boxLength/2.0) / 60) * PI/180.0;
        float dtheta = halfBoxLength;
        float dphi = dtheta;

        // calculate theta_cut and phi_cut, in arcsec, given the specified boxLength
        theta_cut[haloIdx].push_back( (PI/2 - dtheta) * 180.0/PI * ARCSEC );
        theta_cut[haloIdx].push_back( (PI/2 + dtheta) * 180.0/PI * ARCSEC );
        phi_cut[haloIdx].push_back( (0 - dphi) * 180.0/PI * ARCSEC );
        phi_cut[haloIdx].push_back( (0 + dphi) * 180.0/PI * ARCSEC );
        
        if(myrank == 0 and printHalo){
            cout << "\n--- Target halo " << haloIdx << " ---" << endl; 
            cout << "theta bounds set to: ";
            cout << theta_cut[haloIdx][0]/ARCSEC << "deg -> " << theta_cut[haloIdx][1]/ARCSEC <<"deg"<< endl;
            cout << "phi bounds set to: ";
            cout << phi_cut[haloIdx][0]/ARCSEC << "deg -> " << phi_cut[haloIdx][1]/ARCSEC <<"deg" << endl;
            cout << "theta-phi bounds result in box width of " << 
                    tan(dtheta) * halo_r * 2 << 
                    " Mpc/h at distance to halo of " << halo_r << " Mpc/h" << endl << 
                    "        " << "= " << dtheta*2*180.0/PI << "deg x " << dphi*2*180.0/PI << 
                    "deg field of view" << endl;
        } 
        
        // Now let's rotate the angular bounds to their true position in the sky
        // Using the Rodrigues rotation formula...
        float tmp_rot_pos[] = {halo_r, 0, 0};
        vector<float> rotated_pos(tmp_rot_pos, tmp_rot_pos+3);
        if(myrank == 0 and printHalo){
            cout << "\nFinding axis of rotation to move (" << 
                    this_halo_pos[0]<< ", " << this_halo_pos[1]<< ", " << this_halo_pos[2]<< ") to (" <<
                    rotated_pos[0] << ", " << rotated_pos[1] << ", " << rotated_pos[2] <<
                    ")" << endl; 
        }

        // get angle and axis of rotation
        normCross(this_halo_pos, rotated_pos, k[haloIdx]);
        B[haloIdx] = vecPairAngle(this_halo_pos, rotated_pos);
        if(myrank == 0 and printHalo){
            cout << "Rotation is " << B[haloIdx]*(180/PI) << "deg about axis k = (" << 
                    k[haloIdx][0]<< ", " << k[haloIdx][1] << ", " << k[haloIdx][2] << ")" << endl;
        }
        
        // get rotation matrix R
        cross_prod_matrix(k[haloIdx], K[haloIdx]);
        rotation_matrix(K[haloIdx], B[haloIdx], R[haloIdx]);

        // invert rotation matrix R
        R_inv[haloIdx] = invert_3x3(R[haloIdx]);
         
        // verbose output
        if(myrank == 0 and verbose == true and printHalo){
            cout << "\nRotation Matrix is " << endl << 
                    "{ " << R[haloIdx][0][0] << ", " << R[haloIdx][0][1] << ", " << 
                            R[haloIdx][0][2] << "}" << endl <<
                    "{ " << R[haloIdx][1][0] << ", " << R[haloIdx][1][1] << ", " << 
                            R[haloIdx][1][2] << "}" << endl <<
                    "{ " << R[haloIdx][2][0] << ", " << R[haloIdx][2][1] << ", " << 
                            R[haloIdx][2][2] << "}" << endl;
        }
        if(myrank == 0 and verbose == true and printHalo){
            cout << "Inverted Rotation Matrix is " << endl << 
                    "{ " << R_inv[haloIdx][0][0] << ", " << R_inv[haloIdx][0][1] << ", " << 
                            R_inv[haloIdx][0][2] << "}" << endl <<
                    "{ " << R_inv[haloIdx][1][0] << ", " << R_inv[haloIdx][1][1] << ", " << 
                            R_inv[haloIdx][1][2] << "}" << endl <<
                    "{ " << R_inv[haloIdx][2][0] << ", " << R_inv[haloIdx][2][1] << ", " << 
                            R_inv[haloIdx][2][2] << "}" << endl;
        }

        // now, we have defined theta_cut and phi_cut above in such a way that it 
        // encapsulates a field of angular size boxLength for the target halo
        // *if it were lying on the equator*. With the cutout sized correctly, we can now
        // "point it" at the true halo posiiton to get the real non-constant angular bounds.
        // We do this by defining four unit vectors (A, B, C, D) directed toward the 
        // corners of the bounded fov on the equator, and rotate them according to the 
        // inverted rotation matrix defined above:
        //                                          C
        //                                        _--_ 
        //    B-----------C                     _-    -_  
        //    |           |    R_inv          _-        -_   
        //    |           |  -------->>   B _-            -_ D
        //    |   fov     |                  -_   fov    _-
        //    |           |                    -_      _-
        //    A-----------D                      -_  _-  <--- boxLength
        //          ^                              --
        //          |__ boxLength                  A
        //
      
        // get theta and phi bounds in radians 
        vector<float> theta_cut_rad;
        theta_cut_rad.push_back( (theta_cut[haloIdx][0] / ARCSEC) * PI/180.0); 
        theta_cut_rad.push_back( (theta_cut[haloIdx][1] / ARCSEC) * PI/180.0); 
        vector<float> phi_cut_rad;
        phi_cut_rad.push_back( (phi_cut[haloIdx][0] / ARCSEC) * PI/180.0); 
        phi_cut_rad.push_back( (phi_cut[haloIdx][1] / ARCSEC) * PI/180.0); 
        
        float fov_corner_A[] = { sin(theta_cut_rad[1]) * cos(phi_cut_rad[1]), 
                                 sin(theta_cut_rad[1]) * sin(phi_cut_rad[1]), 
                                 cos(theta_cut_rad[1]) };
        vector<float> A(fov_corner_A, fov_corner_A+3);
        vector<float> A_rot;  

        float fov_corner_B[] = { sin(theta_cut_rad[0]) * cos(phi_cut_rad[1]), 
                                 sin(theta_cut_rad[0]) * sin(phi_cut_rad[1]), 
                                 cos(theta_cut_rad[0]) };
        vector<float> B(fov_corner_B, fov_corner_B+3);
        vector<float> B_rot;  

        float fov_corner_C[] = { sin(theta_cut_rad[0]) * cos(phi_cut_rad[0]), 
                                 sin(theta_cut_rad[0]) * sin(phi_cut_rad[0]), 
                                 cos(theta_cut_rad[0]) };
        vector<float> C(fov_corner_C, fov_corner_C+3);
        vector<float> C_rot;  
        
        float fov_corner_D[] = { sin(theta_cut_rad[1]) * cos(phi_cut_rad[0]), 
                                 sin(theta_cut_rad[1]) * sin(phi_cut_rad[0]), 
                                 cos(theta_cut_rad[1]) };
        vector<float> D(fov_corner_D, fov_corner_D+3);
        vector<float> D_rot;  
        
        A_rot = matVecMul(R_inv[haloIdx], A);
        B_rot = matVecMul(R_inv[haloIdx], B); 
        C_rot = matVecMul(R_inv[haloIdx], C); 
        D_rot = matVecMul(R_inv[haloIdx], D); 

        if(myrank == 0 and verbose == true and printHalo){ 
                               cout << "\nRotated cartesian vectors pointing toward FOV corners are" << endl <<
                               "A = {" << A_rot[0] << ", " << A_rot[1] << ", " << A_rot[2] << "}" << endl <<
                               "B = {" << B_rot[0] << ", " << B_rot[1] << ", " << B_rot[2] << "}" << endl <<
                               "C = {" << C_rot[0] << ", " << C_rot[1] << ", " << C_rot[2] << "}" << endl <<
                               "D = {" << D_rot[0] << ", " << D_rot[1] << ", " << D_rot[2] << "}" << endl;
        }

        // convert vectors A,B,C,D to 2-dimensional spherical coordinate 
        // vectors {theta, phi} in arcsec 
        vector<float> A_sph;
        A_sph.push_back( acos(A_rot[2]/1) * 180.0/PI * ARCSEC);
        A_sph.push_back( atan(A_rot[1]/A_rot[0]) * 180.0/PI * ARCSEC);
        
        vector<float> B_sph;
        B_sph.push_back( acos(B_rot[2]/1) * 180.0/PI * ARCSEC);
        B_sph.push_back( atan(B_rot[1]/B_rot[0]) * 180.0/PI * ARCSEC);
        
        vector<float> C_sph;
        C_sph.push_back( acos(C_rot[2]/1) * 180.0/PI * ARCSEC);
        C_sph.push_back( atan(C_rot[1]/C_rot[0]) * 180.0/PI * ARCSEC);
        
        vector<float> D_sph;
        D_sph.push_back( acos(D_rot[2]/1) * 180.0/PI * ARCSEC);
        D_sph.push_back( atan(D_rot[1]/D_rot[0]) * 180.0/PI * ARCSEC);
        
        if(myrank == 0 and verbose == true and printHalo){ 
                               cout << "\nAngular positions of the FOV corners in degrees are" << endl <<
                               "A = {" << A_sph[0]/ARCSEC << ", " << A_sph[1]/ARCSEC <<
                                          "}" << endl <<
                               "B = {" << B_sph[0]/ARCSEC << ", " << B_sph[1]/ARCSEC <<
                                          "}" << endl <<
                               "C = {" << C_sph[0]/ARCSEC << ", " << C_sph[1]/ARCSEC <<
                                          "}" << endl <<
                               "D = {" << D_sph[0]/ARCSEC << ", " << D_sph[1]/ARCSEC <<
                                          "}" << endl;
        }

        // define rough-cut bounds, in arcsec, to be used to quickly remove particles that certainly 
        // are not in the field of view. After a cut is done in this way, we can treat the remaining 
        // particles more carefully. We do this by simply cutting on the max and min theta and phi 
        // values among the vectors A,B,and C, with a buffer of constant size 10 arcmin (where typical
        // cluster cutouts at modest redshifts come out to have a width of order 1 degree)

        float ang_buffer = 600; // buffer in arcsec

        float tmp_thetas[] = {A_sph[0], B_sph[0], C_sph[0], D_sph[0]};
        vector<float> theta_corners(tmp_thetas, tmp_thetas+4);
        theta_cut_rough[haloIdx].push_back( *min_element(theta_corners.begin(), theta_corners.end()) - ang_buffer );
        theta_cut_rough[haloIdx].push_back( *max_element(theta_corners.begin(), theta_corners.end()) + ang_buffer ); 
        
        float tmp_phis[] = {A_sph[1], B_sph[1], C_sph[1], D_sph[1]};
        vector<float> phi_corners(tmp_phis, tmp_phis+4);
        phi_cut_rough[haloIdx].push_back( *min_element(phi_corners.begin(), phi_corners.end()) - ang_buffer );
        phi_cut_rough[haloIdx].push_back( *max_element(phi_corners.begin(), phi_corners.end()) + ang_buffer ); 
        
        if(myrank == 0 and printHalo){
            cout << "\nrough theta bounds set to: ";
            cout << theta_cut_rough[haloIdx][0]/ARCSEC << "deg -> " << 
                    theta_cut_rough[haloIdx][1]/ARCSEC <<"deg"<< endl;
            cout << "rough phi bounds set to: ";
            cout << phi_cut_rough[haloIdx][0]/ARCSEC << "deg -> " << 
                    phi_cut_rough[haloIdx][1]/ARCSEC <<"deg" << endl;
        }

        // Write out a csv file to this halo's cutout directory to contain halo properties, 
        // cutout specifications,  and run meta data
        
        ofstream props_file;
        ostringstream props_file_name; 
        props_file_name << out_dirs[haloIdx]<< "/properties.csv";
        props_file.open(props_file_name.str().c_str());
        
        props_file << "#halo_redshift" << ", " << "halo_lc_shell" << ", ";
        if(numProps == 4)
            props_file << "sod_halo_mass" << ", " <<  "sod_halo_radius";
        else
            props_file << "fof_halo_mass";
        props_file << ", " << "halo_lc_x" << ", " << "halo_lc_y" << ", " << "halo_lc_z" << ", "
                              "boxRadius_Mpc" << ", " << "boxRadius_arcsec" << "\n";
        for(int i=0; i<this_halo_props.size(); ++i)
            props_file << this_halo_props[i] << ", ";
        for(int i=0; i<this_halo_pos.size(); ++i)
            props_file << this_halo_pos[i] << ", ";
       props_file << atan(halfBoxLength) * halo_r << ", " << halfBoxLength * 180.0/PI * ARCSEC << "\n";
       props_file.close();
    }

    ///////////////////////////////////////////////////////////////
    //
    //                 Loop over step subdirs
    //
    ///////////////////////////////////////////////////////////////

    MPI_Barrier(MPI_COMM_WORLD);

    // perform cutout on data from each lc output step
    size_t max_size = 0;
    int step;
 
    MPI_Datatype particles_mpi_pos = createParticles_pos();
    MPI_Datatype particles_mpi_vel = createParticles_vel();

    vector<double> read_times;
    vector<double> redist_times;
    vector<double> cutout_times; 
    vector<double> write_times; 
    double start;
    double stop;
    double duration;
    
    for (int i=0; i<step_strings.size(); ++i){
   
        // time read in 
        MPI_Barrier(MPI_COMM_WORLD);
        start = MPI_Wtime();
        
        // instances of buffer struct at file header for read in data
        Buffers_read r;
        
        // continue if this is the step at z=0 (lightcone volume zero)
        step =atoi(step_strings[i].c_str());
        if(step == 499){ continue;}

        // find header file
        if(myrank == 0){
            cout << "\n=================================================" << endl;
            cout << "============== Working on step "<< step_strings[i] <<" ==============\n" << endl; 
        }
        string file_name;
        ostringstream file_name_stream;
        file_name_stream << dir_name << subdirPrefix << step_strings[i];

        getLCFile(file_name_stream.str(), file_name);
        file_name_stream << "/" << file_name; 

        ///////////////////////////////////////////////////////////////
        //
        //                        do reading
        //
        ///////////////////////////////////////////////////////////////
        
        // setup gio
        size_t Np = 0;
        unsigned Method = GenericIO::FileIOPOSIX;
        const char *EnvStr = getenv("GENERICIO_USE_MPIIO");
        if(EnvStr && string(EnvStr) == "1"){
            Method = GenericIO::FileIOMPI;  
        }

        // create gio reader, open lightcone file header in new scope
        {
            if(myrank == 0){ cout << "Opening file: " << file_name_stream.str() << endl; }
            GenericIO GIO(MPI_COMM_WORLD, file_name_stream.str(), Method);
            GIO.openAndReadHeader(GenericIO::MismatchRedistribute);

            MPI_Barrier(MPI_COMM_WORLD);
            Np = GIO.readNumElems();
           
            // resize buffers   
            r.x.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
            r.y.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
            r.z.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
            r.a.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
            r.id.resize(Np + GIO.requestedExtraSpace()/sizeof(ID_T));
            if(!positionOnly){
                r.vx.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
                r.vy.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
                r.vz.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
                r.rotation.resize(Np + GIO.requestedExtraSpace()/sizeof(int));
                r.replication.resize(Np + GIO.requestedExtraSpace()/sizeof(int32_t));
            }

            // do reading
            GIO.addVariable("x", r.x, true); 
            GIO.addVariable("y", r.y, true); 
            GIO.addVariable("z", r.z, true); 
            GIO.addVariable("a", r.a, true); 
            GIO.addVariable("id", r.id, true); 
            if(!positionOnly){
                GIO.addVariable("vx", r.vx, true); 
                GIO.addVariable("vy", r.vy, true); 
                GIO.addVariable("vz", r.vz, true); 
                GIO.addVariable("rotation", r.rotation, true); 
                GIO.addVariable("replication", r.replication, true);
            }

            GIO.readData(); 
        }

        // resize again to remove reader extra space
        r.x.resize(Np);
        r.y.resize(Np);
        r.z.resize(Np);
        r.a.resize(Np);
        r.id.resize(Np);
        if(!positionOnly){
            r.vx.resize(Np);
            r.vy.resize(Np);
            r.vz.resize(Np);
            r.rotation.resize(Np);
            r.replication.resize(Np);
        }
        if(myrank == 0){ cout<<"done resizing"<<endl; }
        
        // calc d, theta, and phi per particle
        r.d.resize(Np);
        r.theta.resize(Np);
        r.phi.resize(Np);
                
        for (int n=0; n<Np; ++n) {
            // spherical coordinate transformation
            r.d[n] = (float)sqrt( r.x[n]*r.x[n] + r.y[n]*r.y[n] + r.z[n]*r.z[n]);
            r.theta[n] = acos(r.z[n]/r.d[n]) * 180.0 / PI * ARCSEC;
            
            // prevent NaNs on y-z plane
            if(r.x[n] == 0 && r.y[n] > 0)
                r.phi[n] = 90.0 * ARCSEC;
            else if(r.x[n] == 0 && r.y[n] < 0)
                r.phi[n] = -90.0 * ARCSEC;
            else
                r.phi[n] = atan(r.y[n]/r.x[n]) * 180.0 / PI * ARCSEC;
        }

        MPI_Barrier(MPI_COMM_WORLD);
        stop = MPI_Wtime();
    
        duration = stop - start;
        if(myrank == 0 and timeit == true){ cout << "Read time: " << duration << " s" << endl; }
        read_times.push_back(duration);
         
        ///////////////////////////////////////////////////////////////
        //
        //           evenly distribute particles to all ranks
        //
        ///////////////////////////////////////////////////////////////
        
        // now, we have likely ended up in the situation where most of the data read
        // resides in only a small subset of the ranks in this communicator. This is becuase
        // (I think) only a small fraction of the number of ranks that generated the lightcone 
        // at any particular step intersected that lightcone shell. Tthat information is encoded 
        // in the lightcone output in the fact that all ranks which found no particles interesecting
        // the shell created an empty block in the resultant GIO files.
        // That's no good, because for this cutout code, we don't want some few nodes to be doing 
        // lots of computation while most others do none, so let's scatter the input data evenly 
        // across all ranks
     
        // time redistribution 
        MPI_Barrier(MPI_COMM_WORLD);
        start = MPI_Wtime();
         
        // find number of empty ranks
        vector<size_t> Np_read_per_rank(numranks); 
        MPI_Allgather(&Np, 1, MPI_INT64_T, &Np_read_per_rank[0], 1, MPI_INT64_T, 
                      MPI_COMM_WORLD);
        int num_readNone = count(&Np_read_per_rank[0], &Np_read_per_rank[numranks], 0); 
        
        // and number of total particles per rank with data
        size_t totalNp = 0;
        for(int ri = 0; ri < numranks; ++ri){
            totalNp += Np_read_per_rank[ri];
        }

        if(myrank == 0){
            cout << "Total number of particles is " << totalNp << endl;
            cout << "Redistributing particles to all from " << numranks - num_readNone << 
                    " of " << numranks << " ranks" << endl;
        }   
         
        vector<int> even_redistribute;
        vector<int> redist_send_count(numranks);
        vector<int> redist_recv_count(numranks);
        vector<int> redist_send_offset(numranks);
        vector<int> redist_recv_offset(numranks);
        
        // compute number of particles to send to each other rank
        comp_rank_scatter(Np, even_redistribute, numranks);
        for(int ri = 0; ri < numranks; ++ri){
            redist_send_count[ri] = count(&even_redistribute[0], &even_redistribute[Np], ri);
        }
        
        // get number of particles to recieve from every other rank
        MPI_Alltoall(&redist_send_count[0], 1, MPI_INT, &redist_recv_count[0], 1, MPI_INT, MPI_COMM_WORLD);

        // compute sending+recieving offsets to/from each other rank
        for(int ri=1; ri < numranks; ++ri){
            redist_send_offset[ri] = redist_send_offset[ri-1] + redist_send_count[ri-1];
            redist_recv_offset[ri] = redist_recv_offset[ri-1] + redist_recv_count[ri-1];
        }

        // pack GIO data vectors into particle structs to be distributed by alltoallv
        // ("particle_pos" and "particle_vel" structs defined in util.h) 
        vector<particle_pos> send_particles_pos;
        vector<particle_vel> send_particles_vel;
        vector<particle_pos> recv_particles_pos;
        vector<particle_vel> recv_particles_vel;
     
        for(int n = 0; n < Np; ++n){
            
            particle_pos nextParticle_pos = {r.x[n], r.y[n], r.z[n], r.a[n], r.d[n], r.theta[n], 
                                             r.phi[n], r.id[n], even_redistribute[n]};
            send_particles_pos.push_back(nextParticle_pos);
            
            if(!positionOnly){
                particle_vel nextParticle_vel = {r.vx[n], r.vy[n], r.vz[n], 
                                                 r.rotation[n], r.replication[n], 
                                                 even_redistribute[n]};
                send_particles_vel.push_back(nextParticle_vel);
            }
        }

        recv_particles_pos.resize(redist_recv_offset.back() + redist_recv_count.back());
        if(!positionOnly)
            recv_particles_vel.resize(redist_recv_offset.back() + redist_recv_count.back());

        // now we need to sort our particle data by it's destination rank. As an example;
        // if Np = 12 and numranks = 4, then the above call to comp_rank_scatter will result in
        // even_redistribute = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3}, 
        // which indicates the receiving rank of each particle at position i. In the most recent
        // loop above, we filled the particle struct "myrank" field with the contents of 
        // even_redistribute. So, we can sort the particle objects by that field, in order for our
        // send+offset pair to give the expected result 
        sort(send_particles_pos.begin(), send_particles_pos.end(), comp_rank<particle_pos>);
        if(!positionOnly)
            sort(send_particles_vel.begin(), send_particles_vel.end(), comp_rank<particle_vel>);

        // OK, all read, now to redsitribute the particles evely-ish across ranks
        MPI_Alltoallv(&send_particles_pos[0], &redist_send_count[0], &redist_send_offset[0], particles_mpi_pos,
                      &recv_particles_pos[0], &redist_recv_count[0], &redist_recv_offset[0], particles_mpi_pos, 
                      MPI_COMM_WORLD);
        if(!positionOnly)
            MPI_Alltoallv(&send_particles_vel[0], &redist_send_count[0], &redist_send_offset[0], particles_mpi_vel,
                          &recv_particles_vel[0], &redist_recv_count[0], &redist_recv_offset[0], particles_mpi_vel, 
                          MPI_COMM_WORLD);
        
        // particles now redistributed; find new Np to verify all particles accounted for
        send_particles_pos.clear();
        send_particles_vel.clear();
        Np = redist_recv_offset.back() + redist_recv_count.back(); 
         
        vector<size_t> Np_recv_per_rank(numranks); 
        MPI_Allgather(&Np, 1, MPI_INT64_T, &Np_recv_per_rank[0], 1, MPI_INT64_T, 
                      MPI_COMM_WORLD);

        totalNp = 0;
        for(int ri = 0; ri < numranks; ++ri){
            totalNp += Np_recv_per_rank[ri];
        }
        size_t avg_Np_recv_per_rank = (size_t)(totalNp/numranks);
        if(myrank == 0){
            cout << "Total number of particles after redistribution is " << totalNp << " (about " << 
                    avg_Np_recv_per_rank << " particles per rank)" << endl;
        }   

        MPI_Barrier(MPI_COMM_WORLD);

        // Searching through these particles to perform the cutout later is going to be too 
        // inefficient unless we prepare out particles for a binary search in at least one 
        // dimension. We do this by sorting the recieved particles in ascending order of theta
       
        // arg sort to map new ordering to recv_particles_vel 
        vector<int> theta_argSort(Np);
        int idx=0;
        std::iota(theta_argSort.begin(), theta_argSort.end(), idx++);
        sort(theta_argSort.begin(), theta_argSort.end(), 
             [&](int n, int m){return recv_particles_pos[n].theta < recv_particles_pos[m].theta;} );
         
        // sort recv_particles_pos such that theta is in ascending order
        sort(recv_particles_pos.begin(), recv_particles_pos.end(), comp_by_theta);
        
        MPI_Barrier(MPI_COMM_WORLD);
        stop = MPI_Wtime();
    
        duration = stop - start;
        if(myrank == 0 and timeit == true){ cout << "Redistribution time: " << duration << " s" << endl; }
        redist_times.push_back(duration);
        
        
        ///////////////////////////////////////////////////////////////
        //
        //           Create output files + write buffers
        //
        ///////////////////////////////////////////////////////////////

        // loop over all target halos
       
        // time cutout computation 
        MPI_Barrier(MPI_COMM_WORLD);
        start = MPI_Wtime();
 
        for(int h=0; h<halo_pos.size(); h+=3){
            
            int error = 0; 
            int haloIdx = h/3;
            printHalo = (numHalos < 20) | (haloIdx%20==0) ? 1:0;
            if(myrank == 0 and printHalo){
                cout<< "\n---------- cutout at halo "<< h/3 <<"----------" << endl; 
            }
            
            // instances of buffer struct at file header for output data
            Buffers_write w;

            // open cutout subdirectory for this step...
            ostringstream step_subdir;
            step_subdir << out_dirs[haloIdx] << subdirPrefix << "Cutout" << step_strings[i];
            
            DIR *dir = opendir(step_subdir.str().c_str());
            struct dirent *d;
            char full_path[256];
            char *ext;
            int nf = 0;
            vector<string> files_to_remove;

            // if subdir already exists, make sure it's empty, because overwriting
            // binary files isn't always clean. If 'overwrite' is true, the delete 
            // all binary files in the subdir and continue. 
            // Only have rank 0 do this.
            if(myrank == 0){
                
                // doesn't exist; create the subdir
                if(dir == NULL){
                    mkdir(step_subdir.str().c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IXOTH);
                    if(printHalo){
                        cout << "Created subdir: " << step_subdir.str() << endl;
                    }
                }
                
                // does exist; handle
                else{
                    while((d = readdir(dir)) != NULL){ if(++nf>2){ break;} }
                    
                    // empty but exists; do nothing
                    if(nf<=2){
                        closedir(dir);
                        if(printHalo){
                            cout << "Entered subdir: " << step_subdir.str() << endl;
                        }
                    }

                    // not empty; abort
                    else if(nf > 2 and overwrite == false){
                        closedir(dir);
                        if(printHalo){
                            cout << "\nDirectory " << step_subdir.str() << " is non-empty; skipping halo" << endl;
                        }
                        error = 2;
                    }

                    // not empty; overwrite
                    else if(nf > 2 and overwrite == true){
                        if(printHalo){
                            cout << "\nDirectory " << step_subdir.str() << " is non-empty";
                            cout << " and --overwrite flag was passed; removing binary files here" << endl;
                        }
                        while( (d = readdir(dir)) != NULL){
                            
                            // just to be safe
                            if (0==strcmp(d->d_name, ".") || 
                                0==strcmp(d->d_name, "..")) { continue; }

                            sprintf(full_path, "%s%s%s", step_subdir.str().c_str(), "/", d->d_name);
                            ext = strrchr(d->d_name, '.');

                            // make sure directory only contains binary files
                            if( strcmp(ext, ".bin") != 0){
                                if(printHalo){
                                    cout << "\nDirectory contains non-cutout files! Maybe passed wrong path?" << endl;
                                }
                                error = 1;
                                break;
                            }
                            files_to_remove.push_back(full_path);
                        }
                        
                        if(error == 0){
                            for(int l = 0; l < files_to_remove.size(); ++l){
                                if(verbose == true and (numHalos < 20) | (haloIdx%20==0) ? 1:0){
                                    cout << "\nDELETING " << files_to_remove[l] << endl;
                                }
                                remove(files_to_remove[l].c_str());
                            }
                        }
                        closedir(dir);
                    }
                }    
            }

            // check for potential errors raised above
            // error = 1 is fatal and exits. error = 2 just skips the current halo
            MPI_Bcast(&error, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if(error == 1){ 
                MPI_Finalize();
                exit(EXIT_FAILURE);
            }
            else if(error == 2){ continue; }


            // create binary files for cutout output
            MPI_File id_file;
            MPI_File redshift_file;
            MPI_File x_file;
            MPI_File y_file;
            MPI_File z_file;
            MPI_File vx_file;
            MPI_File vy_file;
            MPI_File vz_file;
            MPI_File rotation_file;
            MPI_File replication_file;
            MPI_File theta_file;
            MPI_File phi_file;

            MPI_Request id_req;
            MPI_Request redshift_req;
            MPI_Request x_req;
            MPI_Request y_req;
            MPI_Request z_req;
            MPI_Request vx_req;
            MPI_Request vy_req;
            MPI_Request vz_req;
            MPI_Request rotation_req;
            MPI_Request replication_req;
            MPI_Request theta_req;
            MPI_Request phi_req;
            
            ostringstream id_file_name;
            ostringstream redshift_file_name;
            ostringstream x_file_name;
            ostringstream y_file_name;
            ostringstream z_file_name;
            ostringstream vx_file_name;
            ostringstream vy_file_name;
            ostringstream vz_file_name;
            ostringstream rotation_file_name;
            ostringstream replication_file_name;
            ostringstream theta_file_name;
            ostringstream phi_file_name; 

            id_file_name << step_subdir.str() << "/id." << step << ".bin";
            redshift_file_name << step_subdir.str() << "/redshift." << step << ".bin";
            x_file_name << step_subdir.str() << "/x."<< step <<".bin";
            y_file_name << step_subdir.str() << "/y."<< step <<".bin";
            z_file_name << step_subdir.str() << "/z."<< step <<".bin";
            vx_file_name << step_subdir.str() << "/vx."<< step <<".bin";
            vy_file_name << step_subdir.str() << "/vy."<< step <<".bin";
            vz_file_name << step_subdir.str() << "/vz."<< step <<".bin";
            rotation_file_name << step_subdir.str() << "/rotation."<< step <<".bin";
            replication_file_name << step_subdir.str() << "/replication."<< step <<".bin";
            theta_file_name << step_subdir.str() << "/theta." << step << ".bin";
            phi_file_name << step_subdir.str() << "/phi." << step << ".bin";

            if(myrank == 0 and printHalo){
                cout<<"opening output files"<<endl;
            }

            MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(id_file_name.str().c_str()), 
                    MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &id_file);
            MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(x_file_name.str().c_str()), 
                    MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &x_file);
            MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(y_file_name.str().c_str()),
                    MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &y_file);
            MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(z_file_name.str().c_str()),
                    MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &z_file);
            MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(redshift_file_name.str().c_str()),
                    MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &redshift_file);
            MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(theta_file_name.str().c_str()),
                    MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &theta_file);
            MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(phi_file_name.str().c_str()),
                    MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &phi_file);
            if(!positionOnly){
                MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(vx_file_name.str().c_str()),
                        MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &vx_file);
                MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(vy_file_name.str().c_str()),
                        MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &vy_file);
                MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(vz_file_name.str().c_str()),
                        MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &vz_file);
                MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(rotation_file_name.str().c_str()),
                        MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &rotation_file);
                MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(replication_file_name.str().c_str()),
                        MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &replication_file);
            }
            
            ///////////////////////////////////////////////////////////////
            //
            //                         Do cutting
            //
            ///////////////////////////////////////////////////////////////
        
            // let's also time the computation per-rank
            clock_t thisRank_start = clock();

            if(myrank == 0 and printHalo){
                cout << "converting positions..." << endl;
            }
            
            int cutout_size = 0;
            
            // define vectors joining fov corners to particle (see comments/diagram above)
            vector<float> AM(2);
            vector<float> BM(2);
        
            // all of this ranks recieved particles were sorted, after read-in, by their 
            // "theta" attribute. So, we can do a binary search for our rough theta bounds
            // to limit our search to an annulus around the sky parallel to the equator...
            particle_pos left_dummy;
            left_dummy.theta = theta_cut_rough[haloIdx][0];
            particle_pos right_dummy;
            right_dummy.theta = theta_cut_rough[haloIdx][1];
            
            auto leftCut_iter = std::lower_bound(recv_particles_pos.begin(), recv_particles_pos.end(), 
                                                 left_dummy, comp_by_theta);
            auto rightCut_iter = std::upper_bound(recv_particles_pos.begin(), recv_particles_pos.end(), 
                                                  right_dummy, comp_by_theta);
            
            int minN = std::distance(recv_particles_pos.begin(), leftCut_iter);
            int maxN = std::distance(recv_particles_pos.begin(), rightCut_iter);
                            
            auto posIdx_iter = std::find(theta_argSort.begin(), theta_argSort.end(), 0);
            int posIdx = std::distance(theta_argSort.begin(), posIdx_iter);
            
            // Now, brute force search on phi to finish rough cut out
            for (int n=minN; n<maxN; ++n) {
                
                float d = recv_particles_pos[n].d; 
                float theta = recv_particles_pos[n].theta;
                float phi = recv_particles_pos[n].phi;

                if (phi > phi_cut_rough[haloIdx][0] && phi < phi_cut_rough[haloIdx][1]) {
                 
                    // of the particles surviving the rough cut, let's do a proper rotation 
                    // on them to find the true cutout memership, and return cluster-centric 
                    // angular coordinates
                    
                    // do coordinate rotation center halo at (r, 90, 0)
                    // B and k are the angle and axis of rotation, respectively,
                    // calculated near the beginning of this function
                    float tmp_v[] = {recv_particles_pos[n].x, recv_particles_pos[n].y, recv_particles_pos[n].z};
                    vector<float> v(tmp_v, tmp_v+3);
                    vector<float> v_rot;
                    v_rot = matVecMul(R[haloIdx], v);

                    // spherical coordinate transformation
                    d = (float)sqrt(v_rot[0]*v_rot[0] + v_rot[1]*v_rot[1] + 
                                    v_rot[2]*v_rot[2]);
                    float v_theta = acos(v_rot[2]/d) * 180.0 / PI * ARCSEC;
                    float v_phi;

                    // prevent NaNs on y-z plane
                    if(v_rot[0] == 0 && v_rot[1] > 0)
                        v_phi = 90.0 * ARCSEC;
                    else if(v_rot[0] == 0 && v_rot[1] < 0)
                        v_phi = -90.0 * ARCSEC;
                    else
                        v_phi = atan(v_rot[1]/v_rot[0]) * 180.0 / PI * ARCSEC; 
                 
                    // do final cut
                    if (v_theta > theta_cut[haloIdx][0] && v_theta < theta_cut[haloIdx][1] && 
                        v_phi > phi_cut[haloIdx][0] && v_phi < phi_cut[haloIdx][1] ) {

                        // get redshift from scale factor
                        float zz = aToZ(recv_particles_pos[n].a);  
                        
                        // spherical corrdinate transform of rotated positions
                        w.theta.push_back(v_theta);
                        w.phi.push_back(v_phi);

                        // other columns
                        w.x.push_back(recv_particles_pos[n].x);
                        w.y.push_back(recv_particles_pos[n].y);
                        w.z.push_back(recv_particles_pos[n].z);
                        w.redshift.push_back(zz);
                        w.id.push_back(recv_particles_pos[n].id);
                        if(!positionOnly){
                            
                            // get index of matching entry in recv_particles_vel
                            auto velIdx_iter = std::find(theta_argSort.begin(), theta_argSort.end(), n);
                            int velIdx = std::distance(theta_argSort.begin(), velIdx_iter);
                            
                            w.vx.push_back(recv_particles_vel[theta_argSort[velIdx]].vx);
                            w.vy.push_back(recv_particles_vel[theta_argSort[velIdx]].vy);
                            w.vz.push_back(recv_particles_vel[theta_argSort[velIdx]].vz);
                            w.rotation.push_back(recv_particles_vel[theta_argSort[velIdx]].rotation);
                            w.replication.push_back(recv_particles_vel[theta_argSort[velIdx]].replication);
                        }
                        cutout_size++;

                        /*
                        // DEBUG
                        // print out individual particle info

                        if(myrank == 1){
                            cout << endl << "Particle " << recv_particles_pos[n].id << ":   " << endl << 
                            "x: " << recv_particles_pos[n].x << endl << 
                            "y: " << recv_particles_pos[n].y << endl << 
                            "z: " << recv_particles_pos[n].z << endl <<
                            "a: " << recv_particles_pos[n].a << endl <<
                            "rs: " << zz << endl <<
                            "theta: " << v_theta << endl << 
                            "phi: " << v_phi << endl; 
                        }
                        */
                    }
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            
            stop = MPI_Wtime();
            duration = stop - start;
            if(myrank == 0 and timeit==true and printHalo){
                cout << "cutout computation time: " << duration << " s" << endl; 
            }
            cutout_times.push_back(duration);
            
            if(verbose == true and timeit == true and printHalo){
                
                // check load balancing (all ranks should have taken more or less the same amount of time here)
                clock_t thisRank_end = clock();
                clock_t thisRank_time = thisRank_end - thisRank_start;
                double thisRank_secs = thisRank_time / (double) CLOCKS_PER_SEC;
                vector<double> allRank_secs(numranks);
                
                MPI_Allgather(&thisRank_secs, 1, MPI_DOUBLE, 
                              &allRank_secs[0], 1, MPI_DOUBLE, MPI_COMM_WORLD);
     
                if(myrank == 0){
                    double min_compTime = 9999;
                    double max_compTime = 0;
                    double mean_compTime;
                    double sq_compTime;
                    double std_compTime;

                    cout << "allRank_secs: [";
                    for(int cc = 0; cc < numranks; ++cc){
                        cout << allRank_secs[cc] << ", ";
                    }
                    cout << endl;        

                    for(int cc = 0; cc < numranks; ++cc){
                        if(allRank_secs[cc] < min_compTime){ min_compTime = allRank_secs[cc];}
                        if(allRank_secs[cc] > max_compTime){ max_compTime = allRank_secs[cc];}
                    }

                    mean_compTime = accumulate(allRank_secs.begin(), allRank_secs.end(), 0.0) / allRank_secs.size();
                    sq_compTime = inner_product(allRank_secs.begin(), allRank_secs.end(), allRank_secs.begin(), 0.0);
                    std_compTime = sqrt(sq_compTime / allRank_secs.size() - mean_compTime * mean_compTime);

                    cout << "Min rank computation time: " << min_compTime << " s" << endl; 
                    cout << "Max rank computation time: " << max_compTime << " s" << endl; 
                    cout << "Mean rank computation time: " << mean_compTime << " s" << endl; 
                    cout << "Std dev rank computation time: " << std_compTime << " s" << endl; 
                }
            }

            ///////////////////////////////////////////////////////////////
            //
            //                   write out
            //
            ///////////////////////////////////////////////////////////////
       
            // time write out 
            MPI_Barrier(MPI_COMM_WORLD);
            start = MPI_Wtime();

            // define MPI file writing offset for the current rank --
            // This offset will be the sum of elements in all lesser ranks,
            // multiplied by the type size for each file    
            w.np_count.clear();
            w.np_count.resize(numranks);
            w.np_offset.clear();
            w.np_offset.push_back(0);

            // get number of elements in each ranks portion of cutout 
            MPI_Allgather(&cutout_size, 1, MPI_INT, 
                          &w.np_count[0], 1, MPI_INT, MPI_COMM_WORLD);
            
            // compute each ranks writing offset
            for(int j=1; j < numranks; ++j){
                w.np_offset.push_back(w.np_offset[j-1] + w.np_count[j-1]);
            }
            MPI_Barrier(MPI_COMM_WORLD); 
           
            // print out offset vector for verification
            if(myrank == 0 and printHalo){
                if(numranks < 20){
                    cout << "rank object counts: [";
                    for(int m=0; m < numranks; ++m){ cout << w.np_count[m] << ","; }
                    cout << "]" << endl;
                    cout << "rank offsets: [";
                    for(int m=0; m < numranks; ++m){ cout << w.np_offset[m] << ","; }
                    cout << "]" << endl;
                } else {
                   int numEmpty = count(&w.np_count[0], &w.np_count[numranks], 0);
                   cout << numranks - numEmpty << " of " << numranks << 
                   " ranks found members within cutout field of view" << endl;
                }
                cout << "Writing files..." << endl;
            }

            MPI_Offset offset_posvel = sizeof(POSVEL_T) * w.np_offset[myrank];
            MPI_Offset offset_id = sizeof(ID_T) * w.np_offset[myrank];
            MPI_Offset offset_float = sizeof(float) * w.np_offset[myrank];
            MPI_Offset offset_int = sizeof(int) * w.np_offset[myrank];
            MPI_Offset offset_int32 = sizeof(int32_t) * w.np_offset[myrank];

            // write
            MPI_File_seek(id_file, offset_id, MPI_SEEK_SET);
            MPI_File_iwrite(id_file, &w.id[0], w.id.size(), MPI_INT64_T, &id_req);
            MPI_Wait(&id_req, MPI_STATUS_IGNORE);

            MPI_File_seek(x_file, offset_posvel, MPI_SEEK_SET);
            MPI_File_iwrite(x_file, &w.x[0], w.x.size(), MPI_FLOAT, &x_req);
            MPI_Wait(&x_req, MPI_STATUS_IGNORE);

            MPI_File_seek(y_file, offset_posvel, MPI_SEEK_SET);
            MPI_File_iwrite(y_file, &w.y[0], w.y.size(), MPI_FLOAT, &y_req);
            MPI_Wait(&y_req, MPI_STATUS_IGNORE);
            
            MPI_File_seek(z_file, offset_posvel, MPI_SEEK_SET);
            MPI_File_iwrite(z_file, &w.z[0], w.z.size(), MPI_FLOAT, &z_req);
            MPI_Wait(&z_req, MPI_STATUS_IGNORE);
            
            MPI_File_seek(theta_file, offset_posvel, MPI_SEEK_SET);
            MPI_File_iwrite(theta_file, &w.theta[0], w.theta.size(), MPI_FLOAT, &theta_req);
            MPI_Wait(&theta_req, MPI_STATUS_IGNORE);
            
            MPI_File_seek(phi_file, offset_posvel, MPI_SEEK_SET);
            MPI_File_iwrite(phi_file, &w.phi[0], w.phi.size(), MPI_FLOAT, &phi_req);
            MPI_Wait(&phi_req, MPI_STATUS_IGNORE);
            
            MPI_File_seek(redshift_file, offset_posvel, MPI_SEEK_SET);
            MPI_File_iwrite(redshift_file, &w.redshift[0], w.redshift.size(), MPI_FLOAT, &redshift_req);
            MPI_Wait(&redshift_req, MPI_STATUS_IGNORE);
            
            if(!positionOnly){
                MPI_File_seek(vx_file, offset_posvel, MPI_SEEK_SET);
                MPI_File_iwrite(vx_file, &w.vx[0], w.vx.size(), MPI_FLOAT, &vx_req);
                MPI_Wait(&vx_req, MPI_STATUS_IGNORE);

                MPI_File_seek(vy_file, offset_posvel, MPI_SEEK_SET);
                MPI_File_iwrite(vy_file, &w.vy[0], w.vy.size(), MPI_FLOAT, &vy_req);
                MPI_Wait(&vy_req, MPI_STATUS_IGNORE);
                
                MPI_File_seek(vz_file, offset_posvel, MPI_SEEK_SET);
                MPI_File_iwrite(vz_file, &w.vz[0], w.vz.size(), MPI_FLOAT, &vz_req);
                MPI_Wait(&vz_req, MPI_STATUS_IGNORE);
                
                MPI_File_seek(rotation_file, offset_posvel, MPI_SEEK_SET);
                MPI_File_iwrite(rotation_file, &w.rotation[0], w.rotation.size(), 
                                MPI_FLOAT, &rotation_req);
                MPI_Wait(&rotation_req, MPI_STATUS_IGNORE);
                
                MPI_File_seek(replication_file, offset_posvel, MPI_SEEK_SET);
                MPI_File_iwrite(replication_file, &w.replication[0], w.replication.size(), 
                                MPI_FLOAT, &replication_req);
                MPI_Wait(&replication_req, MPI_STATUS_IGNORE);
            }

            MPI_File_close(&id_file);
            MPI_File_close(&x_file);
            MPI_File_close(&y_file);
            MPI_File_close(&z_file);
            MPI_File_close(&redshift_file);
            MPI_File_close(&theta_file);
            MPI_File_close(&phi_file);
            if(!positionOnly){
                MPI_File_close(&vx_file);
                MPI_File_close(&vy_file);
                MPI_File_close(&vz_file);
                MPI_File_close(&rotation_file);
                MPI_File_close(&replication_file);
            }
        
            MPI_Barrier(MPI_COMM_WORLD);
            stop = MPI_Wtime();
        
            duration = stop - start;
            if(myrank == 0 and timeit == true and printHalo){ 
                cout << "write time: " << duration << " s" << endl; 
            }
            write_times.push_back(duration);
        }

    }
    
    if(myrank == 0 and timeit == true){
        
        cout << "\nread_times = np.array([";
        for(int hh = 0; hh < read_times.size(); ++hh){
            cout << read_times[hh];
            if(hh < read_times.size()-1){ cout << ", "; }
        }
        cout << "]" << endl;
        
        cout << "redist_times = np.array([";
        for(int hh = 0; hh < redist_times.size(); ++hh){
            cout << redist_times[hh];
            if(hh < redist_times.size()-1){ cout << ", "; }
        }
        cout << "]" << endl;
        
        cout << "cutout_times = np.array([";
        for(int hh = 0; hh < cutout_times.size(); ++hh){
            cout << cutout_times[hh];
            if(hh < cutout_times.size()-1){ cout << ", "; }
        }
        cout << "]" << endl;
        
        cout << "write_times = np.array([";
        for(int hh = 0; hh < write_times.size(); ++hh){
            cout << write_times[hh];
            if(hh < write_times.size()-1){ cout << ", "; }
        }
        cout << "]" << endl;
    }
}
