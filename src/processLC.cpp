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
    string subdirPrefix;
    int prefixSize;
    if(myrank == 0){ 
        cout << "\nReading directory: " << dir_name << endl;
        
        vector<string> subdirs;
        getLCSubdirs(dir_name, subdirs);
        
        cout << "Found subdirs:" << endl;
        for (vector<string>::const_iterator i = subdirs.begin(); i != subdirs.end(); ++i){
            cout << *i << ' ';
        }
        cout << endl;

        // find the prefix (chars before the step number) in the subdirectory names.
        // It is assumed that all subdirs have the same prefix.
        for(string::size_type j = 0; j < subdirs[0].size(); ++j){
            if( isdigit(subdirs[0][j]) > 0){
                subdirPrefix = subdirs[0].substr(0, j);
                break;
            }
        }
        cout << "Subdir prefix is: " << subdirPrefix << endl;
        prefixSize = subdirPrefix.size();
    }
    MPI_Bcast(&prefixSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(myrank != 0){ subdirPrefix.resize(prefixSize); }
    MPI_Bcast(const_cast<char*>(subdirPrefix.data()), prefixSize, MPI_CHAR, 0, MPI_COMM_WORLD);

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
        int fname_size;
        ostringstream file_name_stream;
        file_name_stream << dir_name << subdirPrefix << step_strings[i]; 
        if(myrank == 0){
            getLCFile(file_name_stream.str(), file_name);
            fname_size = file_name.size();
        } 
        
        MPI_Bcast(&fname_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if(myrank != 0){ file_name.resize(fname_size); }
        MPI_Bcast(const_cast<char*>(file_name.data()), fname_size, MPI_CHAR, 0, MPI_COMM_WORLD);
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
               int numranks, bool verbose, bool timeit, bool overwrite, bool positionOnly, 
               bool forceWriteProps, bool propsOnly){


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
        printHalo = (numHalos < 20) | (haloIdx%100==0) ? 1:0;
        
        // get next three values in halo_pos
        float tmp_pos[] = {halo_pos[h], halo_pos[h+1], halo_pos[h+2]};
        vector<float> this_halo_pos(tmp_pos, tmp_pos+3);

        // get next three or four values in halo_props 
        // (redshift, step, sod_mass, sod_radius, sod_concentration, sod_concentration_error) or 
        // (redshift, step, fof_mass) 
        // last 3 elements of tmp_props will be garbage if numProps == 3 (massDef == 'fof'),
        // but that's okay. 
        int numProps = halo_props.size() / numHalos;
        float tmp_props[] = {halo_props[numProps*haloIdx], halo_props[numProps*haloIdx+1], 
                             halo_props[numProps*haloIdx+2], halo_props[numProps*haloIdx+3],
                             halo_props[numProps*haloIdx+4], halo_props[numProps*haloIdx+5]};
        if(numProps != 3 and numProps != 6){
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
        // cutout specifications, and run meta data. By default, skip this step if that property
        // file already exists. 
        // If the code had been updated to change the content of that file, 
        // and you want old outputs to be overwritten, then just change forceWrite to true.
        if(myrank == 0){ 
            ofstream props_file;
            ostringstream props_file_name; 
            props_file_name << out_dirs[haloIdx]<< "/properties.csv";

            if(does_file_exist(props_file_name.str()) == false or forceWriteProps == true){
                props_file.open(props_file_name.str().c_str());
                
                props_file << "#halo_redshift" << ", " << "halo_lc_shell" << ", ";
                if(numProps == 6)
                    props_file << "sod_halo_mass" << ", " <<  "sod_halo_radius" << ", " 
                               << "sod_halo_cdelta" << ", " << "sod_halo_cdelta_error";
                else
                    props_file << "fof_halo_mass";
                props_file << ", " << "halo_lc_x" << ", " << "halo_lc_y" << ", " << "halo_lc_z" << ", "
                           << "boxRadius_Mpc" << ", " << "boxRadius_arcsec" << "\n";

                for(int i=0; i<this_halo_props.size(); ++i)
                    props_file << this_halo_props[i] << ", ";
                for(int i=0; i<this_halo_pos.size(); ++i)
                    props_file << this_halo_pos[i] << ", ";
                props_file << atan(halfBoxLength) * halo_r << ", " << halfBoxLength * 180.0/PI * ARCSEC << "\n";
                 
                props_file.close();
                if(printHalo){ cout << "wrote halo info to properties.csv" << endl; }
            
            }else{ 
                if(printHalo){ cout << "properties.csv already exists; skipping write" << endl; }
            }
        }
    }

    if(propsOnly == true){
        return;
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
    vector<double> sort_times;
    vector<double> cutout_times; 
    vector<double> write_times; 
    double start;
    double stop;
    double duration;
    
    for (int i=0; i<step_strings.size(); ++i){
   
        // time read in 
        MPI_Barrier(MPI_COMM_WORLD);
        start = MPI_Wtime();
         
        // continue if this is the step at z=0 (lightcone volume zero)
        step =atoi(step_strings[i].c_str());
        if(step == 499){ continue;}

        // find header file
        string file_name;
        int fname_size;
        ostringstream file_name_stream;
        file_name_stream << dir_name << subdirPrefix << step_strings[i]; 
        
        if(myrank == 0){
            cout << "\n=================================================" << endl;
            cout << "============== Working on step "<< step_strings[i] <<" ==============\n" << endl; 
           
            if(myrank == 0){ cout << "getting lc file..." << endl; } 
            
            getLCFile(file_name_stream.str(), file_name);
            fname_size = file_name.size();
        } 
            
        MPI_Barrier(MPI_COMM_WORLD); 
        if(myrank == 0){ cout << "brodcasting file name..." << endl; } 
        MPI_Barrier(MPI_COMM_WORLD); 
        
        MPI_Bcast(&fname_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if(myrank != 0){ file_name.resize(fname_size); }
        MPI_Bcast(const_cast<char*>(file_name.data()), fname_size, MPI_CHAR, 0, MPI_COMM_WORLD);
        file_name_stream << "/" << file_name;
        
       
        ///////////////////////////////////////////////////////////////
        //
        //         do reading + spherical coordinate transform
        //
        ///////////////////////////////////////////////////////////////
       
        // instances of buffer struct at file header for particle data after redisribution
        Buffers_read recv_particles;
        size_t Np = 0;
        
        { // start read-in scope; read buffers r destroyed after here

            // instances of buffer struct at file header for read in data
            Buffers_read r;
            
            // setup gio
            MPI_Barrier(MPI_COMM_WORLD); 
            if(myrank == 0){ cout << "setting up gio..." << endl; } 
            MPI_Barrier(MPI_COMM_WORLD); 
            
            unsigned Method = GenericIO::FileIOPOSIX;
            const char *EnvStr = getenv("GENERICIO_USE_MPIIO");
            if(EnvStr && string(EnvStr) == "1"){
                Method = GenericIO::FileIOMPI;  
            }
            
            MPI_Barrier(MPI_COMM_WORLD); 
            if(myrank == 0){ cout << "done setting up gio..." << endl; } 
            MPI_Barrier(MPI_COMM_WORLD); 

            // create gio reader, open lightcone file header in new scope
            {
                MPI_Barrier(MPI_COMM_WORLD); 
                if(myrank == 0){ cout << "Opening file: " << file_name_stream.str() << endl; }
                MPI_Barrier(MPI_COMM_WORLD); 
                
                GenericIO GIO(MPI_COMM_WORLD, file_name_stream.str(), Method);
                GIO.openAndReadHeader(GenericIO::MismatchRedistribute);

                MPI_Barrier(MPI_COMM_WORLD);
                Np = GIO.readNumElems();
               
                // resize buffers
                resize_read_buffers(r, Np, positionOnly, GIO.requestedExtraSpace());

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
            resize_read_buffers(r, Np, positionOnly);
            if(myrank == 0){ cout<<"done resizing"<<endl; }
            
            // calc d, theta, and phi per particle
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
            
            // ************ debug ***************** 
            if(myrank == 0){ cout << endl;}
            MPI_Barrier(MPI_COMM_WORLD);
            if(numranks == 2){
                float max_theta_rank = -1e9;
                float min_theta_rank = 1e9;
                float avg_theta_rank = 0;
                float avg_phi_rank = 0;
                for(int qq=0; qq < Np; qq++){
                    if(r.theta[qq] > max_theta_rank)
                        max_theta_rank = r.theta[qq];
                    if(r.theta[qq] < min_theta_rank)
                        min_theta_rank = r.theta[qq];
                    avg_theta_rank = avg_theta_rank + r.theta[qq];
                    avg_phi_rank = avg_phi_rank + r.phi[qq];
                }
                avg_theta_rank = avg_theta_rank / r.theta.size();
                avg_phi_rank = avg_phi_rank / r.theta.size();
                
                if(myrank==0){    
                    cout << "particles at rank " << myrank << ": " << Np << endl;
                    cout << "[min, avg, max] theta at rank " << myrank << ": [" << 
                            min_theta_rank << ", " << avg_theta_rank << ", " << max_theta_rank << 
                            "]" << endl;
                    cout << "avg phi at rank " << myrank << ": " << avg_phi_rank << endl << endl;
                }
                MPI_Barrier(MPI_COMM_WORLD);
                if(myrank==1){    
                    cout << "particles at rank " << myrank << ": " << Np << endl;
                    cout << "[min, avg, max] theta at rank " << myrank << ": [" << 
                            min_theta_rank << ", " << avg_theta_rank << ", " << max_theta_rank << 
                            "]" << endl;
                    cout << "avg phi at rank " << myrank << ": " << avg_phi_rank << endl << endl;
                }
                
            }
            // ************ debug *****************  


            ///////////////////////////////////////////////////////////////
            //
            //           evenly distribute particles to all ranks
            //
            ///////////////////////////////////////////////////////////////
            
            // now, we have likely ended up in the situation where most of the data read
            // resides in only a small subset of the ranks in this communicator. This is becuase
            // (I think) only a small fraction of the number of ranks that generated the lightcone 
            // at any particular step intersected that lightcone shell. That information is encoded 
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
            if(myrank == 0){cout << "0" << endl;}
            comp_rank_scatter(Np, even_redistribute, numranks);
            for(int ri = 0; ri < numranks; ++ri){
                redist_send_count[ri] = count(&even_redistribute[0], &even_redistribute[Np], ri);
            }
            
            // get number of particles to recieve from every other rank
            if(myrank == 0){cout << "1" << endl;}
            MPI_Alltoall(&redist_send_count[0], 1, MPI_INT, &redist_recv_count[0], 1, MPI_INT, MPI_COMM_WORLD);

            // compute sending+recieving offsets to/from each other rank
            if(myrank == 0){cout << "2" << endl;}
            for(int ri=1; ri < numranks; ++ri){
                redist_send_offset[ri] = redist_send_offset[ri-1] + redist_send_count[ri-1];
                redist_recv_offset[ri] = redist_recv_offset[ri-1] + redist_recv_count[ri-1];
            }
            
            // resize reciving buffers
            if(myrank == 0){cout << "3" << endl;}
            int tot_num_recv = redist_recv_offset.back() + redist_recv_count.back();
            resize_read_buffers(recv_particles, tot_num_recv, positionOnly);
     
            // OK, all read, now to redsitribute the particles evely-ish across ranks
            if(myrank == 0){cout << "4" << endl;}
            POSVEL_T *fcols_send_pos[] = {&r.x[0], &r.y[0], &r.z[0], &r.d[0], &r.theta[0], &r.phi[0], &r.a[0]};
            POSVEL_T *fcols_recv_pos[] = {&recv_particles.x[0], &recv_particles.y[0], &recv_particles.z[0], 
                                          &recv_particles.d[0], &recv_particles.theta[0], &recv_particles.phi[0], 
                                          &recv_particles.a[0]};

            if(myrank == 0){cout << "5" << endl;}
            POSVEL_T *fcols_send_vel[] = {&r.vx[0], &r.vy[0], &r.vz[0]};
            POSVEL_T *fcols_recv_vel[] = {&recv_particles.vx[0], &recv_particles.vy[0], &recv_particles.vz[0]};
            
            if(myrank == 0){cout << "6" << endl;}
            for(int coln; coln < 7; coln++){
                MPI_Alltoallv(fcols_send_pos[coln], &redist_send_count[0], &redist_send_offset[0], MPI_FLOAT,
                              fcols_recv_pos[coln], &redist_recv_count[0], &redist_recv_offset[0], MPI_FLOAT, 
                              MPI_COMM_WORLD);
            }
            MPI_Alltoallv(&r.id[0], &redist_send_count[0], &redist_send_offset[0], MPI_INT64_T,
                          &recv_particles.id[0], &redist_recv_count[0], &redist_recv_offset[0], MPI_INT64_T, 
                          MPI_COMM_WORLD);

            if(myrank == 0){cout << "7" << endl;}
            if(!positionOnly){
                for(int coln; coln < 3; coln++){
                    MPI_Alltoallv(fcols_send_vel[coln], &redist_send_count[0], &redist_send_offset[0], MPI_FLOAT,
                                  fcols_recv_vel[coln], &redist_recv_count[0], &redist_recv_offset[0], MPI_FLOAT, 
                                  MPI_COMM_WORLD);
                }
                MPI_Alltoallv(&r.replication[0], &redist_send_count[0], &redist_send_offset[0], MPI_INT32_T,
                              &recv_particles.replication[0], &redist_recv_count[0], &redist_recv_offset[0], 
                              MPI_INT32_T, MPI_COMM_WORLD);
                MPI_Alltoallv(&r.rotation[0], &redist_send_count[0], &redist_send_offset[0], MPI_INT,
                              &recv_particles.rotation[0], &redist_recv_count[0], &redist_recv_offset[0], 
                              MPI_INT, MPI_COMM_WORLD);
            }
            
            if(myrank == 0){cout << "8" << endl;}
            // particles now redistributed; find new Np to verify all particles accounted for
            Np = recv_particles.id.size(); 
             
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
            stop = MPI_Wtime(); 
            duration = stop - start;
            if(myrank == 0 and timeit == true){ cout << "Redistribution time: " << duration << " s" << endl; }
            redist_times.push_back(duration);
            
            // ************ debug *****************  
            if(myrank == 0){ cout << endl;}
            MPI_Barrier(MPI_COMM_WORLD);
            if(numranks == 2){
                float max_theta_rank = -1e9;
                float min_theta_rank = 1e9;
                float avg_theta_rank = 0;
                float avg_phi_rank = 0;
                for(int qq=0; qq < Np; qq++){
                    if(recv_particles.theta[qq] > max_theta_rank)
                        max_theta_rank = recv_particles.theta[qq];
                    if(recv_particles.theta[qq] < min_theta_rank)
                        min_theta_rank = recv_particles.theta[qq];
                    avg_theta_rank = avg_theta_rank + recv_particles.theta[qq];
                    avg_phi_rank = avg_phi_rank + recv_particles.phi[qq];
                }
                avg_theta_rank = avg_theta_rank / recv_particles.theta.size();
                avg_phi_rank = avg_phi_rank / recv_particles.theta.size();
                
                if(myrank==0){    
                    cout << "particles at rank " << myrank << ": " << Np << endl;
                    cout << "[min, avg, max] theta at rank " << myrank << ": [" << 
                            min_theta_rank << ", " << avg_theta_rank << ", " << max_theta_rank << 
                            "]" << endl;
                    cout << "avg phi at rank " << myrank << ": " << avg_phi_rank << endl << endl;
                }
                MPI_Barrier(MPI_COMM_WORLD);
                if(myrank==1){    
                    cout << "particles at rank " << myrank << ": " << Np << endl;
                    cout << "[min, avg, max] theta at rank " << myrank << ": [" << 
                            min_theta_rank << ", " << avg_theta_rank << ", " << max_theta_rank << 
                            "]" << endl;
                    cout << "avg phi at rank " << myrank << ": " << avg_phi_rank << endl << endl;
                }
                
            }
            // ************ debug *****************  
         
            MPI_Barrier(MPI_COMM_WORLD);
        
        } // end of read-in scope (read buffers r destroyed here)


        ///////////////////////////////////////////////////////////////
        //
        //           Sort particles by theta for binary search
        //
        ///////////////////////////////////////////////////////////////

        // Searching through these particles to perform the cutout later is going to be too 
        // inefficient unless we prepare out particles for a binary search in at least one 
        // dimension. We do this by sorting the recieved particles in ascending order of theta
        
        // time sort 
        MPI_Barrier(MPI_COMM_WORLD);
        start = MPI_Wtime();
        
        // arg sort to map new ordering to recv_particles_vel 
        vector<int> theta_argSort(Np);
        int idx=0;
        std::iota(theta_argSort.begin(), theta_argSort.end(), idx++);
        stable_sort(theta_argSort.begin(), theta_argSort.end(), 
             [&](int n, int m){return recv_particles.theta[n] < recv_particles.theta[m];} );
        
        // sort recv_particles_pos such that theta is in ascending order
        stable_sort(recv_particles.theta.begin(), recv_particles.theta.end());
        
        MPI_Barrier(MPI_COMM_WORLD);
        stop = MPI_Wtime(); 
        duration = stop - start;
        if(myrank == 0 and timeit == true){ cout << "Particle sort time: " << duration << " s" << endl; }
        sort_times.push_back(duration);
        

        ///////////////////////////////////////////////////////////////
        //
        //                 Loop over all target halos
        //
        ///////////////////////////////////////////////////////////////
 
        MPI_Barrier(MPI_COMM_WORLD);
        for(int h=0; h<halo_pos.size(); h+=3){
            
            int error = 0; 
            int haloIdx = h/3;
            printHalo = (numHalos < 20) | (haloIdx%100==0) ? 1:0;
            if(myrank == 0 and printHalo){
                cout<< "\n---------- cutout at halo "<< h/3 <<"----------" << endl; 
            }
        
            
            ///////////////////////////////////////////////////////////////
            //
            //           Create output files + write buffers
            //
            ///////////////////////////////////////////////////////////////
            
            // instances of buffer struct at file header for output data
            Buffers_write w;

            // open cutout subdirectory for this step...
            // if step subdir already exists, make sure it's empty, because overwriting
            // binary files isn't always clean. If 'overwrite' is true, the delete all 
            // binary files in the subdir and continue. 
            // Only have rank 0 do this.
            ostringstream step_subdir;
            step_subdir << out_dirs[haloIdx] << subdirPrefix << "Cutout" << step_strings[i];
            if(myrank == 0){ 
                error = prepStepSubdir(step_subdir.str(), overwrite, printHalo, verbose);
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
            MPI_File id_file, x_file, y_file, z_file, vx_file, vy_file, vz_file,
                     redshift_file, theta_file, phi_file,
                     rotation_file, replication_file;
            
            MPI_Request id_req, x_req, y_req, z_req, vx_req, vy_req, vz_req,
                        redshift_req, theta_req, phi_req,
                        rotation_req, replication_req;
            
            ostringstream id_file_name, x_file_name, y_file_name, z_file_name, 
                          vx_file_name, vy_file_name, vz_file_name,
                          redshift_file_name, theta_file_name, phi_file_name,
                          rotation_file_name, replication_file_name;

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

            
            ///////////////////////////////////////////////////////////////
            //
            //                         Do cutting
            //
            ///////////////////////////////////////////////////////////////
        
            // time cutout computation 
            MPI_Barrier(MPI_COMM_WORLD);
            start = MPI_Wtime();
        
            // let's also time the computation per-rank
            clock_t thisRank_start = clock();
            clock_t thisRank_end;

            if(myrank == 0 and printHalo){
                cout << "converting positions..." << endl;
            }
            
            int cutout_size = 0;
            int rough_cutout_size = 0;
            
            // define vectors joining fov corners to particle (see comments/diagram above)
            vector<float> AM(2);
            vector<float> BM(2);
        
            // all of this rank's recieved particles were sorted, after read-in, by their 
            // "theta" attribute. So, we can do a binary search for our rough theta bounds
            // to limit our following brute force search to an annulus around the sky 
            // parallel to the equator...
            auto leftCut_iter = std::lower_bound(recv_particles.theta.begin(), recv_particles.theta.end(), 
                                                 theta_cut_rough[haloIdx][0]);
            auto rightCut_iter = std::lower_bound(recv_particles.theta.begin(), recv_particles.theta.end(), 
                                                  theta_cut_rough[haloIdx][1]);
            
            int minN = std::distance(recv_particles.theta.begin(), leftCut_iter);
            int maxN = std::distance(recv_particles.theta.begin(), rightCut_iter);
            
            // run with 2 ranks
            if(numranks == 2){
                if(myrank == 0){
                    bool ordered  = true;
                    for(int i = 0; i < (recv_particles.theta.size()-1); ++i){
                        if(recv_particles.theta[i] > recv_particles.theta[i+1])
                            ordered = false;
                    }
                    cout << endl;
                    cout << "rank 0 ORDERED: " << ordered << endl;
                    cout << "rank 0 PHI ROUGH MIN: " << phi_cut_rough[haloIdx][0] << endl;
                    cout << "rank 0 PHI ROUGH MAX: " << phi_cut_rough[haloIdx][1] << endl;
                    cout << "rank 0 THETA ROUGH MIN: " << theta_cut_rough[haloIdx][0] << endl;
                    cout << "rank 0 THETA ROUGH MAX: " << theta_cut_rough[haloIdx][1] << endl;
                    cout << "rank 0 PHI MIN: " << phi_cut[haloIdx][0] << endl;
                    cout << "rank 0 PHI MAX: " << phi_cut[haloIdx][1] << endl;
                    cout << "rank 0 THETA MIN: " << theta_cut[haloIdx][0] << endl;
                    cout << "rank 0 THETA MAX: " << theta_cut[haloIdx][1] << endl;
                    cout << "rank 0 MAX IDX: " << maxN << endl;
                    cout << "rank 0 MIN IDX: " << minN << endl;
                }
                MPI_Barrier(MPI_COMM_WORLD);
                if(myrank == 1){
                    bool ordered  = true;
                    for(int i = 0; i < (recv_particles.theta.size()-1); ++i){
                        if(recv_particles.theta[i] > recv_particles.theta[i+1])
                            ordered = false;
                    }
                    cout << endl;
                    cout << "rank 1  ORDERED: " << ordered << endl;
                    cout << "rank 1 MAX IDX: " << maxN << endl;
                    cout << "rank 1 MIN IDX: " << minN << endl;
                    cout << endl;
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }
                            
            // Now, brute force search on phi to finish rough cutout
            for (int idx=minN; idx<maxN; ++idx){ 
                
                if(myrank == 0){
                    cout << "RANK " << myrank << " AT IDX: " << idx << endl; 
                }
                MPI_Barrier(MPI_COMM_WORLD);
                if(myrank == 1){
                    cout << "RANK " << myrank << " AT IDX: " << idx << endl; 
                } 
                MPI_Barrier(MPI_COMM_WORLD);
                
                int n = theta_argSort[idx]; 
                float phi = recv_particles.phi[n];

                if (phi > phi_cut_rough[haloIdx][0] && phi < phi_cut_rough[haloIdx][1]) {
                 
                    // of the particles surviving the rough cut, let's do a proper rotation 
                    // on them to find the true cutout membership, and return cluster-centric 
                    // angular coordinates
                    
                    // do coordinate rotation center halo at (r, 90, 0)
                    // B and k are the angle and axis of rotation, respectively,
                    // calculated near the beginning of this function
                    
                    //if(myrank == 0 and idx%500 == 0){
                    //    cout << "[x,y,z,t,p] = " << recv_particles.x[n] << ", " <<
                    //                                recv_particles.y[n] << ", " <<
                    //                                recv_particles.z[n] << ", " <<
                    //                                recv_particles.theta[idx] << ", " <<
                    //                                phi << "]" << endl;
                    //}

                    float tmp_v[] = {recv_particles.x[n], recv_particles.y[n], recv_particles.z[n]};
                    vector<float> v(tmp_v, tmp_v+3);
                    vector<float> v_rot;
                    v_rot = matVecMul(R[haloIdx], v);

                    // spherical coordinate transformation
                    float d = (float)sqrt(v_rot[0]*v_rot[0] + v_rot[1]*v_rot[1] + v_rot[2]*v_rot[2]);
                    float v_theta = acos(v_rot[2]/d) * 180.0 / PI * ARCSEC;
                    float v_phi;

                    // prevent NaNs on y-z plane
                    if(v_rot[0] == 0 && v_rot[1] > 0)
                        v_phi = 90.0 * ARCSEC;
                    else if(v_rot[0] == 0 && v_rot[1] < 0)
                        v_phi = -90.0 * ARCSEC;
                    else
                        v_phi = atan(v_rot[1]/v_rot[0]) * 180.0 / PI * ARCSEC; 
                    rough_cutout_size++;

                    // do final cut
                    if (v_theta > theta_cut[haloIdx][0] && v_theta < theta_cut[haloIdx][1] && 
                        v_phi > phi_cut[haloIdx][0] && v_phi < phi_cut[haloIdx][1] ) {

                        // get redshift from scale factor
                        float zz = aToZ(recv_particles.a[n]);
                        
                        // spherical corrdinate transform of rotated positions
                        w.theta.push_back(v_theta);
                        w.phi.push_back(v_phi);

                        // other columns
                        w.x.push_back(recv_particles.x[n]);
                        w.y.push_back(recv_particles.y[n]);
                        w.z.push_back(recv_particles.z[n]);
                        w.redshift.push_back(zz);
                        w.id.push_back(recv_particles.id[n]);
                        if(!positionOnly){ 
                            w.vx.push_back(recv_particles.vx[n]);
                            w.vy.push_back(recv_particles.vy[n]);
                            w.vz.push_back(recv_particles.vz[n]);
                            w.rotation.push_back(recv_particles.rotation[n]);
                            w.replication.push_back(recv_particles.replication[n]);
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
            thisRank_end = clock();

            MPI_Barrier(MPI_COMM_WORLD);
            
            stop = MPI_Wtime();
            duration = stop - start;
            if(myrank == 0 and timeit==true and printHalo){
                cout << "cutout computation time: " << duration << " s" << endl; 
            }
            cutout_times.push_back(duration);
            
            if(verbose == true and timeit == true and printHalo){
                
                // check load balancing (all ranks should have taken more or less the same amount of time here)
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
                    cout << "]" << endl;        

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
            //                          write out
            //
            ///////////////////////////////////////////////////////////////

            // time write out 
            MPI_Barrier(MPI_COMM_WORLD);
            start = MPI_Wtime();

            // define MPI file writing offset for the current rank --
            // This offset will be the sum of elements in all lesser ranks,
            // multiplied by the type size for each file    
            w.np_rough_count.clear();
            w.np_rough_count.resize(numranks);
            w.np_count.clear();
            w.np_count.resize(numranks);
            w.np_offset.clear();
            w.np_offset.push_back(0);

            // get number of elements in each ranks portion of cutout
            
            if(myrank == 0){
                cout << "RANK " << myrank << " CUTOUT COUNT: " << cutout_size << endl; 
                cout << "RANK " << myrank << " ROUGH CUTOUT COUNT: " << rough_cutout_size << endl; 
                cout << "RANK " << myrank << " CUTOUT COUNT SIZE: " << w.x.size() << endl;
            }
            MPI_Barrier(MPI_COMM_WORLD);
            if(myrank == 1){
                cout << "RANK " << myrank << " CUTOUT COUNT: " << cutout_size << endl; 
                cout << "RANK " << myrank << " ROUGH CUTOUT COUNT: " << rough_cutout_size << endl; 
                cout << "RANK " << myrank << " CUTOUT COUNT SIZE: " << w.x.size() << endl;
            } 

            MPI_Allgather(&cutout_size, 1, MPI_INT, 
                          &w.np_count[0], 1, MPI_INT, MPI_COMM_WORLD);
            MPI_Allgather(&rough_cutout_size, 1, MPI_INT, 
                          &w.np_rough_count[0], 1, MPI_INT, MPI_COMM_WORLD);
            int tot_cutout_size = 0; for(int i=0;i<numranks;++i){tot_cutout_size += w.np_count[i];}
            int tot_rough_cutout_size = 0; for(int i=0;i<numranks;++i){tot_rough_cutout_size += w.np_rough_count[i];}
            if(myrank == 0){
                cout << "(" << tot_rough_cutout_size << ") " << tot_cutout_size << 
                        " total particles in (rough) cutout" << endl;
            }
            
            // compute each ranks writing offset
            for(int j=1; j < numranks; ++j){
                w.np_offset.push_back(w.np_offset[j-1] + w.np_count[j-1]);
            }
            MPI_Barrier(MPI_COMM_WORLD); 
           
            // print out offset vector for verification
            if(myrank == 0 and printHalo){
                if(numranks < 20){
                    cout << "rank object rough counts: [";
                    for(int m=0; m < numranks; ++m){ cout << w.np_rough_count[m] << ","; }
                    cout << "]" << endl;
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

            // write... 
            MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(id_file_name.str().c_str()), 
                    MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &id_file);
            MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(x_file_name.str().c_str()), 
                    MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &x_file);
            MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(y_file_name.str().c_str()),
                    MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &y_file);
            MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(z_file_name.str().c_str()),
                    MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &z_file);
            MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(theta_file_name.str().c_str()),
                    MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &theta_file);
            MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(phi_file_name.str().c_str()),
                    MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &phi_file);
            MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(redshift_file_name.str().c_str()),
                    MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &redshift_file);
            
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
            
            MPI_File_close(&id_file);
            MPI_File_close(&x_file);
            MPI_File_close(&y_file);
            MPI_File_close(&z_file);
            MPI_File_close(&theta_file);
            MPI_File_close(&phi_file);
            MPI_File_close(&redshift_file);
            
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

        // print timing arrays out in a form ready for copy&paste into python...
        
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
        
        cout << "sort_times = np.array([";
        for(int hh = 0; hh < sort_times.size(); ++hh){
            cout << sort_times[hh];
            if(hh < sort_times.size()-1){ cout << ", "; }
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
