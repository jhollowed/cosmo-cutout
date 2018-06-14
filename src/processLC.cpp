#include "processLC.h"

using namespace std;
using namespace gio;

//////////////////////////////////////////////////////

struct Buffers_read {

    // Buffers to fill with data read from input LC
    vector<POSVEL_T> x;
    vector<POSVEL_T> y;
    vector<POSVEL_T> z;
    vector<POSVEL_T> vx;
    vector<POSVEL_T> vy;
    vector<POSVEL_T> vz;
    vector<POSVEL_T> a;
    vector<ID_T> id;
    vector<int> step;
    vector<int> rotation;
    vector<int32_t> replication;
};

struct Buffers_write {

    // Buffers to fill with data to write out to cut out
    vector<POSVEL_T> x;
    vector<POSVEL_T> y;
    vector<POSVEL_T> z;
    vector<POSVEL_T> vx;
    vector<POSVEL_T> vy;
    vector<POSVEL_T> vz;
    vector<POSVEL_T> a;
    vector<ID_T> id;
    vector<int> step;
    vector<int> rotation;
    vector<int32_t> replication;
    vector<float> theta;
    vector<float> phi;
    
    // Buffers to fill with MPI file writing offset values
    vector<int> np_count; // length of output data vecotrs for each rank
    vector<int> np_offset; // cumulative sum of np_count
};


//////////////////////////////////////////////////////
//
//                Cutout function
//                  Use Case 1
//           Custom theta - phi bounds
//
//////////////////////////////////////////////////////

void processLC(string dir_name, string out_dir, vector<string> step_strings, 
        vector<float> theta_cut, vector<float> phi_cut, int rank, int numranks){

    ///////////////////////////////////////////////////////////////
    //
    //                          Setup
    //
    ///////////////////////////////////////////////////////////////

    // instances of buffer structs at file header
    Buffers_read r;
    Buffers_write w;
    
    w.np_count.resize(numranks);
    w.np_offset.resize(numranks);

    // find all lc sub directories for each step in step_strings
    if(rank == 0){ cout << "\nReading directory: " << dir_name << endl; }
    vector<string> subdirs;
    getLCSubdirs(dir_name, subdirs);
    if(rank==0){ 
        cout << "Found subdirs:" << endl;
        for (vector<string>::const_iterator i = subdirs.begin(); i != subdirs.end(); ++i){
            cout << *i << ' ';
        }
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
    if(rank == 0){ cout << "Subdir prefix is: " << subdirPrefix << endl; }

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
        if(rank == 0){
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
            if(rank == 0){ cout << "Opening file: " << file_name_stream.str() << endl; }
            GenericIO GIO(MPI_COMM_WORLD, file_name_stream.str(), Method);
            GIO.openAndReadHeader(GenericIO::MismatchRedistribute);

            MPI_Barrier(MPI_COMM_WORLD);
            Np = GIO.readNumElems();
            if(rank == 0){
                cout << "Number of elements in lc step at rank " << rank << ": " << 
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
            r.step.resize(Np + GIO.requestedExtraSpace()/sizeof(int));
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
            GIO.addVariable("step", r.step, true); 
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
        r.step.resize(Np);
        r.rotation.resize(Np);
        r.replication.resize(Np);
        if(rank == 0){ cout<<"done resizing"<<endl; }

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
                ostringstream badDir;
                badDir << "Directory " << step_subdir.str() << " is non-empty";
                throw runtime_error(badDir.str());
            }
            if(rank == 0){ cout << "Entered subdir: " << step_subdir.str() << endl; }
        }
        // Otherwise, create the subdir
        else{
            mkdir(step_subdir.str().c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IXOTH);
            if(rank == 0){ cout << "Created subdir: " << step_subdir.str() << endl; }
        }

        // create binary files for cutout output
        MPI_File id_file;
        MPI_File a_file;
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
        MPI_Request a_req;
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
        ostringstream a_file_name;
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
        a_file_name << step_subdir.str() << "/a." << step << ".bin";
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

        if(rank == 0){ cout<<"starting to open files"<<endl; }

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
        MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(a_file_name.str().c_str()),
                MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &a_file);
        MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(rotation_file_name.str().c_str()),
                MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &rotation_file);
        MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(replication_file_name.str().c_str()),
                MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &replication_file);
        MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(theta_file_name.str().c_str()),
                MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &theta_file);
        MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(phi_file_name.str().c_str()),
                MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &phi_file);
        
        if(rank == 0){ cout<<"done opening files"<<endl; }

        ///////////////////////////////////////////////////////////////
        //
        //                         Do cutting
        //
        ///////////////////////////////////////////////////////////////

        if(rank == 0){ cout << "Converting positions..." << endl; }

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
                    w.a.push_back(r.a[n]);
                    w.id.push_back(r.id[n]);
                    w.rotation.push_back(r.rotation[n]);
                    w.replication.push_back(r.replication[n]);
                }
            }
        }

        ///////////////////////////////////////////////////////////////
        //
        //                   write out
        //
        ///////////////////////////////////////////////////////////////

        // define MPI file writing offset for the current rank --
        // This offset will be the sum of elements in all lesser ranks,
        // multiplied by the type size for each file    
        MPI_Barrier(MPI_COMM_WORLD);
        w.np_count.clear();
        w.np_offset.clear();
        w.np_offset[0] = 0;
        int cutout_size = int(w.a.size());
        
        // get number of elements in each ranks portion of cutout
        MPI_Allgather(&cutout_size, 1, MPI_INT, &w.np_count[0], 1, MPI_INT, 
                MPI_COMM_WORLD);
        
        // compute each ranks writing offset
        for(int j=1; j < numranks; ++j){
            w.np_offset[j] = w.np_offset[j-1] + w.np_count[j-1];
        }
        
        // print out offset vector for verification
        if(rank == 0){
            cout << "rank object counts: [";
            for(int m=0; m < numranks; ++m){ cout << w.np_count[m] << ","; }
            cout << "]" << endl;
            cout << "rank offsets: [";
            for(int m=0; m < numranks; ++m){ cout << w.np_offset[m] << ","; }
            cout << "]" << endl;
            cout << "Writing files..." << endl;
        }

        MPI_Offset offset_posvel = sizeof(POSVEL_T) * w.np_offset[rank];
        MPI_Offset offset_id = sizeof(ID_T) * w.np_offset[rank];
        MPI_Offset offset_float = sizeof(float) * w.np_offset[rank];
        MPI_Offset offset_int = sizeof(int) * w.np_offset[rank];
        MPI_Offset offset_int32 = sizeof(int32_t) * w.np_offset[rank];

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
        
        MPI_File_seek(a_file, offset_posvel, MPI_SEEK_SET);
        MPI_File_iwrite(a_file, &w.a[0], w.a.size(), MPI_FLOAT, &a_req);
        MPI_Wait(&a_req, MPI_STATUS_IGNORE);
        
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
        MPI_File_close(&a_file);
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
        vector<float> halo_pos, float boxLength, int rank, int numranks){

    ///////////////////////////////////////////////////////////////
    //
    //                          Setup
    //
    ///////////////////////////////////////////////////////////////

    // find all lc sub directories for each step in step_strings
    if(rank == 0){ cout << "\nReading directory: " << dir_name << endl; }
    vector<string> subdirs;
    getLCSubdirs(dir_name, subdirs);
    if(rank==0){ 
        cout << "Found subdirs:" << endl;
        for (vector<string>::const_iterator i = subdirs.begin(); i != subdirs.end(); ++i){
            cout << *i << ' ';
        }
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
    if(rank == 0){ cout << "Subdir prefix is: " << subdirPrefix << endl; }

    ///////////////////////////////////////////////////////////////
    //
    //                  Start coordinate rotation
    //
    ///////////////////////////////////////////////////////////////

    if(rank == 0){
        cout<< "\n\n---------- Setting up for coordinate rotation ----------" << endl; 
    }

    // do coordinate rotation to center each input halo at (r, 90, 0) in spherical coords...
    // halo_pos is a vector of halo positions, each having three components
    // k is a vector of vectors (the rotation axis per target halo)
    // B is a vector (the rotation angle per target halo)
    // theta is a vector of vectors (the altitude angular bounds per target halo)
    // phi is a vector of vectors (the azimuthal angular bounds per target halo)
    
    int numHalos = halo_pos.size()/3;
    vector<vector<float> > k(numHalos);
    vector<float> B(numHalos);
    vector<vector<float> > theta_cut(numHalos);
    vector<vector<float> > phi_cut(numHalos);

    for(int h=0; h<halo_pos.size(); h+=3){ 
        
        int haloIdx = h/3;
        
        // get next three values in halo_pos
        float tmp_pos[] = {halo_pos[h], halo_pos[h+1], halo_pos[h+2]};
        vector<float> this_halo_pos(tmp_pos, tmp_pos+3);
        
        // find distance magnitude (new rotated halo position)
        float halo_r = (float)sqrt(this_halo_pos[0]*this_halo_pos[0] + 
                                        this_halo_pos[1]*this_halo_pos[1] + 
                                        this_halo_pos[2]*this_halo_pos[2]);
        float tmp[] = {halo_r, 0, 0};
        vector<float> rotated_pos(tmp, tmp+3);
        if(rank == 0){ cout << "\nFinding axis of rotation to move (" << 
                       this_halo_pos[0]<< ", " << this_halo_pos[1]<< ", " << this_halo_pos[2]<< ") to (" <<
                       rotated_pos[0] << ", " << rotated_pos[1] << ", " << rotated_pos[2] <<
                       ")" << endl; }

        // get angle and axis of rotation -- this only needs to be calculated once for all
        // steps, and it will be used to rotate all other position vectors in the 
        // loops below
        normCross(this_halo_pos, rotated_pos, k[haloIdx]);
        B[haloIdx] = vecPairAngle(this_halo_pos, rotated_pos);

        if(rank == 0){ cout << "Rotation is " << B[haloIdx]*(180/PI) << "� about axis k = (" << 
                       k[haloIdx][0]<< ", " << k[haloIdx][1] << ", " << k[haloIdx][2] << ")" << endl; }

        // calculate theta_cut and phi_cut, in arcsec, given the specified boxLength
        float halfBoxLength = boxLength / 2.0;
        float dtheta = atan(halfBoxLength / halo_r);
        float dphi = dtheta;

        // calculate theta-phi bounds of cutout under coordinate rotation
        theta_cut[haloIdx].push_back( (PI/2 - dtheta) * 180.0/PI * ARCSEC );
        theta_cut[haloIdx].push_back( (PI/2 + dtheta) * 180.0/PI * ARCSEC );
        phi_cut[haloIdx].push_back( (0 - dphi) * 180.0/PI * ARCSEC );
        phi_cut[haloIdx].push_back( (0 + dphi) * 180.0/PI * ARCSEC );
        if(rank == 0){ 
            cout << "theta bounds set to: ";
            cout << theta_cut[haloIdx][0]/ARCSEC << "째 -> " << theta_cut[haloIdx][1]/ARCSEC <<"째"<< endl;
            cout << "phi bounds set to: ";
            cout << phi_cut[haloIdx][0]/ARCSEC << "째 -> " << phi_cut[haloIdx][1]/ARCSEC <<"째" << endl;
            cout << "theta-phi bounds result in box width of " << 
                    tan(dtheta) * halo_r * 2 << 
                    " Mpc at distance to halo of " << halo_r << endl << 
                    "        " << "= " << dtheta*2*180.0/PI << "째x" << dphi*2*180.0/PI << 
                    "째 field of view" << endl;
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
    
        // instances of buffer struct at file header for read in data
        Buffers_read r;

        // continue if this is the step at z=0 (lightcone volume zero)
        step =atoi(step_strings[i].c_str());
        if(step == 499){ continue;}

        // find header file
        if(rank == 0){
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
            if(rank == 0){ cout << "Opening file: " << file_name_stream.str() << endl; }
            GenericIO GIO(MPI_COMM_WORLD, file_name_stream.str(), Method);
            GIO.openAndReadHeader(GenericIO::MismatchRedistribute);

            MPI_Barrier(MPI_COMM_WORLD);
            Np = GIO.readNumElems();
            if(rank == 0){
                cout << "Number of elements in lc step at rank " << rank << ": " << 
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
            r.step.resize(Np + GIO.requestedExtraSpace()/sizeof(int));
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
            GIO.addVariable("step", r.step, true); 
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
        r.step.resize(Np);
        r.rotation.resize(Np);
        r.replication.resize(Np);
        if(rank == 0){ cout<<"done resizing"<<endl; }

        ///////////////////////////////////////////////////////////////
        //
        //           Create output files + start reading
        //
        ///////////////////////////////////////////////////////////////

        // again loop over all target halos
        for(int h=0; h<halo_pos.size(); h+=3){
            if(rank == 0){
                cout<< "\n---------- Working on halo "<< h <<"----------" << endl; 
            }
            int haloIdx = h/3;
            
            // instances of buffer struct at file header for output data
            Buffers_write w;
            w.np_count.resize(numranks);
            w.np_offset.resize(numranks);

            // open cutout subdirectory for this step...
            ostringstream step_subdir;
            step_subdir << out_dirs[haloIdx] << subdirPrefix << "Cutout" << step_strings[i];
            
            DIR *dir = opendir(step_subdir.str().c_str());
            struct dirent *d;
            int nf = 0;
            
            // if subdir already exists, make sure it's empty, because overwriting
            // binary files isn't always clean
            if(dir != NULL){
                while((d = readdir(dir)) != NULL){ if(++nf>2){ break;} }
                closedir(dir);

                if(nf > 2){
                    ostringstream badDir;
                    badDir << "Directory " << step_subdir.str() << " is non-empty";
                    throw runtime_error(badDir.str());
                }
                if(rank == 0){ cout << "Entered subdir: " << step_subdir.str() << endl; }
            }
            // Otherwise, create the subdir
            else{
                mkdir(step_subdir.str().c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IXOTH);
                if(rank == 0){ cout << "Created subdir: " << step_subdir.str() << endl; }
            }

            // create binary files for cutout output
            MPI_File id_file;
            MPI_File a_file;
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
            MPI_Request a_req;
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
            ostringstream a_file_name;
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
            a_file_name << step_subdir.str() << "/a." << step << ".bin";
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

            if(rank == 0){ cout<<"starting to open files"<<endl; }

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
            MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(a_file_name.str().c_str()),
                    MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &a_file);
            MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(rotation_file_name.str().c_str()),
                    MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &rotation_file);
            MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(replication_file_name.str().c_str()),
                    MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &replication_file);
            MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(theta_file_name.str().c_str()),
                    MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &theta_file);
            MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(phi_file_name.str().c_str()),
                    MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &phi_file);

            if(rank == 0){ cout<<"done opening files"<<endl; }
            
            ///////////////////////////////////////////////////////////////
            //
            //                         Do cutting
            //
            ///////////////////////////////////////////////////////////////

            if(rank == 0){ cout << "Converting positions..." << endl; }

            for (int n=0; n<Np; ++n) {

                // limit cutout to first octant for speed
                if (r.x[n] > 0.0 && r.y[n] > 0.0 && r.z[n] > 0.0){

                    // do coordinate rotation center halo at (r, 90, 0)
                    // B and k are the angle and axis of rotation, respectively,
                    // calculated near the beginning of this function
                    float tmp[] = {r.x[n], r.y[n], r.z[n]};
                    vector<float> v(tmp, tmp+3);
                    vector<float> v_rot;
                    rotate(k[haloIdx], B[haloIdx], v, v_rot);

                    // spherical coordinate rotation
                    float d = (float)sqrt(v_rot[0]*v_rot[0] + v_rot[1]*v_rot[1] + 
                                               v_rot[2]*v_rot[2]);
                    float v_theta = acos(v_rot[2]/d) * 180.0 / PI * ARCSEC;
                    float v_phi = atan(v_rot[1]/v_rot[0]) * 180.0 / PI * ARCSEC;

                    // do cut and push back data for objects in cutout
                    if (v_theta > theta_cut[haloIdx][0] && v_theta < theta_cut[haloIdx][1] && 
                            v_phi > phi_cut[haloIdx][0] && v_phi < phi_cut[haloIdx][1]) {
                        
                        // spherical corrdinate transform of rotated positions
                        w.theta.push_back(v_theta);
                        w.phi.push_back(v_phi);

                        // other columns
                        w.x.push_back(r.x[n]);
                        w.y.push_back(r.y[n]);
                        w.z.push_back(r.z[n]);
                        w.vx.push_back(r.vx[n]);
                        w.vy.push_back(r.vy[n]);
                        w.vz.push_back(r.vz[n]);
                        w.a.push_back(r.a[n]);
                        w.id.push_back(r.id[n]);
                        w.rotation.push_back(r.rotation[n]);
                        w.replication.push_back(r.replication[n]);
                    }
                }
            }

            ///////////////////////////////////////////////////////////////
            //
            //                   write out
            //
            ///////////////////////////////////////////////////////////////

            // define MPI file writing offset for the current rank --
            // This offset will be the sum of elements in all lesser ranks,
            // multiplied by the type size for each file    
            MPI_Barrier(MPI_COMM_WORLD);
            w.np_count.clear();
            w.np_offset.clear();
            w.np_offset[0] = 0;
            int cutout_size = int(w.a.size());
            
            // get number of elements in each ranks portion of cutout
            MPI_Allgather(&cutout_size, 1, MPI_INT, &w.np_count[0], 1, MPI_INT, 
                    MPI_COMM_WORLD);
            
            // compute each ranks writing offset
            for(int j=1; j < numranks; ++j){
                w.np_offset[j] = w.np_offset[j-1] + w.np_count[j-1];
            }
           
            // print out offset vector for verification
            if(rank == 0){
                cout << "rank object counts: [";
                for(int m=0; m < numranks; ++m){ cout << w.np_count[m] << ","; }
                cout << "]" << endl;
                cout << "rank offsets: [";
                for(int m=0; m < numranks; ++m){ cout << w.np_offset[m] << ","; }
                cout << "]" << endl;
                cout << "Writing files..." << endl;
            }

            MPI_Offset offset_posvel = sizeof(POSVEL_T) * w.np_offset[rank];
            MPI_Offset offset_id = sizeof(ID_T) * w.np_offset[rank];
            MPI_Offset offset_float = sizeof(float) * w.np_offset[rank];
            MPI_Offset offset_int = sizeof(int) * w.np_offset[rank];
            MPI_Offset offset_int32 = sizeof(int32_t) * w.np_offset[rank];

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
            
            MPI_File_seek(a_file, offset_posvel, MPI_SEEK_SET);
            MPI_File_iwrite(a_file, &w.a[0], w.a.size(), MPI_FLOAT, &a_req);
            MPI_Wait(&a_req, MPI_STATUS_IGNORE);
            
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
            MPI_File_close(&a_file);
            MPI_File_close(&theta_file);
            MPI_File_close(&phi_file);
            MPI_File_close(&rotation_file);
            MPI_File_close(&replication_file);
        }
    }
}
