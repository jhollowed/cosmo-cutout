#include "processLC.h"

using namespace std;
using namespace gio;

//////////////////////////////////////////////////////

struct Buffers {
    
    // From LC output
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

    // New data columns
    vector<POSVEL_T> x_out;
    vector<POSVEL_T> y_out;
    vector<POSVEL_T> z_out;
    vector<POSVEL_T> vx_out;
    vector<POSVEL_T> vy_out;
    vector<POSVEL_T> vz_out;
    vector<POSVEL_T> a_out;
    vector<ID_T> id_out;
    vector<int> step_out;
    vector<int> rotation_out;
    vector<int32_t> replication_out;
    vector<float> theta;
    vector<float> phi;
    vector<float> thetaRot;
    vector<float> phiRot;
};


//////////////////////////////////////////////////////
//
//                Cutout function
//             	    Use Case 1
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
               vector<float> halo_pos, float boxLength, int rank, int numranks){

    ///////////////////////////////////////////////////////////////
    //
    //                          Setup
    //
    ///////////////////////////////////////////////////////////////
    
    Buffers b; // arrays to hold quantities
    vector<int> valid_indices; // array to hold indices of all objects within the cutout
    vector<int> np_count; // sum of pIdx for each rank
    vector<int> np_offset; // cumulative sum of np_count
    np_count.resize(numranks);
    np_offset.resize(numranks);

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
       
        // create gio reader and open file header in new scope
        size_t Np = 0;
	unsigned Method = GenericIO::FileIOPOSIX;
	const char *EnvStr = getenv("GENERICIO_USE_MPIIO");
	if(EnvStr && string(EnvStr) == "1")
	    Method = GenericIO::FileIOMPI;	

        {
        cout << "Opening file: " << file_name_stream.str() << endl;
        GenericIO GIO(MPI_COMM_WORLD, file_name_stream.str(), Method);
        GIO.openAndReadHeader(GenericIO::MismatchRedistribute);
	
	MPI_Barrier(MPI_COMM_WORLD);
        Np = GIO.readNumElems();
	cout << "Np: " << Np << endl;
	
	// resize buffers	
	b.x.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
	b.y.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
	b.z.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
	b.vx.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
	b.vy.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
	b.vz.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
	b.a.resize(Np + GIO.requestedExtraSpace()/sizeof(POSVEL_T));
	b.id.resize(Np + GIO.requestedExtraSpace()/sizeof(ID_T));
	b.step.resize(Np + GIO.requestedExtraSpace()/sizeof(int));
	b.rotation.resize(Np + GIO.requestedExtraSpace()/sizeof(int));
	b.replication.resize(Np + GIO.requestedExtraSpace()/sizeof(int32_t));
	
	// do reading
        GIO.addVariable("x", &b.x[0]); 
        GIO.addVariable("y", &b.y[0]); 
        GIO.addVariable("z", &b.z[0]); 
        GIO.addVariable("vx", &b.vx[0]); 
        GIO.addVariable("vy", &b.vy[0]); 
        GIO.addVariable("vz", &b.vz[0]); 
        GIO.addVariable("a", &b.a[0]); 
        GIO.addVariable("step", &b.step[0]); 
        GIO.addVariable("id", &b.id[0]); 
        GIO.addVariable("rotation", &b.rotation[0]); 
        GIO.addVariable("replication", &b.replication[0]);

	GIO.readData(); 
	}
	
	// resize again to remove reader extra space
        b.x.resize(Np);
        b.y.resize(Np);
        b.z.resize(Np);
        b.vx.resize(Np);
        b.vy.resize(Np);
        b.vz.resize(Np);
        b.a.resize(Np);
        b.id.resize(Np);
        b.step.resize(Np);
        b.rotation.resize(Np);
        b.replication.resize(Np);
        cout<<"done resizing"<<endl;
         
        ///////////////////////////////////////////////////////////////
        //
        //           Create output files + start reading
        //
        ///////////////////////////////////////////////////////////////
        
        // create cutout subdirectory for this step
        ostringstream step_subdir;
        step_subdir << out_dir << subdirPrefix << "Cutout" << step_strings[i];
        mkdir(step_subdir.str().c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IXOTH);
        cout << "Created subdir: " << step_subdir.str() << endl;
        
	// create binary files for cutout output
        MPI_File id_file;
	MPI_File theta_file;
        MPI_File phi_file;
        MPI_File a_file;
        MPI_File x_file;
        MPI_File y_file;
        MPI_File z_file;
        MPI_File vx_file;
        MPI_File vy_file;
        MPI_File vz_file;
        MPI_File rotation_file;
        MPI_File replication_file;
        MPI_File thetaRot_file;
        MPI_File phiRot_file;
        
	MPI_Request id_req;
	MPI_Request theta_req;
        MPI_Request phi_req;
        MPI_Request a_req;
        MPI_Request x_req;
        MPI_Request y_req;
        MPI_Request z_req;
        MPI_Request vx_req;
        MPI_Request vy_req;
        MPI_Request vz_req;
        MPI_Request rotation_req;
        MPI_Request replication_req;
        MPI_Request thetaRot_req;
        MPI_Request phiRot_req;

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
	MPI_File_open(MPI_COMM_WORLD, id_file_name.str().c_str(), 
		      MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &id_file);
	MPI_File_open(MPI_COMM_WORLD, x_file_name.str().c_str(), 
		      MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &x_file);
	MPI_File_open(MPI_COMM_WORLD, y_file_name.str().c_str(), 
		      MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &y_file);
	MPI_File_open(MPI_COMM_WORLD, z_file_name.str().c_str(), 
		      MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &z_file);
	MPI_File_open(MPI_COMM_WORLD, vx_file_name.str().c_str(), 
		      MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &vx_file);
	MPI_File_open(MPI_COMM_WORLD, vy_file_name.str().c_str(), 
		      MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &vy_file);
	MPI_File_open(MPI_COMM_WORLD, vz_file_name.str().c_str(), 
		      MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &vz_file);
	MPI_File_open(MPI_COMM_WORLD, a_file_name.str().c_str(), 
		      MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &a_file);
	MPI_File_open(MPI_COMM_WORLD, rotation_file_name.str().c_str(), 
		      MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &rotation_file);
	MPI_File_open(MPI_COMM_WORLD, replication_file_name.str().c_str(), 
		      MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &replication_file);
	MPI_File_open(MPI_COMM_WORLD, theta_file_name.str().c_str(), 
		      MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &theta_file);
	MPI_File_open(MPI_COMM_WORLD, phi_file_name.str().c_str(), 
		      MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &phi_file);
	MPI_File_open(MPI_COMM_WORLD, thetaRot_file_name.str().c_str(), 
		      MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &thetaRot_file);
	MPI_File_open(MPI_COMM_WORLD, phiRot_file_name.str().c_str(), 
		      MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &phiRot_file);
        cout<<"done opening files"<<endl;


        ///////////////////////////////////////////////////////////////
        //
        //                         Do cutting
        //
        ///////////////////////////////////////////////////////////////
        
         
	cout << "Converting positions..." << endl;

	for (int n=0; n<Np; ++n) {
	    
	    // limit cutout to first octant for speed
	    if (b.x[n] > 0.0 && b.y[n] > 0.0 && b.z[n] > 0.0){
	
		// do coordinate rotation center halo at (r, 90, 0)
		// B and k are the angle and axis of rotation, respectively,
		// calculated near the beginning of this function
		float tmp[] = {b.x[n], b.y[n], b.z[n]};
		vector<float> v(tmp, tmp+3);
		vector<float> v_rot;
		rotate(k, B, v, v_rot);

		// spherical coordinate rotation
		float r = (float)sqrt(v_rot[0]*v_rot[0] + v_rot[1]*v_rot[1] + 
				      v_rot[2]*v_rot[2]);
		float v_theta = acos(v_rot[2]/r) * 180.0 / PI * ARCSEC;
		float v_phi = atan(v_rot[1]/v_rot[0]) * 180.0 / PI * ARCSEC;
			
		// do cut and push back data for objects in cutout
		if (v_theta > theta_cut[0] && v_theta < theta_cut[1] && 
		    v_phi > phi_cut[0] && v_phi < phi_cut[1]) {
		    
		    // spherical coordiante transform of original (pre-rotate) position  
		    b.theta.push_back(acos(b.z[n]/r) * 180.0 / PI * ARCSEC);
		    b.phi.push_back(atan(b.y[n]/b.x[n]) * 180.0 / PI * ARCSEC);
		    
		    // spherical corrdinate transform of rotated positions
		    b.thetaRot.push_back(v_theta);
	 	    
                    // other columns
	 	    b.x_out.push_back(b.x[n]);
	 	    b.y_out.push_back(b.y[n]);
	 	    b.z_out.push_back(b.z[n]);
	 	    b.vx_out.push_back(b.vx[n]);
	 	    b.vy_out.push_back(b.vy[n]);
	 	    b.vz_out.push_back(b.vz[n]);
	 	    b.a_out.push_back(b.a[n]);
	 	    b.id_out.push_back(b.id[n]);
	 	    b.rotation_out.push_back(b.rotation[n]);
	 	    b.replication_out.push_back(b.replication[n]);
                }
            }
        }
        
	///////////////////////////////////////////////////////////////
        //
        //                 	 write out
        //
        ///////////////////////////////////////////////////////////////

	// define MPI file writing offset for the current rank --
	// This offset will be the sum of elements in all lesser ranks,
	// multiplied by the type size for each file	
	MPI_Barrier(MPI_COMM_WORLD);	
	np_offset[0] = 0;
	int cutout_size = int(b.a_out.size());

	MPI_Alltoall(&cutout_size, 1, MPI_INT, &np_count[0], 1, MPI_INT, 
		     MPI_COMM_WORLD);
	
	for(int j=1; j <= numranks; ++j){
            np_offset[j] = np_offset[j-1] + np_count[j-1];
	}

	MPI_Offset offset_posvel = sizeof(POSVEL_T)*np_offset[rank];
	MPI_Offset offset_id = sizeof(ID_T)*np_offset[rank];
	MPI_Offset offset_float = sizeof(float)*np_offset[rank];
	MPI_Offset offset_int = sizeof(int)*np_offset[rank];
	MPI_Offset offset_int32 = sizeof(int32_t)*np_offset[rank];
	
	// write
	MPI_File_seek(id_file, offset_id, MPI_SEEK_SET);
	MPI_File_iwrite(id_file, &b.id_out[0], b.id_out.size(), MPI_INT64_T, &id_req);
	MPI_Wait(&id_req, MPI_STATUS_IGNORE);
	
	MPI_File_seek(x_file, offset_posvel, MPI_SEEK_SET);
	MPI_File_iwrite(x_file, &b.x_out[0], b.x_out.size(), MPI_FLOAT, &x_req);
	MPI_Wait(&x_req, MPI_STATUS_IGNORE);
	
	MPI_File_seek(a_file, offset_posvel, MPI_SEEK_SET);
	MPI_File_iwrite(a_file, &b.a_out[0], b.a_out.size(), MPI_FLOAT, &a_req);
	MPI_Wait(&a_req, MPI_STATUS_IGNORE);
		 
        MPI_File_close(&id_file);
        MPI_File_close(&x_file);
        MPI_File_close(&a_file);
    }
}
