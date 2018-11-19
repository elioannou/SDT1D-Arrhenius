#include "./sourceFiles/core.hpp"
#include "./sourceFiles/hllc.hpp"
#include "./sourceFiles/limiters.hpp"
#include "./sourceFiles/chemistry.hpp"


/*
Comments:
* MAKE SURE THERE IS A CORRESPONDING INITIAL FILE AND OUTPUT FOLDER.
* with non-dimensionalisation and arrhenious chemistry
* The non-dimensionalised state has rho=1, u=0, p=1.
*
*	WORKING! Fixed minor bug with wall heating. It was too big.
* 

*/

int main(){
 

	int M;
	float target_time,x0;
	double domainlength=1;
	double rho_il,u_il,p_il,c0_il,c3_il,c8_il,rho_ir,u_ir,p_ir,c0_ir,c3_ir,c8_ir;
	clock_t runtime;

	//INPUTING DATA
 
	std::ifstream initial_data;
	std::string simulation_name;
	simulation_name="detonation";
	initial_data.open(("./initial_data/initial_data_"+simulation_name).c_str());

	if (initial_data.is_open()){
		initial_data >> domainlength;
		initial_data >> M;
		initial_data >> x0;
		initial_data >> target_time;
		initial_data >> rho_il;
		initial_data >> u_il;
		initial_data >> p_il;
		initial_data >> c0_il;
		initial_data >> c3_il;
		initial_data >> c8_il;
		initial_data >> rho_ir;
		initial_data >> u_ir;
		initial_data >> p_ir;
		initial_data >> c0_ir;
		initial_data >> c3_ir;
		initial_data >> c8_ir;
	}

	else {
		std::cout << "ERROR: Inital Data file open error " << std::endl;
		exit(0);
	}
	//std::cout << "Enter time target" << std::endl;
	//std::cin >> target_time;

	//INITIALIZING DATA

	double Smax, dt, dx=domainlength/double(M), time=0;
	// PRIMITIVE VARIABLES
	double *rho= new double[M+4];
	double *u= new double[M+4];
	double *p= new double[M+4];
	double **c= new double*[M+4];
	for ( int i=0; i<M+4; i++){
		  c[i]= new double[9];
	}
	//CONSERVED VARIABLES
	double *U1= new double[M+4];
	double *U2= new double[M+4];
	double *U3= new double[M+4];
	double **U4= new double*[M+4];
	for ( int i=0; i<M+4; i++){
		U4[i]= new double[9];
	}
	//AUXILIARY VARIABLES
	double *e= new double[M+4];
	double *T=new double[M+4];  
	//D VECTORS
	double *D1= new double[M+4];//Use only M+2 elements
	double *D2= new double[M+4];
	double *D3= new double[M+4];
	double **D4= new double*[M+4];
	for ( int i=0; i<M+4; i++){
		  D4[i]= new double[9];
	}
	//BOUNDARY EXTRAPOLATED VALUES
	double *U1_l= new double[M+4];//Use only M+2 elements
	double *U2_l= new double[M+4];
	double *U3_l= new double[M+4];
	double **U4_l= new double*[M+4];
	for ( int i=0; i<M+4; i++){
		U4_l[i]= new double[9];
	}
	double *U1_r= new double[M+4];
	double *U2_r= new double[M+4];
	double *U3_r= new double[M+4];
	double **U4_r= new double*[M+4];
	for ( int i=0; i<M+4; i++){
		U4_r[i]= new double[9];
	}
	//AND BOUNDARY FLUXES
	double *F1_l= new double[M+4];//use only M+2 elements
	double *F2_l= new double[M+4];
	double *F3_l= new double[M+4];
	double** F4_l= new double*[M+4];
	for ( int i=0; i<M+4; i++){
		F4_l[i]= new double[9];
	}
	double *F1_r= new double[M+4];
	double *F2_r= new double[M+4];
	double *F3_r= new double[M+4];
	double** F4_r= new double*[M+4];
	for ( int i=0; i<M+4; i++){
		F4_r[i]= new double[9];
	}
	//FLUXES
	double *F1= new double[M+4];//use only M+1 elements
	double *F2= new double[M+4];
	double *F3= new double[M+4];
	double** F4= new double*[M+4];
	for ( int i=0; i<M+4; i++){
		F4[i]= new double[9];
	}

	int i,j,k=0; //counters. k is for timesteps
	
	for (i=0;i < x0*M+2;i++){
		rho[i]=rho_il;
		u[i]=u_il;
		p[i]=p_il;
		for (j=0; j < 9; j++){
			c[i][j]=0;
		}
		c[i][0]=c0_il;
		c[i][3]=c3_il;
		c[i][8]=c8_il;
	}

	for (i=x0*M+2;i < M+4;i++){
	    rho[i]=rho_ir;
	    u[i]=u_ir;
	    p[i]=p_ir;
	    for (j=0; j < 9; j++){
			c[i][j]=0;
		}
		c[i][0]=c0_ir;
		c[i][3]=c3_ir;
		c[i][8]=c8_ir;
	}
	  
	for (i=0;i<M+4;i++){
		F1[i]=F2[i]=F3[i]=0;
	    U1[i]=U2[i]=U3[i]=e[i]=T[i]=0;
	    D1[i]=D2[i]=D3[i]=0;
	    U1_l[i]=U2_l[i]=U3_l[i]=U1_r[i]=U2_r[i]=U3_r[i]=0;
	    for(j=0;j<9;j++){
			U4[i][j]=0;
			F4[i][j]=0;
			D4[i][j]=0;
			U4_l[i][j]=0;
			U4_r[i][j]=0;
	  }
	}

	// OUTPUTING INITIAL DATA on terminal and on file sim_data_
	
	std::cout << "\n\n|||||||||||||||||||||||||||||||||***SIMULATION: " << simulation_name << "***||||||||||||||||||||||||||||||||||" << std::endl;
	std::cout << "Domain length = " << domainlength << ". No. of cells = " << M << std::endl;
	std::cout << "Separation diaphragm at x0 = " << x0 << std::endl;
	std::cout << "Simulation target time = " << target_time << std::endl;
	std::cout << "|||||||||||||||||||||||||||||||***INITIAL CONDITIONS***||||||||||||||||||||||||||||||||||" << std::endl;
	std::cout << "LEFT\t" << "\t\t\tRIGHT" << std::endl;
	std::cout << "Density_l = " << rho[int(x0*M/2)] << "\t\t Density_r = " << rho[int(x0*M+(1-x0)*M/2)] << std::endl; 
	std::cout << "Velocity_l = " << u[int(x0*M/2)] << "\t\t Velocity_r = " << u[int(x0*M+(1-x0)*M/2)] << std::endl;
	std::cout <<  "Pressure_l = " << p[int(x0*M/2)] << "\t\t Pressure_r = " << p[int(x0*M+(1-x0)*M/2)] << std::endl; 
	std::cout <<  "H2_l = " << c[int(x0*M/2)][0] << "\t\t\t H2_r = " << c[int(x0*M+(1-x0)*M/2)][0] << std::endl;
	std::cout <<  "O2_l = " << c[int(x0*M/2)][3] << "\t\t\t O2_r = " << c[int(x0*M+(1-x0)*M/2)][3] << std::endl;
	std::cout <<  "Ar_l = " << c[int(x0*M/2)][8] << "\t\t\t Ar_r = " << c[int(x0*M+(1-x0)*M/2)][8] << std::endl;
	std::cout << "ps. above data come from cells no. " <<  int(x0*M/2) << " for left side and no. " << int(x0*M+(1-x0)*M/2) << " for right side." <<std::endl;
	std::cout << "|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n" << std::endl;

	  std::ofstream simulation_data;

    simulation_data.open(("./output/"+simulation_name+"/sim_data_"+simulation_name).c_str());
	if (simulation_data.is_open()) {
	simulation_data << "|||||||||||||||||||||||||||||||||***SIMULATION: " << simulation_name << "***||||||||||||||||||||||||||||||||||" << std::endl;
	simulation_data << "Domain length = " << domainlength << ". No. of cells = " << M << std::endl;
	simulation_data << "Separation diaphragm at x0 = " << x0 << std::endl;
	simulation_data << "Simulation target time = " << target_time << std::endl;
	simulation_data << "|||||||||||||||||||||||||||||||***INITIAL CONDITIONS***||||||||||||||||||||||||||||||||||" << std::endl;
	simulation_data << "LEFT\t" << "\t\t\tRIGHT" << std::endl;
	simulation_data << "Density_l = " << rho[int(x0*M/2)] << "\t\t Density_r = " << rho[int(x0*M+(1-x0)*M/2)] << std::endl; 
	simulation_data << "Velocity_l = " << u[int(x0*M/2)] << "\t\t Velocity_r = " << u[int(x0*M+(1-x0)*M/2)] << std::endl;
	simulation_data <<  "Pressure_l = " << p[int(x0*M/2)] << "\t\t Pressure_r = " << p[int(x0*M+(1-x0)*M/2)] << std::endl; 
	simulation_data <<  "H2_l = " << c[int(x0*M/2)][0] << "\t\t\t H2_r = " << c[int(x0*M+(1-x0)*M/2)][0] << std::endl;
	simulation_data <<  "O2_l = " << c[int(x0*M/2)][3] << "\t\t\t O2_r = " << c[int(x0*M+(1-x0)*M/2)][3] << std::endl;
	simulation_data <<  "Ar_l = " << c[int(x0*M/2)][8] << "\t\t\t Ar_r = " << c[int(x0*M+(1-x0)*M/2)][8] << std::endl;
	simulation_data << "ps. above data come from cells no. " <<  int(x0*M/2) << " for left side and no. " << int(x0*M+(1-x0)*M/2) << " for right side." <<std::endl;
	simulation_data << "|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;
	}

	else std::cout << "Unable to open file to output simulation data (at the start of the simulation)" << std::endl;
	simulation_data.close();

	//std::cout << "PRIMITIVE"<<std::endl;
	//print4Variables(M,dx,rho,u,p,c);

	int nsp=9; //number of species
	
	for (i=0;i<M+4;i++) converttoConserved(rho[i],u[i],p[i],c[i],U1[i],U2[i],U3[i],U4[i]);

	//std::cout << "CONSERVED"<<std::endl;
	//print4Variables(M,dx,U1,U2,U3,U4);

	// SETTINGS FOR OUTPUT
	int out_counter=1;
	double out_interval=0.05;
	bool out=false;
	
	 //TIME LOOP

	while ( time < target_time ) {
		k++;

		//WALL HEATING
  
		if (k == 100){
			for  (i=2;i<6;i++){
				// rho[i]=rho[i+4]; 
				 rho[i]=rho[i+5] - 0.0002; //copy neighbouring points only up to i<6
				// rho[i]=(0.319137)*((i-1)*dx) + 0.996467; //linear fit of points 6-17
				// rho[i]= -13.5875*((i-1)*dx)*((i-1)*dx) + 0.61876*((i-1)*dx) + 0.995034; //quadratic fit of points 6-7 also bad. better if fit extends to i<19
				rho[0]=rho[3];
				rho[1]=rho[2]; 	
				converttoConserved(rho[i],u[i],p[i],c[i],U1[i],U2[i],U3[i],U4[i]);
			}
			time=0.0;
		}    

   
		//CALCULATE D SLOPE VECTORS

		slopeLimiter(D1,D2,D3,D4,U1,U2,U3,U4,M,k,0); //use conserved variables. 0 for minmod 1 for superbee

		//CALCULATE BOUNDARY EXTRAPOLATED VALUES
    
		for (i=1;i<M+3;i++){
			U1_l[i]=U1[i]-0.5*D1[i];
			U2_l[i]=U2[i]-0.5*D2[i];
			U3_l[i]=U3[i]-0.5*D3[i];
			for (j=0;j<nsp;j++){
				U4_l[i][j]=U4[i][j]-0.5*D4[i][j];
			}
			U1_r[i]=U1[i]+0.5*D1[i];
			U2_r[i]=U2[i]+0.5*D2[i];
			U3_r[i]=U3[i]+0.5*D3[i];
			for (j=0;j<nsp;j++){
				U4_r[i][j]=U4[i][j]+0.5*D4[i][j];
			}
		}
    

		for (i=0;i<M+4;i++){
			if (U1_l[i] != U1_l[i] ) std::cout << "ERROR:(after extrapolating bounadry values) U1 rho_l is a NaN" << std::endl;
			if (U2_l[i] != U2_l[i] ) std::cout << "ERROR:(after extrapolating bounadry values) momentum_l is a NaN" << std::endl;
			if (U3_l[i] != U3_l[i] ) std::cout << "ERROR:(after extrapolating bounadry values) energy_l is a NaN" << std::endl;
			//if (U4_l[i] != U4_l[i] ) std::cout << "ERROR:(after extrapolating bounadry values) U4 rho*c_l is a NaN" << std::endl;
			if ( U1_l[i] < 0 ) std::cout << "ERROR:(after extrapolating bounadry values) U1 rho_l is negative at timestep = " << k << " and x = " << i*dx   << std::endl;
			if ( U3_l[i] < 0 ) std::cout << "ERROR:(after extrapolating bounadry values) energy_l is negative at timestep = " << k << " and i = " << i << std::endl;
			//if ( U4_l[i] < 0 ) std::cout << "ERROR:(after extrapolating bounadry values) U4 rho*c_l is negative at timestep = " << j << " and i = " << i << std::endl;
			if (U1_r[i] != U1_r[i] ) std::cout << "ERROR:(after extrapolating bounadry values) rho_r is a NaN" << std::endl;
			if (U2_r[i] != U2_r[i] ) std::cout << "ERROR:(after extrapolating bounadry values) momentum_r is a NaN" << std::endl;
			if (U3_r[i] != U3_r[i] ) std::cout << "ERROR:(after extrapolating bounadry values) energy_r is a NaN" << std::endl;
			//if (U4_r[i] != U4_r[i] ) std::cout << "ERROR:(after extrapolating bounadry values) U4 rho*c_r is a NaN" << std::endl;
			if ( U1_r[i] < 0 ) std::cout << "ERROR:(after extrapolating bounadry values) rho_r is negative at timestep = " << k << " and x = " << i*dx   << std::endl;
			if ( U3_r[i] < 0 ) std::cout << "ERROR:(after extrapolating bounadry values) energy_r is negative at timestep = " << k << " and i = " << i << " Energy = " << U3[i] << std::endl;
			//if ( U4_r[i] < 0 ) std::cout << "ERROR:(after extrapolating bounadry values) U4 rho*c_r is negative at timestep = " << j << " and i = " << i << std::endl;
		}
 
		//CALCULATE Smax    

		double umax=0;
		Smax=0;
		for (i=0;i<M+4;i++){
			umax=fabs(u[i])+soundSpeed(p[i],rho[i]);
			if (umax > Smax) {
			Smax=umax;
			}
		}

		//SETTING CFL CONDITION AND OUTPUT TIMES

		double cfl=0.95;
		//if (time < 0.02) cfl=0.5;

		if (k < 6) dt=0.2*cfl*dx/Smax;
		else dt=cfl*dx/Smax;

		if ( time > 0.90 && time < 0.90+0.001 ) { // BAD HACK HERE
			out_interval=0.01;
			out_counter=91;
		}

		if ( (time+dt) > out_counter*out_interval) {
			dt=out_counter*out_interval-time;
			out_counter++;
			out=true;
		}

		time+=dt;

		//std::cout << "dt = " << dt << std::endl;
    
		if ( dt <= 0 ){
			std::cout << "ERROR: dt = " << dt << " negative or zero at timestep = " << k << std::endl; 
			exit(1);
		}

		//EVOLVE BOUNDARY EXTRAPOLATED VALUES
 
		for (i=1;i<M+3;i++){
			calculateFluxes(U1_l[i],U2_l[i],U3_l[i],U4_l[i],F1_l[i],F2_l[i],F3_l[i],F4_l[i]);
			calculateFluxes(U1_r[i],U2_r[i],U3_r[i],U4_r[i],F1_r[i],F2_r[i],F3_r[i],F4_r[i]);
		}
    
		for (i=1;i<M+3;i++){
			U1_l[i]+=0.5*(dt/dx)*(F1_l[i]-F1_r[i]);
			U2_l[i]+=0.5*(dt/dx)*(F2_l[i]-F2_r[i]);
			U3_l[i]+=0.5*(dt/dx)*(F3_l[i]-F3_r[i]);
			for (j=0;j<nsp;j++){
				U4_l[i][j]+=0.5*(dt/dx)*(F4_l[i][j]-F4_r[i][j]);
			}	
			U1_r[i]+=0.5*(dt/dx)*(F1_l[i]-F1_r[i]);
			U2_r[i]+=0.5*(dt/dx)*(F2_l[i]-F2_r[i]);
			U3_r[i]+=0.5*(dt/dx)*(F3_l[i]-F3_r[i]);
			for (j=0;j<nsp;j++){
				U4_r[i][j]+=0.5*(dt/dx)*(F4_l[i][j]-F4_r[i][j]);
			}
		}

		for (i=0;i<M+4;i++){
			if (U1_l[i] != U1_l[i] ) std::cout << "ERROR:(after evolving bounadry values) rho_l is a NaN" << std::endl;
			if (U2_l[i] != U2_l[i] ) std::cout << "ERROR:(after evolving bounadry values) momentum_l is a NaN" << std::endl;
			if (U3_l[i] != U3_l[i] ) std::cout << "ERROR:(after evolving bounadry values) energy_l is a NaN" << std::endl;
			if ( U1_l[i] < 0 ) std::cout << "ERROR:(after evolving bounadry values) rho_l is negative at timestep = " << k << " and x = " << i   << std::endl;
			if ( U3_l[i] < 0 ) std::cout << "ERROR:(after evolving bounadry values) energy_l is negative at timestep = " << k << " and x = " << i << " Energy = " << U3_l[i] << std::endl;
			if (U1_r[i] != U1_r[i] ) std::cout << "ERROR:(after evolving bounadry values) rho_r is a NaN" << std::endl;
			if (U2_r[i] != U2_r[i] ) std::cout << "ERROR:(after evolving bounadry values) momentum_r is a NaN" << std::endl;
			if (U3_r[i] != U3_r[i] ) std::cout << "ERROR:(after evolving bounadry values) energy_r is a NaN" << std::endl;
			if ( U1_r[i] < 0 ) std::cout << "ERROR:(after evolving bounadry values) rho_r is negative at timestep = " << k << " and x = " << i << std::endl;
			if ( U3_r[i] < 0 ) std::cout << "ERROR:(after evolving bounadry values) energy_r is negative at timestep = " << k << " and x = " << i << " Energy = " << U3_r[i] << std::endl;
		}

//This probably not needed

		for (i=0;i<M+4;i++)  converttoPrimitive(rho[i], u[i], p[i], c[i], U1[i], U2[i], U3[i], U4[i]);
    
		for (i=0;i<M+4;i++){
			if (rho[i] != rho[i] ) std::cout << "ERROR:(before hllc) rho is nan at timestep = " << k << " and x = " << i*dx   << std::endl;
			if (u[i] != u[i] ) std::cout << "ERROR:(before hllc) velocity is nan at timestep = " << k << " and x = " << i*dx   << std::endl;
			if (p[i] != p[i] )std::cout << "ERROR:(before hllc) pressure is nan at timestep = " << k << " and x = " << i*dx   << std::endl; 
			if ( rho[i] <= 0 )std::cout << "ERROR:(before hllc) rho is negative at timestep = " << k << " and x = " << i*dx   << std::endl; 
			if ( p[i] <= 0 ) std::cout << "ERROR:(before hllc) pressure is negative at timestep = " << k << " and x = " << i*dx  << " Pressure = " << p[i] << std::endl; 
		}

		

    	//std::cout << "CONSERVED(before hllc) LEFT"<<std::endl;
    	//print4Variables(M,dx,U1_l,U2_l,U3_l,U4_l);

    	//std::cout << "CONSERVED(before hllc) RIGHT"<<std::endl;
		//print4Variables(M,dx,U1_r,U2_r,U3_r,U4_r);

		//CALCULATE FLUXES ARRAY

		for (i=1;i<M+2;i++){ 
			hllcFlux(U1_r[i],U2_r[i],U3_r[i],U4_r[i],U1_l[i+1],U2_l[i+1],U3_l[i+1],U4_l[i+1],F1[i],F2[i],F3[i],F4[i],i,k);
		}

  
		//ADVANCING IN TIME

		for (i=2;i<M+2;i++){
			U1[i] += (dt/dx)*(F1[i-1]-F1[i]);
			U2[i] += (dt/dx)*(F2[i-1]-F2[i]);
			U3[i] += (dt/dx)*(F3[i-1]-F3[i]);
			for (j=0;j<nsp;j++) {
				U4[i][j] += (dt/dx)*(F4[i-1][j]-F4[i][j]);
			}
		}

		for (i=0;i<M+4;i++){
			if (U1[i] != U1[i] ) std::cout << "ERROR:(after advancing in time) rhoU is a NaN" << std::endl;
			if (U2[i] != U2[i] ) std::cout << "ERROR:(after advancing in time) momentum is a NaN" << std::endl;
			if (U3[i] != U3[i] ) std::cout << "ERROR:(after advancing in time) energy is a NaN" << std::endl;
			if ( U1[i] < 0 ) std::cout << "ERROR:(after advancing in time) rhoU is negative at timestep = " << k << " and x = " << i*dx   << std::endl;
			if ( U3[i] < 0 ) std::cout << "ERROR:(after advancing in time) energy is negative at timestep = " << k << " and i = " << i << " Energy = " << U3[i] << std::endl;
		}

		//EVOLVING BECAUSE OF SOURCE TERM
    
		for (i=0;i<M+4;i++)  converttoPrimitive(rho[i], u[i], p[i], c[i], U1[i], U2[i], U3[i], U4[i]); //c[i] and U4[i] are arrays

		for (i=2;i<M+2;i++){
			T[i]=p[i]/rho[i];
		}

		if ( k >= 100) arrheniusODE(c,T,dt,M,k);

		for (i=2;i<M+2;i++){
			p[i]=rho[i]*T[i];     
		}
		
		// UPDATE BOUNDARY CONDITIONS

		//(REFLECTIVE)
		rho[1]=rho[2];
		u[1]=-u[2];
		p[1]=p[2];
		for (j=0;j<nsp;j++) {
			c[1][j]=c[2][j];
			c[0][j]=c[3][j];
		}
		rho[0]=rho[3];
		u[0]=-u[3];
		p[0]=p[3];
		//(TRANSMISSIVE)
		rho[M+2]=rho[M+1];
		u[M+2]=u[M+1];
		p[M+2]=p[M+1];
		for (j=0;j<nsp;j++) {
			c[M+2][j]=c[M+1][j];
			c[M+3][j]=c[M][j];
		}
		rho[M+3]=rho[M];
		u[M+3]=u[M];
		p[M+3]=p[M];


		for (i=0;i<M+4;i++){
			if (rho[i] != rho[i] ) std::cout << "ERROR:(after source term) rho is nan at timestep = " << k << " and x = " << i*dx   << std::endl;
			if (u[i] != u[i] ) std::cout << "ERROR:(after source term) velocity is nan at timestep = " << k << " and x = " << i*dx   << std::endl;
			if (p[i] != p[i] )std::cout << "ERROR:(after source term) pressure is nan at timestep = " << k << " and x = " << i*dx   << std::endl; 
			if ( rho[i] < 0 )	std::cout << "ERROR:(after source term) rho is negative at timestep = " << k << " and x = " << i*dx   << std::endl; 
			if ( p[i] < 0 ) std::cout << "ERROR:(after source term) pressure is negative at timestep = " << k << " and x = " << i*dx  << " Pressure = " << p[i] << std::endl; 
		}
    
		//updating conservative variables. is it needed?
 
		for (i=0;i<M+4;i++)  converttoConserved(rho[i], u[i], p[i], c[i], U1[i], U2[i], U3[i], U4[i]);
 
		for (i=0;i<M+4;i++){
			if (U1[i] != U1[i] ) std::cout << "ERROR:(after source term) rhoU is a NaN" << std::endl;
			if (U2[i] != U2[i] ) std::cout << "ERROR:(after source term) momentum is a NaN" << std::endl;
			if (U3[i] != U3[i] ) std::cout << "ERROR:(after source term) energy is a NaN" << std::endl;
			//if (U4[i] != U4[i] ) std::cout << "ERROR:(after source term) rho*c is a NaN. rho[" << i << "] = " << rho[i] << " c[" << i << "] = " << c[i] <<  std::endl;
			if ( U1[i] < 0 ) std::cout << "ERROR:(after source term) rhoU is negative at timestep = " << k << " and x = " << i*dx   << std::endl;
			if ( U3[i] < 0 ) std::cout << "ERROR:(after source term) energy is negative at timestep = " << k << " and i = " << i << " Energy = " << U3[i] << std::endl;
			//if ( U4[i] < 0 ) std::cout << "ERROR:(after source term) rho*c is negative at timestep = " << j << " rho[" << i << "] = " << rho[i] << " c[" << i << "] = " << c[i] <<  std::endl;
		}


		//PERIODIC EXPORTING
		if ( out == true || k==100 || k==99) {
			for (i=0;i<M+4;i++)  e[i]=p[i]/rho[i]/0.4; //Calculating internal energy. (is it dedimensionalised?)
			runtime = clock();
			exportResults(simulation_name,M,k,time,dx,rho,u,p,c,e,T,runtime);
			std::cout << "Progress: " << 100*time/target_time << "% (running time: " << runtime/CLOCKS_PER_SEC/60.0 << " minutes)" << ". Output at t = " << time << std::endl;
			out=false;
		}
  
	}

	std::cout << "\n|||||||||||||||||||||||||||||||||***SIMULATION: " << simulation_name << "***||||||||||||||||||||||||||||||||||" << std::endl;
	std::cout << "DONE in " << runtime/CLOCKS_PER_SEC << " seconds. (" << k << " timesteps). Output in folder with simulation name." << std::endl;
	std::cout << "|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n" << std::endl;
  
	//VARIABLES
	delete[] rho;
	delete[] p;
	delete[] u;
	for (i=0;i<M+4;i++) delete[] c[i];
	delete[] c;
	delete[] U1;
	delete[] U2;
	delete[] U3;
	for (i=0;i<M+4;i++) delete[] U4[i];
	delete[] U4;
	delete[] e;
	delete[] T;
	//FLUXES
	delete[] F1;
	delete[] F2;
	delete[] F3;
	for (i=0;i<M+4;i++) delete[] F4[i];
	delete[] F4;
	//SLOPES
	delete[] D1;
	delete[] D2;
	delete[] D3;
	for (i=0;i<M+4;i++) delete[] D4[i];
	delete[] D4;
	//BOUNDARY EXTRAPOLATED VALUES
	delete[] U1_l;
	delete[] U2_l;
	delete[] U3_l;
    for (i=0;i<M+4;i++) delete[] U4_l[i];
	delete[] U4_l;
	delete[] U1_r;
	delete[] U2_r;
	delete[] U3_r;
    for (i=0;i<M+4;i++) delete[] U4_r[i];
	delete[] U4_r;
	//AND BOUNDARY FLUXES
	delete[] F1_l;
	delete[] F2_l;
	delete[] F3_l;
	for (i=0;i<M+4;i++) delete[] F4_l[i];
	delete[] F4_l;
	delete[] F1_r;
	delete[] F2_r;
	delete[] F3_r;
	for (i=0;i<M+4;i++) delete[] F4_r[i];
	delete[] F4_r;
	
	return 0;
}
