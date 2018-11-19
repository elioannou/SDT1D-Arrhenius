// #include "cantera/zerodim.h"
// #include "cantera/IdealGasMix.h"

void arrheniusODE(double** c,double T[],double DT,int size, int timestep)
{
	const double epsilon=(1.0/14.99),Q=1.24,gamma=1.4;
	double tr,dt,beta,c_n,T_n;
	const int s=10;

	for (int i=2;i<size+2;i++){

		if ( c[i][0] > 1e-7){

			//de-dimensionalizing time
			tr=( T[i]*T[i]*exp( (1.0/epsilon)*( (1.0/T[i]) - 1) ) ) / ( c[i][0]*gamma );
			if ( tr != tr ) std::cout << "ODE ERROR: tr nan at timestep = " << timestep << " position " << i << " c[" << i << "] = " << c[i][0] << " T[" << i << "] = " << T[i] << std::endl;
			//LOCAL TIMESTEP
			dt= (DT/tr)/s;
			if ( dt != dt ) std::cout << "ODE ERROR: dt nan at timestep = " << timestep << " position " << i<< std::endl;
			//useful constant
			beta= gamma*Q*c[i][0]/T[i];
			if ( beta != beta ) std::cout << "ODE ERROR: beta nan at timestep = " << timestep << " position " << i << std::endl;
			//intializing
			T_n=1;
			c_n=1;

			//solver loop
			for (int j=0;j<s;j++){
				c_n=c_n*exp( -dt*epsilon*T[i]*(1.0/beta)*exp( ( 1.0/(epsilon*T[i]) ) * ( 1-(1/T_n) ) ) );
				T_n = ( 1 + beta*(1-c_n) );
				if ( T_n != T_n ) std::cout << "ODE ERROR: Tn nan at iteration " << j << " and timestep = " << timestep << " position " << i << std::endl;
				if ( c_n != c_n ) std::cout << "ODE ERROR: cn nan at iteration " << j << " and timestep = " << timestep << " position " << i << std::endl;    
		 	}

			//dimensionalizing
			T[i]=T[i]*T_n;
			c[i][0]=c[i][0]*c_n;
		}

	}
}
/*
void canteraAdvance(Cantera::IdealGasMix &gas,double t) //UNITS are K, Pa, massfractions, seconds
{	
	if ( t > 0) {
			
		Cantera::Reactor r;
		r.insert(gas);
		
		Cantera::ReactorNet sim;
		sim.addReactor(&r, false);
		
		sim.advance(t);
	}
	
	else std::cout << "ERROR: Negative time" << std::endl;
    
}
*/
