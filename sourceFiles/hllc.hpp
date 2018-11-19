void hllcFlux(const double U1_l,const double U2_l,const double U3_l,const double U4_l[], const double U1_r,const double U2_r,const double U3_r,const double U4_r[], double &F1, double &F2, double &F3, double F4[], int x, int timestep)
{
	const double gamma=1.4;

	//COMPUTING PRIMITIVE VARIABLES
  
	double rho_l,u_l,p_l,rho_r,u_r,p_r;
	double c_l[9],c_r[9];

	converttoPrimitive(rho_l,u_l,p_l,c_l,U1_l,U2_l,U3_l,U4_l);

	if (rho_l != rho_l ) std::cout << "ERROR:(in start of  hllc) rho_l is nan at timestep = " << timestep << " and x = " << x   << std::endl;
	if (u_l != u_l ) std::cout << "ERROR:(in start of  hllc) velocity_l is nan at timestep = " << timestep << " and x = " << x << " U2_l = " << U2_l << std::endl;
	if (p_l != p_l )std::cout << "ERROR:(in start of  hllc) pressure_l is nan at timestep = " << timestep << " and x = " << x   << std::endl; 
	if ( rho_l <= 0 )std::cout << "ERROR:(in start of  hllc) rho_l is negative or ZERO at timestep = " << timestep << " and x = " << x   << " U1_l = " << U1_l << std::endl; 
	if ( p_l < 0 ) std::cout << "ERROR:(in start of  hllc) pressure_l is negative at timestep = " << timestep << " and x = " << x  << " Pressure_l = " << p_l << std::endl;
 
	converttoPrimitive(rho_r,u_r,p_r,c_r,U1_r,U2_r,U3_r,U4_r);

	if (rho_r != rho_r ) std::cout << "ERROR:(in start of  hllc) rho_r is nan at timestep = " << timestep << " and x = " << x   << std::endl;
	if (u_r != u_r ) std::cout << "ERROR:(in start of  hllc) velocity_r is nan at timestep = " << timestep << " and x = " << x   << std::endl;
	if (p_r != p_r )std::cout << "ERROR:(in start of  hllc) pressure_r is nan at timestep = " << timestep << " and x = " << x << " U1_r = " << U1_l << " U2_r = " << U2_r << " U3_r = " << U3_r << std::endl; 
	if ( rho_r < 0 )	std::cout << "ERROR:(in start of  hllc) rho_r is negative at timestep = " << timestep << " and x = " << x   << std::endl; 
	if ( p_r < 0 ) std::cout << "ERROR:(in start of  hllc) pressure_r is negative at timestep = " << timestep << " and x = " << x  << " Pressure_r = " << p_r << std::endl;

	//std::cout << "rho_l = " << rho_l << " u_l = " << u_l << " p_l = "  << p_l << " c_l = " << c_l << " at timestep = " << timestep << " and at position " << x << std::endl;
	//std::cout << "rho_l = " << rho_l << " u_l = " << u_l << " p_l = "  << p_l << " c_l = " << c_l << " at timestep = " << timestep << " and at position " << x << std::endl;
  
	//COMPUTE LEFT AND RIGHT FLUXES

	double F1_l,F2_l,F3_l,F1_r,F2_r,F3_r;
	double F4_l[9],F4_r[9];
  
	calculateFluxes(U1_l,U2_l,U3_l,U4_l,F1_l,F2_l,F3_l,F4_l);
	calculateFluxes(U1_r,U2_r,U3_r,U4_r,F1_r,F2_r,F3_r,F4_r);

	//SPEED OF SOUND
 
	double a_l,a_r;
	a_l=soundSpeed(p_l,rho_l);
	a_r=soundSpeed(p_r,rho_r);

	if (a_l != a_l) std::cout << "ERROR: a_l is a NaN. at timestep = " << timestep << " p_l = " << p_l << " rho_l = " << rho_l <<  std::endl;
	if (a_r != a_r) std::cout << "ERROR: a_r is a NaN. at timestep = " << timestep << " p_r = " << p_r << " rho_r = " << rho_r <<  std::endl;
 
  
	//CALCULATION PRESSURE IN STAR REGION p*

	double p_pvrs, p_star=0;
	p_pvrs=0.5*(p_l+p_r)-0.5*(u_r-u_l)*0.5*(rho_l+rho_r)*0.5*(a_l+a_r);
	p_star = std::max(0.0,p_pvrs);
	if (p_pvrs != p_pvrs) std::cout << "ERROR: p_pvrs is a NaN. at timestep = " << timestep << " p_l = " << p_l << " rho_l = " << rho_l <<  std::endl;
	if (p_star < 0 ) std::cout << "ERROR: Pressure in star region < 0  at timestep = " << timestep << std::endl;
  
	//COMPUTE S-WAVES SPEEDS
  
	double S_l, S_r, S_star;
  
	double q_l=1,q_r=1;

	if ( p_star > p_l ) q_l=sqrt( 1 + (gamma+1)*(p_star/p_l - 1)/(2*gamma) );
	if ( p_star > p_r ) q_r=sqrt( 1 + (gamma+1)*(p_star/p_r - 1)/(2*gamma) );
	if (q_l != q_l) std::cout << "q_l is a NaN" << std::endl;
	if (q_r != q_r) std::cout << "q_r is a NaN" << std::endl;
  
	S_l=u_l-a_l*q_l;
	S_r=u_r+a_r*q_r;
  
  
	//ALTERNATIVE WAY OF APPROXIMATING SPEEDS
		/* DOESNOT WORK WITH THIS
		double S_max;
		S_max = std::max(fabs(u_l)+a_l,fabs(u_r)+a_r);
		S_l=-S_max;
		S_r=S_max;
		*/
	//if (timestep < 3) std::cout << "S_max = " << S_max  << " at timestep = " << timestep << " and at position " << x << std::endl;

	// std::cout << "u_l = " << u_l << " a_l = " << a_l << /*" q_l = "  << q_l << */" S_l = " << S_l << " at timestep = " << timestep << " and at position " << x << std::endl;
	// std::cout << "u_r = " << u_r << " a_r = " << a_r << /*" q_r = "  << q_r << */" S_r = " << S_r << " at timestep = " << timestep << " and at position " << x << std::endl;

	S_star= (p_r - p_l + rho_l*u_l*(S_l-u_l) - rho_r*u_r*(S_r-u_r)) / (rho_l*(S_l-u_l) - rho_r*(S_r-u_r));

	if ( S_l > S_r ) {
		std::cout << "ERROR! Speed of left wave faster than speed of right wave"  << std::endl;
		std::cout << "S_l = " << S_l << " S_star = " << S_star  << " S_r = " << S_r << " at timestep = " << timestep << " and at position " << x << std::endl;
	}

	//FLUXES

	if ( S_l >= 0 ){
		F1=F1_l;
		F2=F2_l;
		F3=F3_l;
		for (int i=0;i<9;i++) F4[i]=F4_l[i];
		//std::cout << "In region F_l" << std::endl;
	}
	else if (S_l<0 && S_star >= 0){
    
		double coeff,U1_lstar,U2_lstar,U3_lstar;
		double U4_lstar[9];
		
		coeff=rho_l*(S_l-u_l)/(S_l-S_star);

		U1_lstar=coeff;
		U2_lstar=coeff*S_star;
		U3_lstar=coeff*(U3_l/U1_l+(S_star - u_l)*( S_star + p_l/gamma/( rho_l*(S_l-u_l) ) ) );
		for (int i=0;i<9;i++) U4_lstar[i]=coeff*c_l[i];

		F1=F1_l + S_l*(U1_lstar - U1_l);
		F2=F2_l + S_l*(U2_lstar - U2_l);
		F3=F3_l + S_l*(U3_lstar - U3_l);
		for (int i=0;i<9;i++) F4[i]=F4_l[i] + S_l*(U4_lstar[i] - U4_l[i]);		
	}
	else if (S_star < 0 && S_r > 0){
		double coeff,U1_rstar,U2_rstar,U3_rstar;
	    double U4_rstar[9];
	    
	    coeff=rho_r*(S_r-u_r)/(S_r-S_star);

	    U1_rstar=coeff;
	    U2_rstar=coeff*S_star;
	    U3_rstar=coeff*(U3_r/U1_r+(S_star - u_r)*( S_star + p_r/gamma/( rho_r*(S_r-u_r) ) ) );
	    for (int i=0;i<9;i++) U4_rstar[i]=coeff*c_r[i];

	    F1=F1_r + S_r*(U1_rstar - U1_r);
	    F2=F2_r + S_r*(U2_rstar - U2_r);
	    F3=F3_r + S_r*(U3_rstar - U3_r);
	    for (int i=0;i<9;i++) F4[i]=F4_r[i] + S_r*(U4_rstar[i] - U4_r[i]);
	}
	else if (S_r <= 0 ){
		F1=F1_r;
		F2=F2_r;
		F3=F3_r;
		for (int i=0;i<9;i++) F4[i]=F4_r[i];
		//std::cout << "In region F_r" << std::endl;
	}
	else {
		std::cout << "ERROR! : Region of flux(0) not found! at timestep = " << timestep << " and at position " << x << std::endl;
		std::cout << "S_l = " << S_l << " S_star = " << S_star  << " S_r = " << S_r << std::endl;
	}
}
