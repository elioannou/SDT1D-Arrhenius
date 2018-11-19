


double delta=1e-7;

void slopeLimiter(double D1[],double D2[],double D3[],double** D4, double U1[],double U2[],double U3[],double** U4, int size, int timestep,int option)
{
	int i;
	double D1_minus=0,D1_plus=0,D2_minus=0,D2_plus=0,D3_minus=0,D3_plus=0,R1=0,R2=0,R3=0;
	double D4_minus[9],D4_plus[9],R4[9];
	double X1=0,X2=0,X3=0;
	double X4[9];
	for (int j=0;j<9;j++) {
		D4_minus[j]=D4_plus[j]=R4[j]=X4[j]=0;
	}

	for (i=1;i<(size+3);i++){
    
		D1_minus=U1[i]-U1[i-1];
		D1_plus=U1[i+1]-U1[i];

		D2_minus=U2[i]-U2[i-1];
		D2_plus=U2[i+1]-U2[i];   
	 
		D3_minus=U3[i]-U3[i-1];
		D3_plus=U3[i+1]-U3[i];

		for (int j=0;j<9;j++) {
			D4_minus[j]=U4[i][j]-U4[i-1][j];
			D4_plus[j]=U4[i+1][j]-U4[i][j];
		}    
		
		if ( fabs(D1_minus) < delta && fabs(D1_plus) < delta ) R1=1;
		else if (fabs(D1_minus) > delta && fabs(D1_plus) < delta) {
			if (D1_plus < 0) R1=-D1_minus/delta;
			else R1=D1_minus/delta;
		}
		else R1=D1_minus/D1_plus;

		if ( fabs(D2_minus) < delta && fabs(D2_plus) < delta ) R2 = 1;
		else if ( fabs(D2_minus) > delta && fabs(D2_plus) < delta ) {
		if (D2_plus < 0) R2 = -D2_minus/delta;
      else R2 = D2_minus/delta;
    }
    else R2=D2_minus/D2_plus;

    if ( fabs(D3_minus) < delta && fabs(D3_plus) < delta ) R3 = 1;
    else if ( fabs(D3_minus) > delta && fabs(D3_plus) < delta ) {
      if (D3_plus < 0) R3 = -D3_minus/delta;
      else R3 = D3_minus/delta;
    }
    else R3=D3_minus/D3_plus;

	for (int j=0;j<9;j++) {
	    if ( fabs(D4_minus[j]) < delta && fabs(D4_plus[j]) < delta ) R4[j] = 1;
	    else if ( fabs(D4_minus[j]) > delta && fabs(D4_plus[j]) < delta ){
	      if (D4_plus[j] <0) R4[j] = -D4_minus[j]/delta;
	      else R4[j] = D4_minus[j]/delta;
	    }
	    else R4[j]=D4_minus[j]/D4_plus[j];
	}
    //MINBEE
    if (option == 0){

		if ( R1 <= 0) X1=0;
		else if ( R1 > 0 && R1 <= 1 ) X1=R1;
		else if ( R1 > 1 ) {
			X1=2/(1+R1);
			if ( X1 > 1 ) X1=1;  
		}
		else std::cout << "X1 not assigned" << " at timestep = " << timestep << " i = " << i << " D1_minus = " <<  D1_minus << " D1_plus = " << D1_plus << std::endl;
 
		if ( R2 <= 0) X2=0;
		else if ( R2 > 0 && R2 <= 1 ) X2=R2;
		else if ( R2 > 1 ) {
			X2=2/(1+R2);
			if ( X2 > 1 ) X2=1;  
		}
		else std::cout << "X2 not assigned" << " at timestep = " << timestep << " i = " << i << " D2_minus = " <<  D2_minus << " D2_plus = " << D2_plus << std::endl;

		if ( R3 <= 0) X3=0;
		else if ( R3 > 0 && R3 <= 1 ) X3=R3;
		else if ( R3 > 1 ) {
			X3=2/(1+R3);
			if ( X3 > 1 ) X3=1;  
		}
		else std::cout << "X3 not assigned" << " at timestep = " << timestep << " i = " << i << " D3_minus = " <<  D3_minus << " D3_plus = " << D3_plus << std::endl;

		for (int j=0;j<9;j++) {
			if ( R4[j] <= 0) X4[j]=0;
			else if ( R4[j] > 0 && R4[j] <= 1 ) X4[j]=R4[j];
			else if ( R4[j] > 1 ) {
				X4[j]=2/(1+R4[j]);
				if ( X4[j] > 1 ) X4[j]=1;  
			}
			else std::cout << "X4 not assigned" << " at timestep = " << timestep << " i = " << i << " D4_minus = " <<  D4_minus[j] << " D4_plus = " << D4_plus[j] << std::endl;
		}

    }
    //SUPERBEE
    if (option == 1){
    
      if ( R1 <= 0) X1=0;
      else if ( (R1 > 0) && (R1 < 0.5) ) X1=2*R1;
      else if ( (R1 >= 0.5) && (R1 <= 1) ) X1=1;
      else if ( R1 > 1 ) {
	X1=2/(1+R1);
	if ( X1 > R1 ) X1=R1;
	if ( X1 > 2) X1=2;
      }
      else std::cout << "X1 not assigned" << " at timestep = " << timestep << " i = " << i << " D1_minus = " <<  D1_minus << " D1_plus = " << D1_plus  << std::endl;

      if ( R2 <= 0) X2=0;
      else if ( (R2 > 0) && (R2 < 0.5) ) X2=2*R2;
      else if ( (R2 >= 0.5) && (R2 <= 1) ) X2=1;
      else if ( R2 > 1 ) {
	X2=2/(1+R2);
	if ( X2 > R2 ) X2=R2;
	if (X2 > 2) X2=2;
      }
      else std::cout << "X2 not assigned" << " at timestep = " << timestep << " i = " << i << " D2_minus = " <<  D2_minus << " D2_plus = " << D2_plus << std::endl;

      if ( R3 <= 0) X3=0;
      else if ( R3 > 0 && R3 < 0.5 ) X3=2*R3;
      else if ( (R3 >= 0.5) && (R3 <= 1) ) X3=1;
      else if ( R3 > 1 ){
			X3=2/(1+R3);
			if ( X3 > R3 ) X3=R3;
			if (X3 > 2) X3=2;
		}
      else std::cout << "X3 not assigned" << " at timestep = " << timestep << " i = " << i << " D3_minus = " <<  D3_minus << " D3_plus = " << D3_plus << std::endl;

	for (int j=0;j<9;j++) {
		if ( R4[j] <= 0) X4[j]=0;
		else if ( R4[j] > 0 && R4[j] < 0.5 ) X4[j]=2*R4[j];
		else if ( (R4[j] >= 0.5) && (R4[j] <= 1) ) X4[j]=1;
		else if ( R4[j] > 1 ) {
			X4[j]=2/(1+R4[j]);
			if ( X4[j] > R4[j] ) X4[j]=R4[j];
			if (X4[j] > 2) X4[j]=2;    
		}
      else std::cout << "X4 not assigned" << " at timestep = " << timestep << " i = " << i << " R4 = " << R4[j] << " D4_minus = " <<  D4_minus[j] << " D4_plus = " << D4_plus[j] << std::endl;   
    }

    D1[i]=X1*(0.5*D1_minus+0.5*D1_plus);
    D2[i]=X2*(0.5*D2_minus+0.5*D2_plus);
    D3[i]=X3*(0.5*D3_minus+0.5*D3_plus);
    for (int j=0;j<9;j++) D4[i][j]=X4[j]*(0.5*D4_minus[j]+0.5*D4_plus[j]);
    
	}
	}
}
