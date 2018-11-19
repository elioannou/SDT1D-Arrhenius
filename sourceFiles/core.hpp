#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <ctime>
#include <string>
#include <sstream>

double soundSpeed(const double p, const double rho)
{
	//double gamma=1.4;
	double a;
	a=sqrt(p/rho);
	return a;
}

void converttoConserved(const double W1,const double W2,const double W3, const double W4[], double &U1, double &U2, double &U3, double U4[])
{
	const double gamma=1.4;
	U1=W1;
	U2=W1*W2;
	U3=0.5*W1*W2*W2+(W3/(gamma*(gamma-1)));
	for (int i=0;i<9;i++) U4[i]=W1*W4[i];
}

void converttoPrimitive(double &W1, double &W2, double &W3, double W4[], const double U1, const double U2, const double U3, const double U4[])
{
	const double gamma=1.4;
	W1=U1;
	W2=U2/W1;
	W3=(U3-0.5*U2*U2/U1)*gamma*(gamma-1);
	for (int i=0;i<9;i++) W4[i]=U4[i]/U1;
}

void calculateFluxes(const double U1,const double U2,const double U3, const double U4[], double &F1, double &F2, double &F3, double F4[]) //conserved variables ONLY!
{
	const double gamma=1.4;
	double rho,u,p;
	double c[9];
	converttoPrimitive(rho,u,p,c,U1,U2,U3,U4);
	F1=U2;
	F2=U2*u+p/gamma;
	F3=u*(U3+p/gamma);
	for (int i=0;i<9;i++) F4[i]=u*U4[i];
}

void print4Variables(int M, double dx, double V1[], double V2[], double V3[], double** V4)
{

	for (int i=0;i<M+4;i++) {
			std::cout  << dx*(i-1) << "\t" << V1[i] << "\t" << V2[i] << "\t" << V3[i] << "\t";
			for (int j=0;j<9;j++){
				std::cout << V4[i][j] << "\t";
			}
			std::cout << std::endl;
		}
}

void exportResults(std::string simulation, int M, int timestep, double simtime, double dx, double V1[], double V2[], double V3[], double** V4, double V5[], double V6[],clock_t runtime)
{
	std::string timestamp;			// string which will contain the result
	std::ostringstream stream;   	// stream used for the conversion

	stream << simtime*100;      	// insert the textual representation of 'Number' in the characters in the stream
	timestamp = stream.str(); 		// set 'Result' to the contents of the stream

	std::ofstream periodic_data;

    periodic_data.open(("./output/"+simulation+"/snapshot_"+timestamp+".dat").c_str());
	if (periodic_data.is_open()) {
		periodic_data << "|||||||||||||||||||||||||||||||||***SIMULATION: " << simulation << "***||||||||||||||||||||||||||||||||||" << std::endl;
		periodic_data << "Sanpshot time = " << simtime << ". Timestep = " << timestep << ". Run time = " << runtime/CLOCKS_PER_SEC << " seconds" << std::endl;
		periodic_data << "No. of cells = " << M << std::endl;
		periodic_data << "#Position\t Density\t Velocity\t Pressure\t Specific_Volume I_Energy\t Temperature\t Species" << std::endl;
		// STARTING from 6 to acount for wall heating.
		for (int i=0;i<M+2;i++) {
			periodic_data << dx*(i-1) << "\t" << V1[i] << "\t" << V2[i] << "\t" << V3[i] << "\t" << 1/V1[i] << "\t" << V5[i] << "\t"  << V6[i] << "\t";
			for (int j=0;j<9;j++){
				periodic_data << V4[i][j] << "\t";
			}
			periodic_data << std::endl;
		}
	}

	else {
		std::cout << "ERROR! Unable to open file to output data at time " << simtime << std::endl;
		exit(-1);
	}
  periodic_data.close();
}
