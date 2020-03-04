#pragma once
class cavityVerification
{
public:
	int nx;
	int ny;
	double width;
	double lenght;
	double Re;
	double nu;
	double DP;
	double rho0;
	double cs;
	double rhoIn;
	double rhout;

	double *u_Numeric;
	double *u_Anlytic;
	double *y_sec;

	cavityVerification();
	~cavityVerification();
	void setVariables(int Nx, int Ny);
	void extrctNumricProfile(double * dataField);
	void writeNumericVelcity(double * dataField);
};

