#pragma once
class poissulleVerificationFW
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

		~poissulleVerificationFW();

		poissulleVerificationFW();

		void setVariables(int Nx, int Ny);

		void analyticVelcity();

		void writeAnalyticVelcity();

		void extrctNumricProfile(double * dataField);

		void writeNumericVelcity(double * dataField);

		void erreurL2(double * dataField);

};

