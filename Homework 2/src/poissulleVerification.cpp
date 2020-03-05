#include <iostream>
#include <fstream>
#include "poissulleVerification.h"

poissulleVerification::~poissulleVerification()
{
}
poissulleVerification::poissulleVerification()
{

}

void poissulleVerification::setVariablesFullWay(int Nx, int Ny)
{
	nx = Nx;
	ny = Ny;

	width = ny-2;
	lenght = nx;

	frstNode = 1;
	lastNode = ny-2;

	u_Numeric = new double[ny];
	u_Anlytic = new double[ny];
	y_sec = new double[ny];

	for (int j = frstNode; j <= lastNode; j++) {
		y_sec[j] = ((double)j - width / 2. - 0.5);
	}

	Re = 10;
	nu = 1. / 6;
	DP = -Re * lenght*nu*nu / (width*width*width);
	rho0 = 1.0;
	cs = 1. / sqrt(3.);
	rhoIn = rho0 - 0.5*DP / (cs*cs);
	rhout = rho0 + 0.5*DP / (cs*cs);
}

void poissulleVerification::setVariablesHalfWay(int Nx, int Ny)
{
	nx = Nx;
	ny = Ny;

	width = ny;
	lenght = nx;

	frstNode = 0;
	lastNode = ny - 1;

	u_Numeric = new double[ny];
	u_Anlytic = new double[ny];
	y_sec = new double[ny];

	for (int j = frstNode; j <= lastNode; j++) {
		y_sec[j] = ((double)j - width / 2. + 0.5);
	}

	Re = 10;
	nu = 1. / 6;
	DP = -Re * lenght*nu*nu / (width*width*width);
	rho0 = 1.0;
	cs = 1. / sqrt(3.);
	rhoIn = rho0 - 0.5*DP / (cs*cs);
	rhout = rho0 + 0.5*DP / (cs*cs);
}


void poissulleVerification::analyticVelcity() {

	double scale = ((DP / (8 * nu*lenght)) * (width*width));

	for (int j = frstNode; j <= lastNode; j++) {
		double y = y_sec[j];
		u_Anlytic[j] = (DP / (2 * nu*lenght)) * (width*width / 4. - y * y);
		u_Anlytic[j] = u_Anlytic[j] / scale;
	}

}

void poissulleVerification::writeAnalyticVelcity() {

	analyticVelcity();

	FILE *file2;
	fopen_s(&file2, "u_Analytic.plt", "w");

	fprintf(file2, "u_Analytic \n ");
	for (int j = frstNode; j <= lastNode; j++) {
		fprintf(file2, "%.3f  %.3f \n ", u_Anlytic[j], y_sec[j] / width);
	}
	fclose(file2);
}

void poissulleVerification::extrctNumricProfile(double *dataField) {

	int ntot = nx * ny;

	int nNodeProfil = 0;
	for (int i = 0; i < ntot; i++) {
		if (i%nx == 2) {
			u_Numeric[nNodeProfil] = -dataField[i];
			nNodeProfil++;
		}
	}

	double uCenter = 0.000000000001;
	for (int i = frstNode; i <= lastNode; i++)
		if (u_Numeric[i] > uCenter) uCenter = u_Numeric[i];

	for (int i = frstNode; i <= lastNode; i++)
		u_Numeric[i] = u_Numeric[i] / uCenter;

}


void poissulleVerification::writeNumericVelcity(double *dataField) {

	extrctNumricProfile(dataField);

	FILE *file2;
	fopen_s(&file2, "u_Numeric.plt", "w");

	fprintf(file2, "u_Numeric \n ");
	for (int j = frstNode; j <= lastNode; j++) {
		fprintf(file2, "%.3f  %.3f \n ", u_Numeric[j], y_sec[j] / width);
	}
	fclose(file2);

}


void poissulleVerification::erreurL2(double *dataField) {

	extrctNumricProfile(dataField);
	analyticVelcity();

	double error = 0.;
	for (int j = frstNode; j <= lastNode; j++)
	{
		error = error + (u_Anlytic[j] - u_Numeric[j])*(u_Anlytic[j] - u_Numeric[j]);
	}
	double erreurL2 = sqrt(error / (lastNode - frstNode + 1) );
	std::cout << "Norm of the error :  " << erreurL2 << std::endl;

}







