#include "pch.h"
#include <iostream>
#include <fstream>
#include "poissulleVerificationFW.h"

poissulleVerificationFW::~poissulleVerificationFW()
{
}
poissulleVerificationFW::poissulleVerificationFW()
{
}

void poissulleVerificationFW::setVariables(int Nx, int Ny)
{
	nx = Nx;
	ny = Ny;

	u_Numeric = new double[ny];
	u_Anlytic = new double[ny];
	y_sec = new double[ny];

	width = ny-2;
	lenght = nx;

	Re = 10;
	nu = 1. / 6;
	DP = -Re * lenght*nu*nu / (width*width*width);
	rho0 = 1.0;
	cs = 1. / sqrt(3.);
	rhoIn = rho0 - 0.5*DP / (cs*cs);
	rhout = rho0 + 0.5*DP / (cs*cs);
}

void poissulleVerificationFW::analyticVelcity() {

	double scale = ((DP / (8 * nu*lenght)) * (width*width));

	for (int j = 1; j < ny - 1; j++) {
		double y = ((double)j - width / 2. - 0.5);
		u_Anlytic[j] = (DP / (2 * nu*lenght)) * (width*width / 4. - y * y);

		u_Anlytic[j] = u_Anlytic[j] / scale;
		y_sec[j] = y / width;
	}

}

void poissulleVerificationFW::writeAnalyticVelcity() {

	analyticVelcity();

	FILE *file2;
	fopen_s(&file2, "u_Analytic.plt", "w");

	fprintf(file2, "u_Analytic \n ");
	for (int j = 1; j < ny - 1; j++) {
		fprintf(file2, "%.3f  %.3f \n ", u_Anlytic[j], y_sec[j]);
	}
	fclose(file2);
}

void poissulleVerificationFW::extrctNumricProfile(double *dataField) {

	int ntot = nx * ny;

	int nNodeProfil = 0;
	for (int i = 0; i < ntot; i++) {
		if (i%nx == 2) {
			u_Numeric[nNodeProfil] = dataField[i];
			nNodeProfil++;
		}
	}

	double uCenter = 0.000000000001;
	for (int i = 1; i < nNodeProfil - 1; i++)
		if (u_Numeric[i] > uCenter) uCenter = u_Numeric[i];

	for (int i = 1; i < nNodeProfil - 1; i++)
		u_Numeric[i] = u_Numeric[i] / uCenter;

}


void poissulleVerificationFW::writeNumericVelcity(double *dataField) {

	extrctNumricProfile(dataField);

	FILE *file2;
	fopen_s(&file2, "u_Numeric.plt", "w");

	fprintf(file2, "u_Numeric \n ");
	for (int j = 1; j < ny - 1; j++) {
		fprintf(file2, "%.3f  %.3f \n ", u_Numeric[j], y_sec[j]);
	}
	fclose(file2);

}


void poissulleVerificationFW::erreurL2(double *dataField) {

	extrctNumricProfile(dataField);
	analyticVelcity();

	double error = 0.;
	for (int j = 1; j < ny - 1; j++)
	{
		error = error + (u_Anlytic[j] - u_Numeric[j])*(u_Anlytic[j] - u_Numeric[j]);
	}
	double erreurL2 = sqrt(error / (ny - 2));
	std::cout << "Norm of the error :  " << erreurL2 << std::endl;

}







