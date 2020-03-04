#include "cavityVerification.h"
#include <iostream>
#include <fstream>



cavityVerification::cavityVerification()
{
}


cavityVerification::~cavityVerification()
{
}

void cavityVerification::setVariables(int Nx, int Ny)
{
	nx = Nx;
	ny = Ny;

	u_Numeric = new double[ny];
	u_Anlytic = new double[ny];
	y_sec = new double[ny];

	width = ny - 2;
	lenght = nx;

	Re = 10;
	nu = 1. / 6;
	DP = -Re * lenght*nu*nu / (width*width*width);
	rho0 = 1.0;
	cs = 1. / sqrt(3.);
	rhoIn = rho0 - 0.5*DP / (cs*cs);
	rhout = rho0 + 0.5*DP / (cs*cs);
}


void cavityVerification::extrctNumricProfile(double *dataField) {

	int nNodeProfil = 0;

	int i = nx/2;
	for (int j = 0; j < ny; j++) {
		int iNode = j * nx + i;
		y_sec[nNodeProfil] = nNodeProfil;
		u_Numeric[nNodeProfil] = dataField[iNode];
		nNodeProfil++;
	}

}


void cavityVerification::writeNumericVelcity(double *dataField) {

	extrctNumricProfile(dataField);

	FILE *file2;
	fopen_s(&file2, "u_Numeric.plt", "w");

	fprintf(file2, "u_Numeric \n ");
	for (int j = 0; j < ny; j++) {
		fprintf(file2, "%.3f  %.3f \n ", u_Numeric[j], y_sec[j]);
	}
	fclose(file2);

}








