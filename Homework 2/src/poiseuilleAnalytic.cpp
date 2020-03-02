#include "pch.h"
#include "Domain.h"
#include "lattice.h"

#include <stdio.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdio.h>


using namespace std;

void poiseuilleAnalytic(Domain mesh, double width, double DP, double lenght, double nu) {

	FILE *file2;
	fopen_s(&file2, "u_Analytic.plt", "w");

	double *u_Analytic = new double[mesh.ny];
	double scale = ((DP / (8 * nu*lenght)) * (width*width));

	for (int j = 0; j < mesh.ny; j++) {
		double y = ( (double)j - width/2. + 0.5);
		u_Analytic[j] = (DP / (2 * nu*lenght)) * (width*width/4. - y * y);

		u_Analytic[j] = u_Analytic[j] /*/ scale*/;


		fprintf(file2, "%.3f  %.3f \n ", u_Analytic[j] , y/ width);

	}

	fclose(file2);



}




void poiseuilleError(Domain mesh, double width, double DP, double lenght, double nu, Lattice lat) {

	int nx = mesh.nx;
	int ny = mesh.ny;
	int ntot = mesh.ntot;

	double *ux_Numeric = new double[nx];


	int cnt = 0;
	for (int i = 0; i < ntot; i++) {
		if (nx <= i && i < 2 * nx) {
			ux_Numeric[cnt] = -lat.u_[i][1];
			cnt++;
		}
	}

	double uCenter = 0.000000000001;
	for (int i = 0; i < cnt; i++) {
		if (ux_Numeric[i] > uCenter) uCenter = ux_Numeric[i];
	}

	for (int i = 0; i < cnt; i++) {
		double y = ((double)i - width / 2. + 0.5);
		ux_Numeric[i] = ux_Numeric[i] /*/ uCenter*/;
	}




	double *u_Analytic = new double[mesh.ny];
	double scale = ((DP / (8 * nu*lenght)) * (width*width));

	for (int j = 0; j < mesh.ny; j++) {
		double y = ((double)j - width / 2. + 0.5);
		u_Analytic[j] = (DP / (2 * nu*lenght)) * (width*width / 4. - y * y);

		u_Analytic[j] = u_Analytic[j] /*/ scale*/;
	}

	

	double erreur = 0.;
	for (int j = 0; j < nx; j++)
	{
		erreur = erreur + (u_Analytic[j] - ux_Numeric[j])*(u_Analytic[j] - ux_Numeric[j]);
	}
	double erreurL2 = sqrt(1. / nx*erreur);
	std::cout << "Norm of the error :  " << erreurL2 << std::endl;



}

