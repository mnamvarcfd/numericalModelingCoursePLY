#include "domain.h"
#include "lattice.h"

void periodic(Domain domain, Lattice lat) {
	
	int nx = domain.getNx();
	int ny = domain.getNy();

	int iNodeTop;
	int iNodeBot;
	for (int j = 0; j < nx ; j++) {

		iNodeTop = (ny - 1) * nx + j;
		iNodeBot = j;

		lat.f_[iNodeBot][2] = lat.f_[iNodeTop][2];
		lat.f_[iNodeBot][6] = lat.f_[iNodeTop][6];
		lat.f_[iNodeBot][5] = lat.f_[iNodeTop][5];

		lat.f_[iNodeTop][4] = lat.f_[iNodeBot][4];
		lat.f_[iNodeTop][7] = lat.f_[iNodeBot][7];
		lat.f_[iNodeTop][8] = lat.f_[iNodeBot][8];            

	}

}



void zouHeBottom(Domain domain, Lattice lat, double rho) {

	int nx = domain.getNx();

	int iNodeBot;
	for (int j = 1; j < nx - 1; j++) {
		iNodeBot = j;

		double f0 = lat.f0_[iNodeBot][0];
		double f1 = lat.f0_[iNodeBot][1];
		double f2 = lat.f0_[iNodeBot][2];
		double f3 = lat.f0_[iNodeBot][3];
		double f4 = lat.f0_[iNodeBot][4];
		double f5 = lat.f0_[iNodeBot][5];
		double f6 = lat.f0_[iNodeBot][6];
		double f7 = lat.f0_[iNodeBot][7];
		double f8 = lat.f0_[iNodeBot][8];


		double u = 1.0 - (f0 + f1 + f3 + 2 * (f8 + f4 + f7)) / rho;

		lat.f_[iNodeBot][2] = f4 + (2. / 3)*rho*u;
		lat.f_[iNodeBot][5] = f7 - 0.5*(f1 - f3) + (1. / 6)*rho*u;
		lat.f_[iNodeBot][6] = f8 + 0.5*(f1 - f3) + (1. / 6)*rho*u;
	}

}


void zouHeTop(Domain domain, Lattice lat, double rho) {

	int nx = domain.getNx();
	int ny = domain.getNy();

	int iNodeTop;
	for (int j = 1; j < nx - 1; j++) {
		iNodeTop = (ny - 1) * nx + j;

		double f0 = lat.f0_[iNodeTop][0];
		double f1 = lat.f0_[iNodeTop][1];
		double f2 = lat.f0_[iNodeTop][2];
		double f3 = lat.f0_[iNodeTop][3];
		double f4 = lat.f0_[iNodeTop][4];
		double f5 = lat.f0_[iNodeTop][5];
		double f6 = lat.f0_[iNodeTop][6];
		double f7 = lat.f0_[iNodeTop][7];
		double f8 = lat.f0_[iNodeTop][8];

		double u = (f0 + f1 + f3 + 2 * (f2 + f5 + f6)) / rho - 1.0;

		lat.f_[iNodeTop][4] = f2 - (2. / 3)*rho*u;
		lat.f_[iNodeTop][7] = f5 + 0.5*(f1 - f3) - (1. / 6)*rho*u;
		lat.f_[iNodeTop][8] = f6 - 0.5*(f1 - f3) - (1. / 6)*rho*u;

	}

}



