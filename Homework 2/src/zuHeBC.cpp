#include "domain.h"
#include "lattice.h"

void zouHeBottom(Domain domain, Lattice lat, double rho) {

	int ny = domain.ny;

	int iNodeBot;
	for (int j = 0; j < ny; j++) {
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


		double u = 1.0 - (f0 + f3 + f7 + 2 * (f4 + f5 + f6)) / rho;

		lat.f0_[iNodeBot][1] = f5 + (2. / 3)*rho*u;
		lat.f0_[iNodeBot][2] = f6 - 0.5*(f3 - f7) + (1. / 6)*rho*u;
		lat.f0_[iNodeBot][8] = f4 + 0.5*(f3 - f7) + (1. / 6)*rho*u;
	}

}


void zouHeTop(Domain domain, Lattice lat, double rho) {

	int nx = domain.nx;
	int ny = domain.ny;

	int iNodeTop;
	for (int j = 0; j < ny; j++) {
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

		double u = (f0 + f3 + f7 + 2 * (f1 + f2 + f8)) / rho - 1.0;

		lat.f0_[iNodeTop][5] = f1 - (2. / 3)*rho*u;
		lat.f0_[iNodeTop][6] = f2 + 0.5*(f3 - f7) - (1. / 6)*rho*u;
		lat.f0_[iNodeTop][4] = f8 - 0.5*(f3 - f7) - (1. / 6)*rho*u;

	}

}



