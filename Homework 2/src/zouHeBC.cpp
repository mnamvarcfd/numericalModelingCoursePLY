#include "domain.h"
#include "lattice.h"

void zouHeBottom(Domain domain, Lattice lat, int iNode, double rho, double *fiOut) {

		double f0 = lat.f0_[iNode][0];
		double f1 = lat.f0_[iNode][1];
		double f2 = lat.f0_[iNode][2];
		double f3 = lat.f0_[iNode][3];
		double f4 = lat.f0_[iNode][4];
		double f5 = lat.f0_[iNode][5];
		double f6 = lat.f0_[iNode][6];
		double f7 = lat.f0_[iNode][7];
		double f8 = lat.f0_[iNode][8];


		double u = 1.0 - (f0 + f1 + f3 + 2 * (f8 + f4 + f7)) / rho;

		fiOut[2] = f4 + (2. / 3)*rho*u;
		fiOut[5] = f7 - 0.5*(f1 - f3) + (1. / 6)*rho*u;
		fiOut[6] = f8 + 0.5*(f1 - f3) + (1. / 6)*rho*u;

}



void zouHeTop(Domain domain, Lattice lat, int iNode, double rho, double *fiOut) {

		double f0 = lat.f0_[iNode][0];
		double f1 = lat.f0_[iNode][1];
		double f2 = lat.f0_[iNode][2];
		double f3 = lat.f0_[iNode][3];
		double f4 = lat.f0_[iNode][4];
		double f5 = lat.f0_[iNode][5];
		double f6 = lat.f0_[iNode][6];
		double f7 = lat.f0_[iNode][7];
		double f8 = lat.f0_[iNode][8];

		double u = (f0 + f1 + f3 + 2 * (f2 + f5 + f6)) / rho - 1.0;

		fiOut[4] = f2 - (2. / 3)*rho*u;
		fiOut[7] = f5 + 0.5*(f1 - f3) - (1. / 6)*rho*u;
		fiOut[8] = f6 - 0.5*(f1 - f3) - (1. / 6)*rho*u;

}