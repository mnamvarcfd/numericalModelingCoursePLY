#include "domain.h"
#include "lattice.h"

void bouncBackBC(Domain domain, Lattice lat) {
	int nx = domain.nx;
	int ny = domain.ny;

	int iNodeTop;
	int iNodeBot;
	for (int j = 0; j < nx; j++) {
		iNodeTop = (ny - 1) * nx + j;
		iNodeBot = j;

		lat.f0_[iNodeTop][6] = lat.f_[iNodeTop][2];
		lat.f0_[iNodeTop][5] = lat.f_[iNodeTop][1];
		lat.f0_[iNodeTop][4] = lat.f_[iNodeTop][8];

		lat.f0_[iNodeBot][2] = lat.f_[iNodeBot][6];
		lat.f0_[iNodeBot][1] = lat.f_[iNodeBot][5];
		lat.f0_[iNodeBot][8] = lat.f_[iNodeBot][4];
	}

}