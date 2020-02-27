#include "domain.h"
#include "lattice.h"

void bouncBackBC(Domain domain, Lattice lat) {
	int nx = domain.getNx();
	int ny = domain.getNy();

	int iNodeLeft;
	int iNodeRigt;
	for (int j = 0; j < ny; j++) {
		iNodeLeft = j * nx;
		iNodeRigt = (j + 1) * nx - 1;

		lat.f0_[iNodeLeft][5] = lat.f_[iNodeLeft][7];
		lat.f0_[iNodeLeft][1] = lat.f_[iNodeLeft][3];
		lat.f0_[iNodeLeft][8] = lat.f_[iNodeLeft][6];

		lat.f0_[iNodeRigt][7] = lat.f_[iNodeRigt][5];
		lat.f0_[iNodeRigt][3] = lat.f_[iNodeRigt][1];
		lat.f0_[iNodeRigt][6] = lat.f_[iNodeRigt][8];
	}

}