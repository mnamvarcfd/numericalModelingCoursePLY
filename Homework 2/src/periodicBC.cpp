#include "domain.h"
#include "lattice.h"

void periodicBC(Domain domain, Lattice lat) {
	
	int nx = domain.nx;
	int ny = domain.ny;

	int iNodeLeft;
	int iNodeRigt;
	for (int j = 1; j < ny-1; j++) {
		iNodeLeft = j * nx;
		iNodeRigt = (j + 1) * nx - 1;

		lat.f0_[iNodeRigt][8] = lat.f0_[iNodeLeft][8];
		lat.f0_[iNodeRigt][7] = lat.f0_[iNodeLeft][7];
		lat.f0_[iNodeRigt][6] = lat.f0_[iNodeLeft][6];

		lat.f0_[iNodeLeft][2] = lat.f0_[iNodeRigt][2];
		lat.f0_[iNodeLeft][3] = lat.f0_[iNodeRigt][3];
		lat.f0_[iNodeLeft][4] = lat.f0_[iNodeRigt][4];
	}


}
