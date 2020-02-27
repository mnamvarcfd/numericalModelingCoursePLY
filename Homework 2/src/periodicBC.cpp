#include "domain.h"
#include "lattice.h"

void periodicBC(Domain domain, Lattice lat) {
	
	int nx = domain.getNx();
	int ny = domain.getNy();

	int iNodeTop;
	int iNodeBot;
	for (int j = 0; j < nx; j++) {

		iNodeTop = (ny - 1) * nx + j;
		iNodeBot = j;

		lat.f0_[iNodeBot][2] = lat.f0_[iNodeTop][2];
		lat.f0_[iNodeBot][6] = lat.f0_[iNodeTop][6];
		lat.f0_[iNodeBot][5] = lat.f0_[iNodeTop][5];

		lat.f0_[iNodeTop][4] = lat.f0_[iNodeBot][4];
		lat.f0_[iNodeTop][7] = lat.f0_[iNodeBot][7];
		lat.f0_[iNodeTop][8] = lat.f0_[iNodeBot][8];            
	}

}
