#include "domain.h"
#include "lattice.h"

void periodicBC(Domain domain, Lattice lat) {
	
	int nx = domain.getNx();
	int ny = domain.getNy();

	int iNodeTop;
	int iNodeBot;
	for (int j = 1; j < nx-1 ; j++) {

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
