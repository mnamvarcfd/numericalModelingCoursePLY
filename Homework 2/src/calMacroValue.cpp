#include "domain.h"
#include "LatticStencil.h"
#include "lattice.h"

void calMacroValue(Domain &domain, LatticStencil stencil, Lattice lat){
	
	double fi;
	double sumFi;
	double moment1x;
	double moment1y;
	double cx;
	double cy;
 
	int ntot = domain.getNTot();

	for (int j = 0; j < ntot; j++)
	{
		sumFi = 0.0;
		moment1x = 0.0;
		moment1y = 0.0;

		for (int k = 0; k < 9; k++) {
			cx = stencil.cx[k];
			cy = stencil.cy[k];

			fi = lat.f0_[j][k];

			sumFi += fi;
			moment1x += fi * cx;
			moment1y += fi * cy;
		}

		lat.rho_[j] = sumFi;

		lat.u_[j][0] = moment1x / sumFi;
		lat.u_[j][1] = moment1y / sumFi;
	}

}