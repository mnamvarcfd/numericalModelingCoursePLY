#include "domain.h"
#include "lattice.h"

void calMacroValue(Domain &domain, double **xi, Lattice lat, double *g, double dt){
	
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
			cx = xi[0][k];
			cy = xi[1][k];

			fi = lat.f_[j][k];

			sumFi += fi;
			moment1x += fi * cx;
			moment1y += fi * cy;
		}

		lat.rho_[j] = sumFi;

		lat.u_[j][0] = (moment1x + 0.5*g[0] * dt) / sumFi; //x
		lat.u_[j][1] = (moment1y + 0.5*g[1] * dt) / sumFi; //y
	}

}