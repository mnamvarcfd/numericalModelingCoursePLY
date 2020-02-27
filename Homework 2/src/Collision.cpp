#include "domain.h"
#include "lattice.h"

void collision(Domain domain, Lattice lat, double tau, double dt, double *omega, double **xi,
	double c1, double c2, double c3, double c4, double c5, double *g) {
	
	int nx = domain.getNx();
	int ny = domain.getNy();

	double cx;
	double cy;
	double ro, ux, uy, udotu;

	int iNode;
	int iPullNode;

	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {

			iNode = j * nx + i;

			ro = lat.rho_[iNode];
			ux = lat.u_[iNode][0];
			uy = lat.u_[iNode][1];
			udotu = ux * ux + uy * uy;

			for (int k = 0; k < 9; k++) {
				cx = xi[0][k];
				cy = xi[1][k];

				double uiei = cx * ux + cy * uy;
				double uiei2 = uiei * uiei;
				double feq = ro * omega[k] * (1 + c1 * uiei + c2 * uiei2 + c3 * udotu);

				double si = c5 * omega[k] * ((cy - uy)*c1 + (cx * ux + cy * uy)*cy*c4) * g[1];

				double fi = lat.f0_[iNode][k];
				lat.f0_[iNode][k] = fi - (fi - feq) / tau + si * dt;
			}

		}
	}

}

