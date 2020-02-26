#include "domain.h"
#include "lattice.h"

void collision(int iNode, Lattice lat, double tau, double dt,double ux, double uy, double *omega, double **xi,
	double c1, double c2, double c3, double c4, double c5, double *g, double udotu, double ro, double *fiOut) {
	
	double cx;
	double cy;

	for (int k = 0; k < 9; k++) {
		cx = xi[0][k];
		cy = xi[1][k];

		double uiei = cx * ux + cy * uy;
		double uiei2 = uiei * uiei;
		double feq = ro * omega[k] * (1 + c1 * uiei + c2 * uiei2 + c3 * udotu);

		double si = c5 * omega[k] * ((cy - uy)*c1 + (cx * ux + cy * uy)*cy*c4) * g[1];

		double fi = lat.f0_[iNode][k];
		fiOut[k] = fi - (fi - feq) / tau + si * dt;
	}
}

