#include "lattice.h"
#include "domain.h"

void streaming(int iNode, double **Cxy, Domain domain, Lattice lat, int k, double fi) {

	int nx = domain.getNx();
	int ny = domain.getNy();

	int i = iNode % nx;
	int j = iNode / nx;


	double cx = Cxy[0][k];
	double cy = Cxy[1][k];

	int pop = (j + cy)*nx + i + cx;

	lat.f_[pop][k] = fi;
}
