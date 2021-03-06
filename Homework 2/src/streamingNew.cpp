#include "lattice.h"
#include "domain.h"
#include <iostream>

int streamTo(int iNode, int k, double cx, double cy, Domain domain) {

	int nx = domain.getNx();
	int ny = domain.getNx();

	int i = iNode % nx;
	int j = iNode / nx;


	int pop = (j + cy)*nx + i + cx;

	return pop;
}


void streamingNew(int pop, Lattice lat, int k, double fi) {
	lat.f_[pop][k] = fi;
}


