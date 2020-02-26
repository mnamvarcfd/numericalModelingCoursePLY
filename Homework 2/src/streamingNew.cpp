#include "lattice.h"
#include "domain.h"
#include <iostream>

int streamTo(int iNode, int k, double cx, double cy, Domain domain) {

	int nx = domain.getNx();
	int ny = domain.getNy();

	int i = iNode % nx;
	int j = iNode / nx;


	if (cx < 0.)cx = -1;
	if (cx > 0.)cx = 1;

	if (cy < 0.)cy = -1;
	if (cy > 0.)cy = 1;

	//std::cout << cx << std::endl;

	int pop = (j + cy)*nx + i + cx;

	return pop;
}


int bottomStreamTo(int iNode, int k, double cx, double cy, Domain domain) {

	int nx = domain.getNx();
	int ny = domain.getNy();

	int i = iNode % nx;
	int j = iNode / nx;


	int pop = (j + cy)*nx + i + cx;

	if (k == 7 || k == 4 || k == 8) pop = (ny - 1) * nx + i;

	return pop;
}

int topStreamTo(int iNode, int k, double cx, double cy, Domain domain) {

	int nx = domain.getNx();
	int ny = domain.getNy();

	int i = iNode % nx;
	int j = iNode / nx;


	int pop = (j + cy)*nx + i + cx;

	if (k == 6 || k == 2 || k == 5) pop = i;
	return pop;
}

void streamingNew(int pop, Lattice lat, int k, double fi) {
	lat.f_[pop][k] = fi;
}


