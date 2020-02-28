#include "lattice.h"
#include "domain.h"

void streaming(Domain domain, Lattice lat) {


	int nx = domain.getNx();
	int ny = domain.getNy();


	int iNode;
	int iPullNode;

	for (int j = 0; j < ny; j++) {
		for (int i = nx - 1; i > 0; i--) {

			iNode = j * nx + i;
			iPullNode = j * nx + i - 1;
			lat.f0_[iNode][1] = lat.f0_[iPullNode][1];
		}

		for (int i = 0; i < nx - 1; i++) {
			iNode = j * nx + i;
			iPullNode = j * nx + i + 1;
			lat.f0_[iNode][3] = lat.f0_[iPullNode][3];
		}
	}

	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny - 1; j++) {
			iNode = j * nx + i;
			iPullNode = (j + 1) * nx + i;
			lat.f0_[iNode][4] = lat.f0_[iPullNode][4];
		}

		for (int j = ny - 1; j > 0; j--) {
			iNode = j * nx + i;
			iPullNode = (j - 1) * nx + i;
			lat.f0_[iNode][2] = lat.f0_[iPullNode][2];
		}
	}

	for (int j = ny - 1; j > 0; j--) {
		for (int i = 0; i < nx - 1; i++) {
			iNode = j * nx + i;
			iPullNode = (j - 1) * nx + i + 1;
			lat.f0_[iNode][6] = lat.f0_[iPullNode][6];
		}

		for (int i = nx - 1; i > 0; i--) {

			iNode = j * nx + i;
			iPullNode = (j - 1) * nx + i - 1;
			lat.f0_[iNode][5] = lat.f0_[iPullNode][5];
		}
	}

	for (int j = 0; j < ny - 1; j++) {
		for (int i = 0; i < nx - 1; i++) {
			iNode = j * nx + i;
			iPullNode = (j + 1) * nx + i + 1;
			lat.f0_[iNode][7] = lat.f0_[iPullNode][7];
		}

		for (int i = nx - 1; i > 0; i--) {
			iNode = j * nx + i;
			iPullNode = (j + 1) * nx + i - 1;
			lat.f0_[iNode][8] = lat.f0_[iPullNode][8];
		}
	}
}
