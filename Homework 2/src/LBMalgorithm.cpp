#include "LBMalgorithm.h"

LBMalgorithm::~LBMalgorithm()
{
}


LBMalgorithm::LBMalgorithm(Grid &grid, double tauIn, double cs, LatticStencil &stencil_)
{
	stencil = stencil_;
	nNode = grid.getNNode();
	nx = grid.getNx();
	ny = grid.getNy();
	c1 = 1. / (cs*cs);
	c2 = 1. / (2.*cs*cs*cs*cs);
	c3 = -1. / (2.*cs*cs);
	tau = tauIn;

	R = grid.R;
	U = grid.U;
	V = grid.V;
	fVal = grid.fVal;
}

int LBMalgorithm::getIndex(int i, int j) {
	return j * nx + i;
}

void LBMalgorithm::init() {

	for (int i = 0; i < nNode; i++) {
		U[i] = 0.0;
		V[i] = 0.0;
		R[i] = 1.0;

		for (int k = 0; k < 9; k++)
			fVal[i][k] = stencil.wt[k];
	}
}


void LBMalgorithm::calMacroValue() {

	double fi;
	double sumFi;
	double moment1x;
	double moment1y;
	double cx;
	double cy;

	for (int j = 0; j < nNode; j++)
	{
		sumFi = 0.0;
		moment1x = 0.0;
		moment1y = 0.0;

		for (int k = 0; k < 9; k++) {
			cx = stencil.cx[k];
			cy = stencil.cy[k];

			fi = fVal[j][k];

			sumFi += fi;
			moment1x += fi * cx;
			moment1y += fi * cy;
		}

		R[j] = sumFi;
		U[j] = moment1x / sumFi;
		V[j] = moment1y / sumFi;
	}

}


void LBMalgorithm::collision() {

	double uiei, uiei2, feq, fi, cx, cy, w, ro, ux, uy, udotu;
	int iNode;

	for (int j = 1; j < ny-1; j++) {
		for (int i = 0; i < nx; i++) {

			iNode = getIndex(i, j);

			ro = R[iNode];
			ux = U[iNode];
			uy = V[iNode];
			udotu = ux * ux + uy * uy;

			for (int k = 0; k < 9; k++) {
				cx = stencil.cx[k];
				cy = stencil.cy[k];
				w = stencil.wt[k];

				uiei = cx * ux + cy * uy;
				uiei2 = uiei * uiei;
				feq = ro * w * (1 + c1 * uiei + c2 * uiei2 + c3 * udotu);

				fi = fVal[iNode][k];
				fVal[iNode][k] = fi - (fi - feq) / tau;
			}

		}
	}

}


void LBMalgorithm::streaming() {

	int iNode;
	int iPullNode;

	for (int j = 0; j < ny; j++) {
		for (int i = nx - 1; i > 0; i--) {
			iNode = getIndex(i,j);
			iPullNode = getIndex(i-1, j);
			fVal[iNode][3] = fVal[iPullNode][3];
		}

		for (int i = 0; i < nx - 1; i++) {
			iNode = getIndex(i, j);
			iPullNode = getIndex(i+1, j);
			fVal[iNode][7] = fVal[iPullNode][7];
		}
	}

	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny - 1; j++) {
			iNode = getIndex(i, j);
			iPullNode = getIndex(i, j+1);
			fVal[iNode][5] = fVal[iPullNode][5];
		}

		for (int j = ny - 1; j > 0; j--) {
			iNode = getIndex(i, j);
			iPullNode = getIndex(i, j-1);
			fVal[iNode][1] = fVal[iPullNode][1];
		}
	}

	for (int j = ny - 1; j > 0; j--) {
		for (int i = 0; i < nx - 1; i++) {
			iNode = getIndex(i, j);
			iPullNode = getIndex(i+1, j-1);
			fVal[iNode][8] = fVal[iPullNode][8];
		}

		for (int i = nx - 1; i > 0; i--) {
			iNode = getIndex(i, j);
			iPullNode = getIndex(i-1, j-1);
			fVal[iNode][2] = fVal[iPullNode][2];
		}
	}

	for (int j = 0; j < ny - 1; j++) {
		for (int i = 0; i < nx - 1; i++) {
			iNode = getIndex(i, j);
			iPullNode = getIndex(i+1, j+1);
			fVal[iNode][6] = fVal[iPullNode][6];
		}

		for (int i = nx - 1; i > 0; i--) {
			iNode = getIndex(i, j);
			iPullNode = getIndex(i-1, j+1);
			fVal[iNode][4] = fVal[iPullNode][4];
		}
	}
}