#include "BoundaryCondition.h"

BoundaryCondition::BoundaryCondition(Grid &grid)
{
	nx = grid.getNx();
	ny = grid.getNy();

	fVal = grid.fVal;
}

BoundaryCondition::~BoundaryCondition()
{
}

void BoundaryCondition::halfWayBouncBackBot() {

	int iNodeBot;
	for (int j = 1; j < nx - 1; j++) {
		iNodeBot = j;

		fVal[iNodeBot][2] = fVal[iNodeBot][6];
		fVal[iNodeBot][1] = fVal[iNodeBot][5];
		fVal[iNodeBot][8] = fVal[iNodeBot][4];
	}
}

void BoundaryCondition::halfWayBouncBackTop() {

	int iNodeTop;
	for (int j = 1; j < nx - 1; j++) {
		iNodeTop = (ny - 1) * nx + j;

		fVal[iNodeTop][6] = fVal[iNodeTop][2];
		fVal[iNodeTop][5] = fVal[iNodeTop][1];
		fVal[iNodeTop][4] = fVal[iNodeTop][8];
	}
}

void BoundaryCondition::fullWayBouncBackBot() {

	int iNodeTop;
	int iNodeBot;
	for (int j = 1; j < nx-1; j++) {
		iNodeBot = j;

		fVal[iNodeBot][2] = fVal[iNodeBot][6];
		fVal[iNodeBot][1] = fVal[iNodeBot][5];
		fVal[iNodeBot][8] = fVal[iNodeBot][4];
	}
}

void BoundaryCondition::fullWayBouncBackTop() {

	int iNodeTop;
	for (int j = 1; j < nx - 1; j++) {
		iNodeTop = (ny - 1) * nx + j;

		fVal[iNodeTop][6] = fVal[iNodeTop][2];
		fVal[iNodeTop][5] = fVal[iNodeTop][1];
		fVal[iNodeTop][4] = fVal[iNodeTop][8];
	}
}

void BoundaryCondition::periodic() {

	int iNodeLeft;
	int iNodeRigt;
	for (int j = 1; j < ny - 1; j++) {
		iNodeLeft = j * nx;
		iNodeRigt = (j + 1) * nx - 1;

		fVal[iNodeRigt][8] = fVal[iNodeLeft][8];
		fVal[iNodeRigt][7] = fVal[iNodeLeft][7];
		fVal[iNodeRigt][6] = fVal[iNodeLeft][6];

		fVal[iNodeLeft][2] = fVal[iNodeRigt][2];
		fVal[iNodeLeft][3] = fVal[iNodeRigt][3];
		fVal[iNodeLeft][4] = fVal[iNodeRigt][4];
	}
}

void BoundaryCondition::zouHeLeft(double rho) {

	int iNodeLeft;
	for (int j = 1; j < ny-1; j++) {
		iNodeLeft = j * nx;

		f0 = fVal[iNodeLeft][0];
		f1 = fVal[iNodeLeft][1];
		f2 = fVal[iNodeLeft][2];
		f3 = fVal[iNodeLeft][3];
		f4 = fVal[iNodeLeft][4];
		f5 = fVal[iNodeLeft][5];
		f6 = fVal[iNodeLeft][6];
		f7 = fVal[iNodeLeft][7];
		f8 = fVal[iNodeLeft][8];

		double u = 1.0 - (f0 + f1 + f5 + 2 * (f6 + f7 + f8)) / rho;

		fVal[iNodeLeft][3] = f7 + (2. / 3)*rho*u;
		fVal[iNodeLeft][2] = f6 + 0.5*(f5 - f1) + (1. / 6)*rho*u;
		fVal[iNodeLeft][4] = f8 - 0.5*(f5 - f1) + (1. / 6)*rho*u;
	}

	zouHeBottomLeftCorner(rho);

	zouHeTopLeftCorner(rho);
}

void BoundaryCondition::zouHeRigt(double rho) {

	int iNodeRigt;
	for (int j = 1; j < ny-1; j++) {
		iNodeRigt = (j + 1) * nx - 1;

		f0 = fVal[iNodeRigt][0];
		f1 = fVal[iNodeRigt][1];
		f2 = fVal[iNodeRigt][2];
		f3 = fVal[iNodeRigt][3];
		f4 = fVal[iNodeRigt][4];
		f5 = fVal[iNodeRigt][5];
		f6 = fVal[iNodeRigt][6];
		f7 = fVal[iNodeRigt][7];
		f8 = fVal[iNodeRigt][8];

		double u = (f0 + f1 + f5 + 2 * (f2 + f3 + f4)) / rho - 1.0;

		fVal[iNodeRigt][7] = f3 - (2. / 3)*rho*u;
		fVal[iNodeRigt][6] = f2 + 0.5*(f1 - f5) - (1. / 6)*rho*u;
		fVal[iNodeRigt][8] = f4 - 0.5*(f1 - f5) - (1. / 6)*rho*u;
	}

	zouHeBottomRightCorner(rho);

	zouHeTopRightCorner(rho);
}


void BoundaryCondition::zouHeBottomLeftCorner(double rho) {

	int iNode=0;

	f0 = fVal[iNode][0];
	f1 = fVal[iNode][1];
	f2 = fVal[iNode][2];
	f3 = fVal[iNode][3];
	f4 = fVal[iNode][4];
	f5 = fVal[iNode][5];
	f6 = fVal[iNode][6];
	f7 = fVal[iNode][7];
	f8 = fVal[iNode][8];

	fVal[iNode][1] = f5;
	fVal[iNode][2] = f6;
	fVal[iNode][3] = f7;
	fVal[iNode][4] = 0.5*(rho - (f0 + 2 * f5 + 2 * f6 + 2 * f7));
	fVal[iNode][8] = fVal[iNode][4];
}


void BoundaryCondition::zouHeTopLeftCorner(double rho) {

	int iNode = (ny - 1) * nx;

	f0 = fVal[iNode][0];
	f1 = fVal[iNode][1];
	f2 = fVal[iNode][2];
	f3 = fVal[iNode][3];
	f4 = fVal[iNode][4];
	f5 = fVal[iNode][5];
	f6 = fVal[iNode][6];
	f7 = fVal[iNode][7];
	f8 = fVal[iNode][8];

	fVal[iNode][2] = 0.5*(rho - f0) - (f1 + f7 + f8);
	fVal[iNode][3] = f7;
	fVal[iNode][4] = f8;
	fVal[iNode][5] = f1;
	fVal[iNode][6] = fVal[iNode][2];
}
 
void BoundaryCondition::zouHeBottomRightCorner(double rho) {

	int iNode = nx - 1;

	f0 = fVal[iNode][0];
	f1 = fVal[iNode][1];
	f2 = fVal[iNode][2];
	f3 = fVal[iNode][3];
	f4 = fVal[iNode][4];
	f5 = fVal[iNode][5];
	f6 = fVal[iNode][6];
	f7 = fVal[iNode][7];
	f8 = fVal[iNode][8];
	 
	fVal[iNode][1] = f5;
	fVal[iNode][2] = 0.5*(rho - f0) - (f3 + f4 + f5);
	fVal[iNode][6] = fVal[iNode][2];
	fVal[iNode][7] = f3;
	fVal[iNode][8] = f4;
}
 
void BoundaryCondition::zouHeTopRightCorner(double rho) {

	int iNode = ny * nx-1;

	f0 = fVal[iNode][0];
	f1 = fVal[iNode][1];
	f2 = fVal[iNode][2];
	f3 = fVal[iNode][3];
	f4 = fVal[iNode][4];
	f5 = fVal[iNode][5];
	f6 = fVal[iNode][6];
	f7 = fVal[iNode][7];
	f8 = fVal[iNode][8];
 
	fVal[iNode][4] = 0.5*(rho - f0) - (f1 + f2 + f3);
	fVal[iNode][5] = f1;
	fVal[iNode][6] = f2;
	fVal[iNode][7] = f3;
	fVal[iNode][8] = fVal[iNode][4];
}



void BoundaryCondition::ZouHeVelTop(double Uwall) {

	int iNode;
	for (int j = 1; j < nx - 1; j++) {
		iNode = (ny - 1) * nx + j;

		f0 = fVal[iNode][0];
		f1 = fVal[iNode][1];
		f2 = fVal[iNode][2];
		f3 = fVal[iNode][3];
		f4 = fVal[iNode][4];
		f5 = fVal[iNode][5];
		f6 = fVal[iNode][6];
		f7 = fVal[iNode][7];
		f8 = fVal[iNode][8];

	    double rho = f0 + 2 * f1 + 2 * f2 + f3 + 2 * f8 + f7;

		fVal[iNode][5] = f1;
		fVal[iNode][6] = 0.5*(2 * f2 + f3 - f7 - rho * Uwall);
		fVal[iNode][4] = f2 + f8 - fVal[iNode][6];

	}

}



void BoundaryCondition::zouHeVelTopLeftCorner(double rho, double ux) {

	int iNode = (ny - 1) * nx;

	f0 = fVal[iNode][0];
	f1 = fVal[iNode][1];
	f2 = fVal[iNode][2];
	f3 = fVal[iNode][3];
	f4 = fVal[iNode][4];
	f5 = fVal[iNode][5];
	f6 = fVal[iNode][6];
	f7 = fVal[iNode][7];
	f8 = fVal[iNode][8];

	fVal[iNode][3] = f7;
	fVal[iNode][5] = f1;
	fVal[iNode][4] = f8+0.5*rho*ux;
	fVal[iNode][6] = (f1 + f3 + f8) - 0.5*rho*(ux - 1);
	fVal[iNode][2] = (f1 + f3 + f8) + 0.5*rho;
}


void BoundaryCondition::zouHeVelTopRigtCorner(double rho, double ux) {

	int iNode = ny * nx - 1;

	f0 = fVal[iNode][0];
	f1 = fVal[iNode][1];
	f2 = fVal[iNode][2];
	f3 = fVal[iNode][3];
	f4 = fVal[iNode][4];
	f5 = fVal[iNode][5];
	f6 = fVal[iNode][6];
	f7 = fVal[iNode][7];
	f8 = fVal[iNode][8];

	fVal[iNode][7] = f3;
	fVal[iNode][5] = f1;
	fVal[iNode][6] = f2 - 0.5*rho*ux;
	fVal[iNode][4] = 0.5*rho*(ux + 1) - (f1 + f2 + f3);
	fVal[iNode][8] = 0.5*rho - (f1 + f2 + f3);
}
