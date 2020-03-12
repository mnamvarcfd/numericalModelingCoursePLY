#include "BoundaryCondition.h"
#include <iostream> 
BoundaryCondition::BoundaryCondition(Grid &grid)
{
	nx = grid.getNx();
	ny = grid.getNy();

	fVal = grid.fVal;
}

BoundaryCondition::~BoundaryCondition()
{
}

void BoundaryCondition::halfWayBouncBackBot(Grid &gridold) {

	int iNodeBot;
	for (int j = 1; j < nx - 1; j++) {
		iNodeBot = j;

		fVal[iNodeBot][2] = gridold.fVal[iNodeBot][6];
		fVal[iNodeBot][1] = gridold.fVal[iNodeBot][5];
		fVal[iNodeBot][8] = gridold.fVal[iNodeBot][4];
	}
}

void BoundaryCondition::halfWayBouncBackTop(Grid &gridold) {

	int iNodeTop;
	for (int j = 1; j < nx - 1; j++) {
		iNodeTop = (ny - 1) * nx + j;

		fVal[iNodeTop][6] = gridold.fVal[iNodeTop][2];
		fVal[iNodeTop][5] = gridold.fVal[iNodeTop][1];
		fVal[iNodeTop][4] = gridold.fVal[iNodeTop][8];
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

void BoundaryCondition::fullWayBouncBackLeft() {

	int iNodeLeft;
	for (int j = 1; j < ny - 1; j++) {
		iNodeLeft = j * nx;

		fVal[iNodeLeft][4] = fVal[iNodeLeft][8];
		fVal[iNodeLeft][3] = fVal[iNodeLeft][7];
		fVal[iNodeLeft][2] = fVal[iNodeLeft][6];
	}
}

void BoundaryCondition::fullWayBouncBackRigt() {

	int iNodeRigt;
	for (int j = 1; j < ny - 1; j++) {
		iNodeRigt = (j + 1) * nx - 1;

		fVal[iNodeRigt][8] = fVal[iNodeRigt][4];
		fVal[iNodeRigt][7] = fVal[iNodeRigt][3];
		fVal[iNodeRigt][6] = fVal[iNodeRigt][2];
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




void BoundaryCondition::ZouHeVelLeft(double ux, double uy) {

	int iNode;
	for (int j = 1; j < ny - 1; j++) {
		iNode = j * nx;

		f0 = fVal[iNode][0];
		f1 = fVal[iNode][1];
		f2 = fVal[iNode][2];
		f3 = fVal[iNode][3];
		f4 = fVal[iNode][4];
		f5 = fVal[iNode][5];
		f6 = fVal[iNode][6];
		f7 = fVal[iNode][7];
		f8 = fVal[iNode][8];

		double rho = (f0 + f1 + f5 + 2 * (f6 + f7 + f8)) / (1 - ux);

		//std::cout << (f0 + f1 + f5 + 2 * (f6 + f7 + f8))  << "====" <<  (1 - ux) << std::endl;

		fVal[iNode][3] = f7 - (2. / 3)*rho*ux;
		fVal[iNode][2] = f6 + 0.5*(f5 - f1) + 0.5*rho * uy + (1. / 6)*rho*ux;
		fVal[iNode][4] = rho * ux + (f6 + f7 + f8) - (fVal[iNode][2] + fVal[iNode][3]);

		//f0 = fVal[iNode][0];
		//f1 = fVal[iNode][1];
		//f2 = fVal[iNode][2];
		//f3 = fVal[iNode][3];
		//f4 = fVal[iNode][4];
		//f5 = fVal[iNode][5];
		//f6 = fVal[iNode][6];
		//f7 = fVal[iNode][7];
		//f8 = fVal[iNode][8];
		//std::cout << (f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8) << "==rho==" << rho << std::endl;
		//std::cout << (f2 + f3 + f4 - f6 - f7 - f8) << "==rho*Uwall==" << rho* ux << std::endl;

		//std::cin >> f0;
	}

}
void BoundaryCondition::ZouHeVelRigt(double ux, double uy) {

	int iNode;
	for (int j = 1; j < ny - 1; j++) {
		iNode = (j + 1) * nx - 1;

		f0 = fVal[iNode][0];
		f1 = fVal[iNode][1];
		f2 = fVal[iNode][2];
		f3 = fVal[iNode][3];
		f4 = fVal[iNode][4];
		f5 = fVal[iNode][5];
		f6 = fVal[iNode][6];
		f7 = fVal[iNode][7];
		f8 = fVal[iNode][8];

		double rho = (f0 + f1 + f5 + 2 * (f2 + f3 + f4)) / (1 + ux);

		fVal[iNode][7] = f3 - (2. / 3)*rho*ux;
		fVal[iNode][6] = f2 + 0.5*(f1 - f5) - 0.5*rho * uy - (1. / 6)*rho*ux;

		fVal[iNode][8] = f2 + f3 + f4 - (f6 + f7) - rho * ux;

		//f0 = fVal[iNode][0];
		//f1 = fVal[iNode][1];
		//f2 = fVal[iNode][2];
		//f3 = fVal[iNode][3];
		//f4 = fVal[iNode][4];
		//f5 = fVal[iNode][5];
		//f6 = fVal[iNode][6];
		//f7 = fVal[iNode][7];
		//f8 = fVal[iNode][8];
		//std::cout << (f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8) << "==rho==" << rho << std::endl;
		//std::cout << (f2 + f3 + f4 - f6 - f7 - f8) << "==rho*Uwall==" << rho* ux << std::endl;


	}

}
void BoundaryCondition::ZouHeVelTop(double ux) {
	double uy = 0.0;
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

	    double rho = (f0 + f3 + f7 + 2 * (f1 + f2 + f8) ) / (1+uy);
		fVal[iNode][5] = f1 - (2. / 3)*rho*uy;
		fVal[iNode][4] = (1. / 3)*rho*uy + 0.5*rho*(1 + ux) - (0.5*f0 + f1 + f2 + f3);
		fVal[iNode][6] = f2 + f3 + fVal[iNode][4] - f7 - f8 - rho * ux;

/*
		f0 = fVal[iNode][0];
		f1 = fVal[iNode][1];
		f2 = fVal[iNode][2];
		f3 = fVal[iNode][3];
		f4 = fVal[iNode][4];
		f5 = fVal[iNode][5];
		f6 = fVal[iNode][6];
		f7 = fVal[iNode][7];
		f8 = fVal[iNode][8];
		std::cout << (f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8) << "==rho==" << rho << std::endl;
		std::cout << (f2 + f3 + f4 - f6 - f7 - f8) << "==rho*Uwall==" << rho* ux << std::endl;*/
		//std::cout << (f1 + f2 + f8 - f4 - f5 - f6) << "=====" << Uwall << std::endl;
	}

}
void BoundaryCondition::ZouHeVelBot(double ux) {
	double uy = 0.0;
	int iNode;
	for (int j = 1; j < nx - 1; j++) {
		iNode = j;

		f0 = fVal[iNode][0];
		f1 = fVal[iNode][1];
		f2 = fVal[iNode][2];
		f3 = fVal[iNode][3];
		f4 = fVal[iNode][4];
		f5 = fVal[iNode][5];
		f6 = fVal[iNode][6];
		f7 = fVal[iNode][7];
		f8 = fVal[iNode][8];

		double rho = ( f0 + f3 + f7 + 2 * (f4 + f5 + f6) ) / (1 - uy);

		fVal[iNode][1] = f5 - (2./3)*rho*uy;
		fVal[iNode][8] = 0.5*rho * (1-ux) - (1./3)*rho*uy - (0.5*f0 + f5 + f6 + f7 +f8);
		fVal[iNode][2] = rho * ux - f3 - f4 + (f6 + f7 + fVal[iNode][8]);

		//f0 = fVal[iNode][0];
		//f1 = fVal[iNode][1];
		//f2 = fVal[iNode][2];
		//f3 = fVal[iNode][3];
		//f4 = fVal[iNode][4];
		//f5 = fVal[iNode][5];
		//f6 = fVal[iNode][6];
		//f7 = fVal[iNode][7];
		//f8 = fVal[iNode][8];
		//std::cout << (f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8) << "==rho==" << rho << std::endl;
		//std::cout << (f2 + f3 + f4 - f6 - f7 - f8) << "==rho*Uwall==" << rho* ux << std::endl;


	}

}


void BoundaryCondition::zouHeVelTopLeftCorner(double rho, double ux) {
	double uy = 0.0;
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

	fVal[iNode][3] = f7 + (2. / 3)*rho*ux;
	fVal[iNode][5] = f1 - (2. / 3)*rho*uy;
	fVal[iNode][4] = f8 + (1. / 6)*rho*(ux - uy);
	fVal[iNode][6] = 0.5*rho*(1 - ux) - (0.5*f0 + f1 + f7 + f8) + (1./ 3)*rho*uy;
	fVal[iNode][2] = 0.5*rho*(uy+1) - (1. / 3)*rho*ux - (0.5*f0 + f1 + f7 + f8);

	//f0 = fVal[iNode][0];
	//f1 = fVal[iNode][1];
	//f2 = fVal[iNode][2];
	//f3 = fVal[iNode][3];
	//f4 = fVal[iNode][4];
	//f5 = fVal[iNode][5];
	//f6 = fVal[iNode][6];
	//f7 = fVal[iNode][7];
	//f8 = fVal[iNode][8];
	//std::cout << (f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8) << "==rho==" << rho << std::endl;
	//std::cout << (f2 + f3 + f4 - f6 - f7 - f8) << "==rho*ux==" << rho* ux << std::endl;
	//std::cout << (f1 + f2 + f8 - f4 - f5 - f6) << "==rho* uy==" << Uwall << std::endl;
}
void BoundaryCondition::zouHeVelTopRigtCorner(double rho, double ux) {
	double uy = 0.0;
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

	fVal[iNode][7] = f3 - (2. / 3)*rho*ux;
	fVal[iNode][5] = f1 - (2. / 3)*rho*uy;
	fVal[iNode][6] = f2 - (1. / 6)*rho*(ux+uy);
	fVal[iNode][4] = 0.5*rho*(ux + 1) + (1. / 3)*rho*uy - (0.5*f0 + f1 + f2 + f3);
	fVal[iNode][8] = 0.5*rho*(uy + 1) + (1. / 3)*rho*ux - (0.5*f0 + f1 + f2 + f3);


	//f0 = fVal[iNode][0];
	//f1 = fVal[iNode][1];
	//f2 = fVal[iNode][2];
	//f3 = fVal[iNode][3];
	//f4 = fVal[iNode][4];
	//f5 = fVal[iNode][5];
	//f6 = fVal[iNode][6];
	//f7 = fVal[iNode][7];
	//f8 = fVal[iNode][8];

	//std::cout << (f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8) << "==rho==" << rho << std::endl;
	//std::cout << (f2 + f3 + f4 - f6 - f7 - f8) << "==rho*ux==" << rho * ux << std::endl;

	//std::cin >> f0;
}
void BoundaryCondition::zouHeVelBotLeftCorner(double rho, double ux) {
	double uy = 0.0;
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

	fVal[iNode][3] = f7 + (2. / 3)*rho*ux;
	fVal[iNode][1] = f5 + (2. / 3)*rho*uy;
	fVal[iNode][2] = f6 + (1. / 6)*rho*(ux + uy);
	fVal[iNode][8] = 0.5*rho - 0.5*rho*ux - (1. / 3)*rho*uy - (0.5*f0 + f5 + f6 + f7);
	fVal[iNode][4] = 0.5*rho - (1. / 3)*rho*ux - 0.5*rho*uy - (0.5*f0 + f5 + f6 + f7);


	//f0 = fVal[iNode][0];
	//f1 = fVal[iNode][1];
	//f2 = fVal[iNode][2];
	//f3 = fVal[iNode][3];
	//f4 = fVal[iNode][4];
	//f5 = fVal[iNode][5];
	//f6 = fVal[iNode][6];
	//f7 = fVal[iNode][7];
	//f8 = fVal[iNode][8];
	//std::cout << (f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8) << "==rho==" << rho << std::endl;
	//std::cout << (f2 + f3 + f4 - f6 - f7 - f8) << "==rho*ux==" << rho * ux << std::endl;
	//if ((f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8) != rho) {

	//	std::cout << "f0    " << f0 << std::endl;
	//	std::cout << "f0    " << f1 << std::endl;
	//	std::cout << "f0    " << f2 << std::endl;
	//	std::cout << "f0    " << f3 << std::endl;
	//	std::cout << "f0    " << f4 << std::endl;
	//	std::cout << "f0    " << f5 << std::endl;
	//	std::cout << "f0    " << f6 << std::endl;
	//	std::cout << "f0    " << f7 << std::endl;
	//	std::cout << "f0    " << f8 << std::endl;
	//}
}
void BoundaryCondition::zouHeVelBotRigtCorner(double rho, double ux) {

	int iNode = ny * nx - 1;
	double uy = 0.0;

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
	fVal[iNode][1] = f5;
	fVal[iNode][6] = 0.5*rho*(1 - uy) - (0.5*f0 + f3 + f4 + f5);
	fVal[iNode][8] = f4 - 0.5*rho*(ux-uy);
	fVal[iNode][2] = 0.5*rho*(1 + ux) - (0.5*f0 + f3 + f4 + f5);


	//f0 = fVal[iNode][0];
	//f1 = fVal[iNode][1];
	//f2 = fVal[iNode][2];
	//f3 = fVal[iNode][3];
	//f4 = fVal[iNode][4];
	//f5 = fVal[iNode][5];
	//f6 = fVal[iNode][6];
	//f7 = fVal[iNode][7];
	//f8 = fVal[iNode][8];
	//std::cout << (f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8) << "==rho==" << rho << std::endl;
	//std::cout << (f2 + f3 + f4 - f6 - f7 - f8) << "==rho*ux==" << rho * ux << std::endl;
}





