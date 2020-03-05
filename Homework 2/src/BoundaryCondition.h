#pragma once
#include "Grid.h"
class BoundaryCondition
{

	double **fVal;

	int nx;
	int ny;
	double f0;
	double f1;
	double f2;
	double f3;
	double f4;
	double f5;
	double f6;
	double f7;
	double f8;

public:
	BoundaryCondition(Grid &grid);
	~BoundaryCondition();

	void halfWayBouncBackBot();

	void halfWayBouncBackTop();

	void halfWayBouncBack();

	void fullWayBouncBackBot();

	void fullWayBouncBackTop();

	void periodic();

	void zouHeLeft(double rho);

	void zouHeRigt(double rho);

	void zouHeBottomLeftCorner(double rho);

	void zouHeTopLeftCorner(double rho);

	void zouHeBottomRightCorner(double rho);

	void zouHeTopRightCorner(double rho);

	void ZouHeVelTop(double Uwall);

	void zouHeVelTopLeftCorner(double rho, double vel);

	void zouHeVelTopRigtCorner(double rho, double ux);

};

