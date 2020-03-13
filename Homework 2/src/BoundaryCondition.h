#pragma once
#include "Grid.h"
class BoundaryCondition
{

	double **fVal;

	double *rhoTopLeftCorn;	
	double *rhoBotLeftCorn;
	double *rhoBotRigtCorn;
	double *rhoTopRigtCorn;

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

	void halfWayBouncBackBot(Grid & gridold);

	void halfWayBouncBackTop(Grid & gridold);


	void fullWayBouncBackBot();

	void fullWayBouncBackTop();

	void fullWayBouncBackLeft();

	void fullWayBouncBackRigt();

	void periodic();

	void zouHeLeft(double rho);

	void zouHeRigt(double rho);

	void zouHeBottomLeftCorner(double rho);

	void zouHeTopLeftCorner(double rho);

	void zouHeBottomRightCorner(double rho);

	void zouHeTopRightCorner(double rho);

	void ZouHeVelTop(double ux, double uy);
	void ZouHeVelBot(double ux, double uy);
	void ZouHeVelLeft(double ux, double uy);
	void ZouHeVelRigt(double ux, double uy);

	void zouHeVelTopLeftCorner(double ux, double uy);
	void zouHeVelTopRigtCorner(double ux, double uy);
	void zouHeVelBotLeftCorner(double ux, double uy);
	void zouHeVelBotRigtCorner(double ux, double uy);



};

