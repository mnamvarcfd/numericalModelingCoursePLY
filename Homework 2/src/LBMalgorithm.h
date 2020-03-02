#pragma once
#include "Grid.h"
#include "LatticStencil.h"

class LBMalgorithm
{
	int nNode;
	int nx;
	int ny;
	double c1;
	double c2;
	double c3;
	double tau;

	LatticStencil stencil;

	double *R;
	double *U;
	double *V;
	double **fVal;

public:

	~LBMalgorithm();
	LBMalgorithm(Grid &grid, double tau, double cs, LatticStencil &stencil);
	int getIndex(int i, int j);

	void init();
	void calMacroValue();
	void collision();
	void streaming();
};

