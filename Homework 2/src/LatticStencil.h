#pragma once
class LatticStencil
{

public:
	double *cx = new double[9];
	double *cy = new double[9];
	double *wt = new double[9];

	LatticStencil();
	void setStencil();
	~LatticStencil();
};

