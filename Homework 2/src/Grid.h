#pragma once
class Grid
{

private:

	int nx;
	int ny;
	int nNode;

	double dx;

	double xMin;
	double xMax;
	double yMin;
	double yMax;

public:

	double *R;
	double *U;
	double *V;
	double **fVal;

	Grid(int nx_, int ny_, int nNode_, double dx_, double xMin_, double xMax_, double yMin_, double yMax_);
	~Grid();

	double getDx();
	double getXMin();
	double getXMax();
	double getYMin();
	double getYMax();

	int getNx();
	int getNy();
	int getNNode();
};

