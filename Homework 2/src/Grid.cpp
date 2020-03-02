#include "Grid.h"

Grid::Grid(int nx_, int ny_, int nNode_, double dx_, double xMin_, double xMax_, double yMin_, double yMax_)
{
	nx=nx_;
	ny=ny_;
	nNode=nNode_;

	dx=dx_;

	xMin = xMin_;
	xMax = xMax_;
	yMin = yMin_;
	yMax = yMax_;


	R = new double[nNode];
	U = new double[nNode];
	V = new double[nNode];

	fVal = new double*[nNode];
	for (int i = 0; i < nNode; i++)
		fVal[i] = new double[9];

}


Grid::~Grid()
{
}

double Grid::getDx() {    
	return dx; 
}

double Grid::getXMin(){	
	return 	xMin; 
}

double Grid::getXMax() { 
	return 	xMax; 
}

double Grid::getYMin() { 
	return 	yMin; 
}

double Grid::getYMax() { 
	return 	yMax; 
}

int Grid::getNx() { 
	return 	nx; 
}

int Grid::getNy() { return 	
nx; ny; 
}

int Grid::getNNode() { 
	return 	nNode; 
}
