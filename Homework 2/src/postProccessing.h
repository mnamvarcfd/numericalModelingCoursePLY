#pragma once
#include <iostream>
#include <fstream>
#include "Grid.h"

class postProccessing
{
public:
	postProccessing();
	~postProccessing();

	void writeVectorFiedldVTK(Grid &grid);

	void writeResults(Grid & grid);

	double calMom1Res(Grid & grid);

	double totalMomentNew;
	double totalMomentOld;

};

