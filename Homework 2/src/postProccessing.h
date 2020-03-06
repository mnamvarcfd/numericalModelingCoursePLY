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

	void writeScalarFieldVTK(Grid & grid);

	void writeResults(Grid & grid);

	double calMom1Res(Grid & grid);

	double calMassRes(Grid & grid);

	double totalMomentNew;
	double totalMomentOld;

	double totalMassNew;
	double totalMassOld;

};

