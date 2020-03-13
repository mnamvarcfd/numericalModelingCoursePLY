#include <iostream>
#include <sstream>
#include "postProccessing.h"

using namespace std;

postProccessing::postProccessing()
{
	totalMomentNew = 0.0;
	totalMomentOld = 0.0;

	totalMassNew = 0.0;
	totalMassOld = 0.0;
}


postProccessing::~postProccessing()
{
}

void postProccessing::writeVectorFiedldVTK(Grid &grid)
{

	// Establish the file name
	string filename = "velocity.vtk";

	// Declare the output flux and open the file
	ofstream file(filename.c_str(), ios::out | ios::trunc);
	cout << "Saving file: " << filename << endl;

	if (file)  // File has been opened correctly
	{
		file << "# vtk DataFile Version 2.0" << endl;
		file << "Test" << endl;
		file << "ASCII" << endl;
		file << "DATASET STRUCTURED_POINTS" << endl;
		file << "DIMENSIONS" << " " << grid.getNx() << " " << grid.getNy() <<
			" " << 1 << endl;
		file << "ORIGIN" << " " << grid.getXMin() << " " << grid.getYMin() <<
			" " << 0. << endl;
		file << "SPACING" << " " << grid.getDx() << " " << grid.getDx() <<
			" " << 1. << endl;
		file << "POINT_DATA" << " " << grid.getNNode() << endl;
		file << "VECTORS" << " " << "Velocity" << " " << "double" << endl;
		for (int i = 0; i < grid.getNNode(); ++i)
		{
			file << grid.U[i] << " " << grid.V[i] << " " << 0. << endl;
		}


		file.close();  // close the file
	}
}


void postProccessing::writeScalarFieldVTK(Grid &grid)
{

	// Establish the file name
	string filename = "density.vtk";

	// Declare the output flux and open the file
	ofstream file(filename.c_str(), ios::out | ios::trunc);
	cout << "Saving file: " << filename << endl;

	if (file)  // File has been opened correctly
	{
		file << "# vtk DataFile Version 2.0" << endl;
		file << "Test" << endl;
		file << "ASCII" << endl;
		file << "DATASET STRUCTURED_POINTS" << endl;
		file << "DIMENSIONS" << " " << grid.getNx() + 1 << " " << grid.getNy() + 1 <<
			" " << 1 << endl;
		file << "ORIGIN" << " " << grid.getXMin() << " " << grid.getYMin() <<
			" " << 0. << endl;
		file << "SPACING" << " " << grid.getDx() << " " << grid.getDx() <<
			" " << 1. << endl;
		file << "CELL_DATA" << " " << grid.getNNode() << endl;
		file << "SCALARS" << " " << "ScalarValue" << " " << "double" << " " << 1 << endl;
		file << "LOOKUP_TABLE default" << endl;
		for (int i = 0; i < grid.getNNode(); ++i)
		{
			file << grid.R[i] << endl;

		}
		file.close();  // close the file
	}

}

void postProccessing::writeResults(Grid &grid) {

	FILE *file1;
	fopen_s(&file1, "Results.txt", "w");
	fprintf(file1, " i  j       f0           ");
	fprintf(file1, "f1          f2            f3           ");
	fprintf(file1, "f4             f5            f6           ");
	fprintf(file1, "f7            f8  \n");

	int iNode = 0;
	for (int j = 0; j < grid.getNy(); j++) {
		for (int i = 0; i < grid.getNx(); i++) {
			iNode = j * grid.getNx() + i;

			fprintf(file1, "%d  %d ", i, j);


			for (int k = 0; k < 9; k++) {
				fprintf(file1, "%.9f  ", grid.fVal[iNode][k]);
			}
			fprintf(file1, "\n");


			iNode++;
		}
	}


	fprintf(file1, " i  j          ro           u             v \n ");

	iNode = 0;
	for (int j = 0; j < grid.getNy(); j++) {
		for (int i = 0; i < grid.getNx(); i++) {
			iNode = j * grid.getNx() + i;

			fprintf(file1, "%d  %d ", i, j);


			fprintf(file1, "%.9f  ", grid.R[iNode]);
			fprintf(file1, "%.9f  ", grid.U[iNode]);
			fprintf(file1, "%.9f  ", grid.V[iNode]);
			fprintf(file1, "  \n ");

			iNode++;
		}
	}
	fclose(file1);
}


double postProccessing::calMom1Res(Grid &grid) {

	totalMomentNew = 0.;
	for (int j = 0; j < grid.getNNode(); j++)
		totalMomentNew += grid.U[j];

	double residual = abs(totalMomentNew - totalMomentOld) / abs(totalMomentNew*grid.getNNode());
	//std::cout << totalMomentNew << "===calMassRes==" << residual << std::endl;
	
	totalMomentOld = totalMomentNew;

	return residual;
}


double postProccessing::calMassRes(Grid &grid) {

	totalMassNew = 0.;
	for (int j = 0; j < grid.getNNode(); j++)
		totalMassNew += grid.R[j];

	double residual = abs(totalMassNew - totalMassOld) / (totalMassNew*grid.getNNode());
	std::cout << totalMassNew << "===calMassRes==" << residual << std::endl;
	totalMassOld = totalMassNew;

	return residual;
}