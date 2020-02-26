#include "lattice.h"
#include "domain.h"
#include <iostream>

void writeResults(Domain domain, Lattice lat) {

	FILE *file1;
	fopen_s(&file1, "Results.txt", "w");
	fprintf(file1, " i  j       f0           ");
	fprintf(file1, "f1          f2            f3           ");
	fprintf(file1, "f4             f5            f6           ");
	fprintf(file1, "f7            f8  \n");

	int iNode = 0;
	for (int j = 0; j < domain.getNy(); j++) {
		for (int i = 0; i < domain.getNx(); i++) {
			iNode = j * domain.getNx() + i;

			fprintf(file1, "%d  %d ", i, j);


			for (int k = 0; k < 9; k++) {
				fprintf(file1, "%.9f  ", lat.f_[iNode][k]);
			}
			fprintf(file1, "\n");


			iNode++;
		}
	}


	fprintf(file1, " i  j          ro           u             v \n ");

	iNode = 0;
	for (int j = 0; j < domain.getNy(); j++) {
		for (int i = 0; i < domain.getNx(); i++) {
			iNode = j * domain.getNx() + i;

			fprintf(file1, "%d  %d ", i, j);


			fprintf(file1, "%.9f  ", lat.rho_[iNode]);
			fprintf(file1, "%.9f  ", lat.u_[iNode][0]);
			fprintf(file1, "%.9f  ", lat.u_[iNode][1]);
			fprintf(file1, "  \n ");

			iNode++;
		}
	}
	fclose(file1);
}
