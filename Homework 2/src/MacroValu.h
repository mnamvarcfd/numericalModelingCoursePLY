#pragma once
#include "LatticMesh.h"
#include "LatticStencil.h"

#include <stdio.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

class MacroValu {

public:
LatticMesh mesh;
LatticStencil stencil;

double *fVal = new double[mesh.nPDF];
double *fValOld = new double[mesh.nPDF];

double *U = new double[mesh.nNode];
double *V = new double[mesh.nNode];
double *R = new double[mesh.nNode];

double *sU = new double[mesh.nNode];
double *sV = new double[mesh.nNode];
double *sR = new double[mesh.nNode];

double omega;

void initMacroValu(LatticMesh mesh, LatticStencil stencil) {
	double nu = 0.166666666666667;
	omega = 2 / (6 * nu + 1);

	for (int k = 0; k < stencil.nVec; k++) {
		for (int i = k * mesh.nNode; i < (k + 1)*mesh.nNode; i++) {
			fVal[i] = stencil.wt[k];
		}
	}

	for (int i = 0; i < mesh.nNode; i++) {
		U[i] = 0.0;
		V[i] = 0.0;
		R[i] = 1.0;
	}
}


void calMacroValu(LatticMesh mesh, LatticStencil stencil) {

	for (int i = 0; i < mesh.nNode; i++) {
		R[i] = sR[i];
		U[i] = sU[i] / sR[i];
		V[i] = sV[i] / sR[i];

		//cout << R[i] << " " << U[i] << " " << V[i] << endl;

		sR[i] = 0.0;
		sU[i] = 0.0;
		sV[i] = 0.0;
	}
	

		//
		//system("pause");
}





void writeMacroInTecplot(LatticMesh mesh) {

	FILE *file1;
	fopen_s(&file1, "macro.plt", "w");

	fprintf(file1, " TITLE = 'heat'  \n");
	fprintf(file1, " VARIABLES = 'X', 'Y' , 'U', 'V', 'Ro'   \n");
	fprintf(file1, "ZONE T = \"heat zone\", I = %d , J = %d , F = POINT \n", mesh.Nx, mesh.Ny);

	int cnt = 0;
	for (int j = 0; j < mesh.Ny; j++) {
		for (int i = 0; i < mesh.Nx; i++) {

			fprintf(file1, "%.9f  %.9f  %.9f   %.9f  %.9f \n ", mesh.X[cnt], mesh.Y[cnt], U[cnt], V[cnt], R[cnt]);
			cnt++;
		}
	}
	fclose(file1);
}




void writeVelocityProfile(LatticMesh mesh) {

	FILE *file1;
	fopen_s(&file1, "VelocityProfile.plt", "w");

	double *ux_Numeric = new double[mesh.Ny];
	double width = (mesh.Ny - 2) / 2.;

	int cnt = 0;
	for (int i = 0; i < mesh.nNode; i++) {
		if (mesh.X[i]<0.45 && mesh.X[i] > 0.40) {

			double y = ((double)cnt - width);
			ux_Numeric[cnt] = U[i];
			cnt++;
		}
	}

	double uCenter = 0.0;
	for (int i = 0; i < cnt; i++)
		if (ux_Numeric[i] > uCenter) uCenter = ux_Numeric[i];


	for (int i = 1; i < cnt-1; i++) {
			double y = ((double)i - width);
			ux_Numeric[i] = ux_Numeric[i] / uCenter;

			fprintf(file1, "%.9f  %.9f  \n ", ux_Numeric[i], (y -0.5)/ width);

	}
	fclose(file1);

}



void writeResults(LatticMesh mesh) {

	FILE *file1;
	fopen_s(&file1, "Results.txt", "w");
	fprintf(file1, " i  j       f0           ");
	fprintf(file1, "f1          f2            f3           ");
	fprintf(file1, "f4             f5            f6           ");
	fprintf(file1, "f7            f8          ro           ");
	fprintf(file1, "u             v \n ");
	int iNode = 0;
	for (int j = 0; j < mesh.Ny; j++) {
		for (int i = 0; i < mesh.Nx; i++) {

			fprintf(file1, "%d  %d ", i, j);

			for (int k = 0; k < 9; k++) {
				int iPDF = mesh.NodePDF[9 * iNode + k];
				fprintf(file1, "%.9f  ", fVal[iPDF]);
				//cout << fVal[iPDF] << "=-----=" << iPDF << endl;
			}

			fprintf(file1, "%.9f  ", R[iNode]);
			fprintf(file1, "%.9f  ", U[iNode]);
			fprintf(file1, "%.9f  ", V[iNode]);
			fprintf(file1, "  \n ");

			iNode++;
		}
	}
	fclose(file1);
}




};