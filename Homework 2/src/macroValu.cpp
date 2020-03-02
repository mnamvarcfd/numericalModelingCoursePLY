
#include "pch.h"
void macroValu(const int nNode, double cx[], double cy[], double wt[], int PDFpt1[], double fVal[], double R[], double U[], double V[]) {

	double *sU = new double[nNode];
	double *sV = new double[nNode];
	double *sR = new double[nNode];

	double Cx;
	double Cy;
	double w;
	for (int k = 1; k < 9; k++) {

		Cx = cx[k];
		Cy = cy[k];
		w = wt[k];

		for (int i = k * nNode; i < (k + 1)*nNode; i++) {
			int iNode = PDFpt1[i];
			double fi = fVal[i];
			sR[iNode] = sR[iNode] + fi;
			sU[iNode] = sU[iNode] + fi * Cx;
			sV[iNode] = sV[iNode] + fi * Cy;
		}
	}


	for (int i = 0; i < nNode; i++) {
		R[i] = sR[i];
		U[i] = sU[i] / R[i];
		V[i] = sV[i] / R[i];
	}
}


//
//void firstMacro(const int i, double sR[], double sU[], double sV[], double fi, double Cx, double Cy) {
//
//	sR[i] = sR[i] + fi;
//	sU[i] = sU[i] + fi * Cx;
//	sV[i] = sV[i] + fi * Cy;
//}
//
//
//void secondMacro(const int nNode, double sR[], double sU[], double sV[], double R[], double U[], double V[]) {
//
//	for (int i = 0; i < 9 * nNode; i++) {
//		R[i] = sR[i];
//		U[i] = sU[i] / R[i];
//		V[i] = sV[i] / R[i];
//	}
//}
//

