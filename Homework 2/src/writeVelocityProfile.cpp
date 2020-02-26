#include <fstream>
#include "lattice.h"
#include "domain.h"
void writeVelocityProfile(std::string fileName, Lattice lat, Domain domain) {

	int nx = domain.getNx();
	int ny = domain.getNy();
	int ntot = domain.getNTot();

	FILE *file1;
	fopen_s(&file1, fileName.c_str(), "w");

	double *ux_Numeric = new double[nx];
	double width = (nx - 1) / 2.;

	int cnt = 0;
	for (int i = 0; i < ntot; i++) {
		if (nx <= i && i < 2 * nx) {
			ux_Numeric[cnt] = -lat.u_[i][1];
			cnt++;
		}
	}

	double uCenter = 0.000000000001;
	for (int i = 0; i < cnt; i++) {
		if (ux_Numeric[i] > uCenter) uCenter = ux_Numeric[i];
	}

	fprintf(file1, "zone  \n ");
	for (int i = 0; i < cnt; i++) {
		double y = ((double)i - width);
		ux_Numeric[i] = ux_Numeric[i] /*/ uCenter*/;

		fprintf(file1, "%.9f  %.9f  \n ", ux_Numeric[i], (y * 0.5) / width);
	}
	fclose(file1);

}
