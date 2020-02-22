#include "domain.h"
#include "lattice.h"

double equilibrum(double ux, double uy, double udotu, double ro,
	double *omega, int k, double cx, double cy, double c1, double c2, double c3) {

	double uiei;
	double uiei2;
	double Feq;
	double fi;
	double Si;

	uiei = cx * ux + cy * uy;
	uiei2 = uiei * uiei;
	Feq = ro * omega[k] * (1 + c1 * uiei + c2 * uiei2 + c3 * udotu);

	return Feq;
}


double source(double ux, double uy, double *omega, int k, double cx, double cy, 
	double c1,	double c4, double c5, double *g) {

	double Si = c5 * omega[k] * ((cy - uy)*c1 + (cx * ux + cy * uy)*cy*c4) * g[1];

	return Si;
}

double applyCollision(int iNode, Lattice lat, int k, double Feq, double Si, double tau, double dt) {

	double fi = lat.f0_[iNode][k];
	fi = fi - (fi - Feq) / tau + Si * dt;

	return fi;
}




//double collision(int iNode, Lattice lat, double ux, double uy, double udotu, double ro,
//	double *omega, int k, double cx, double cy, double c1, double c2, double c3,
//	double c4, double c5, double cs, double tau, double *g, double dt) {
//
//	double uiei;
//	double uiei2;
//	double Feq;
//	double fi;
//	double Si;
//
//	uiei = cx * ux + cy * uy;
//	uiei2 = uiei * uiei;
//	Feq = ro * omega[k] * (1 + c1 * uiei + c2 * uiei2 + c3 * udotu);
//
//	//Feq = ro * omega[k] * (1 + 3 * uiei + 4.5 * uiei2 - 1.5 * udotu);
//
//	Si = c5 * omega[k] * ((cy - uy)*c1 + (cx * ux + cy * uy)*cy*c4) * g[1];
//	//std::cout << " time is  " << tau << "---" << k << std::endl;
//	//Si = 0.0;
//	fi = lat.f0_[iNode][k];
//	fi = fi - (fi - Feq) / tau + Si * dt;
//
//	return fi;
//}
//
//double collisionWall(int iNode, Lattice lat, int k) {
//	return lat.f0_[iNode][k];
//}


