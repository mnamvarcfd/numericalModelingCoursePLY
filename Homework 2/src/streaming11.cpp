#include "lattice.h"
#include "domain.h"
#include <iostream>


void streaming(int nx, int ny, int j, double cx, double cy, Lattice &lat, double fi, int k) {
	int ii = j % nx;
	int jj = j / nx;

	int pop = (jj + cy)*nx + ii + cx;

	//if(k==1) std::cout << pop << std::endl;
	lat.f_[pop][k] = fi;
}


void streamingBottomPeriodicity(int nx, int ny, int j, double cx, double cy, Lattice &lat, double fi, int k) {
	int ii = j % nx;
	int jj= j / nx;

	int pop = (jj + cy)*nx + ii + cx;
	if(k==7 || k==4 || k==8) return;

	lat.f_[pop][k] = fi;
}

void streamingTopPeriodicity(int nx, int ny, int j, double cx, double cy, Lattice &lat, double fi, int k) {
	int ii = j % nx;
	int jj = j / nx;
	int pop = (jj + cy)*nx + ii + cx;

	if (k == 6 || k == 2 || k == 5) return /*pop = ii*/;

	lat.f_[pop][k] = fi;
}



void streamingLeftWall(int nx, int ny, int j, double cx, double cy, Lattice &lat, double fi, int k) {
	int ii = j % nx;
	int jj = j / nx;
	if (k == 6 || k == 3 || k == 7) return;

	int pop = (jj + cy)*nx + ii + cx;
	//if (pop == 13) std::cout << fi << "---"<< k <<std::endl;

	lat.f_[pop][k] = fi;
}

void streamingRigtWall(int nx, int ny, int j, double cx, double cy, Lattice &lat, double fi, int k) {
	int ii = j % nx;
	int jj;
	jj = j / nx;
	if (k == 5 || k == 1 || k == 8) return;

	int pop = (jj + cy)*nx + ii + cx;
	//if (j == 4) std::cout << pop << "---"<< k <<std::endl;

	lat.f_[pop][k] = fi;
}



void streamingBottomLeftWall(int nx, int ny, int j, double cx, double cy, Lattice &lat, double fi, int k) {
	int ii = j % nx;
	int jj;
	jj = j / nx;
	//std::cout <<  "-========-" << k << std::endl;
	if (k != 1 && k != 2 && k != 5) return;

	int pop = (jj + cy)*nx + ii + cx;

	//std::cout << pop << "---"<< k <<std::endl;
	lat.f_[pop][k] = fi;
}

void streamingTopLeftWall(int nx, int ny, int j, double cx, double cy, Lattice &lat, double fi, int k) {
	int ii = j % nx;
	int jj;
	jj = j / nx;
	if (k != 1 && k != 4 && k != 8) return;

	int pop = (jj + cy)*nx + ii + cx;
	//if (j == 4) std::cout << pop << "---"<< k <<std::endl;

	lat.f_[pop][k] = fi;
}

void streamingTopRigtWall(int nx, int ny, int j, double cx, double cy, Lattice &lat, double fi, int k) {
	int ii = j % nx;
	int jj;
	jj = j / nx;
	if (k != 3 && k != 4 && k != 7) return;

	int pop = (jj + cy)*nx + ii + cx;
	//if (j == 4) std::cout << pop << "---"<< k <<std::endl;

	lat.f_[pop][k] = fi;
}

void streamingBottomRigtWall(int nx, int ny, int j, double cx, double cy, Lattice &lat, double fi, int k) {
	int ii = j % nx;
	int jj;
	jj = j / nx;
	if (k != 3 && k != 6 && k != 2) return;

	int pop = (jj + cy)*nx + ii + cx;
	//if (j == 4) std::cout << pop << "---"<< k <<std::endl;

	lat.f_[pop][k] = fi;
}