
#include <iostream> 
#include <cmath>
#include <time.h> 

// Local includes
#include "Main.h"
#include "poissulleVerificationHW.h";
#include "poissulleVerificationFW.h";
#include "LatticStencil.h"
#include "BoundaryCondition.h"
#include "LBMalgorithm.h"
#include "Grid.h"
#include "postProccessing.h"



int main() {

	int N = 16;
	int nx_ = N; 
	int ny_ = N; 

	int nNode_ = nx_ * ny_; 
	double dx_ = 1. / N; 

	double xMin_ = 0.0; 
	double yMin_ = 0.0;
	 
	double xMax_ = 1.0; 
	double yMax_ = 1.0;


	Grid grid = Grid(nx_, ny_, nNode_, dx_, xMin_, xMax_, yMin_, yMax_);

	LatticStencil stencil;
	stencil.setStencil();


	poissulleVerificationFW poissulle;
	//poissulleVerificationFW poissulle;
	poissulle.setVariables(grid.getNx(), grid.getNy());
	double rhoIn = poissulle.rhoIn;
	double rhout = poissulle.rhout;
	double tau = 2./(6* poissulle.nu+1);
	double rho0 = poissulle.rho0;
	double cs = poissulle.cs;
	int outputFrequency = 100;

	
	LBMalgorithm LBM = LBMalgorithm(grid, tau, cs, stencil);
	LBM.init();

	BoundaryCondition BC = BoundaryCondition(grid);


	postProccessing postProc;


	double velocityConvergence = 1e-8;
	bool convergence = false;
	double t = 0;
	int it = 0;
	double deltaTotalMoment = 0.0;



	while (convergence == false && it<8000)
	{

		LBM.collision();
		//for (int j = 0; j < ny; j++) {
		//	for (int i = 0; i < nx; i++) {
		//		int iNode = j * nx + i;
		//		for (int k = 0; k < 9; k++)
		//			lat.f_[iNode][k] = lat.f0_[iNode][k];
		//	}
		//}


		LBM.streaming();
		//BC.periodic();
	
		BC.zouHeLeft(rhoIn);
		BC.zouHeRigt(rhout);

		BC.fullWayBouncBack();
		//BC.halfWayBouncBack();

		LBM.calMacroValue();


		deltaTotalMoment = postProc.calMom1Res(grid);


		if ((deltaTotalMoment) < velocityConvergence)
		{
			std::cout << "Convergence reached : " << deltaTotalMoment << std::endl;
			std::cout << "Times : " << it << std::endl;
			convergence = true;
		}
		if (it%outputFrequency == 0)
		{
			std::cout << "Time : " << it << " - Convergence : " << (deltaTotalMoment) << std::endl;
			//writeLattice(domain, "Lattice", it, lat);
		}


		it++;
	}



	postProc.writeVectorFiedldVTK(grid);
	poissulle.erreurL2(grid.U);
	poissulle.writeNumericVelcity(grid.U);
	poissulle.writeAnalyticVelcity();


	std::cin >> it;
	return 0;
}

