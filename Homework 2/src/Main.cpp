
#include <iostream> 
#include <cmath>
#include <time.h> 

// Local includes
#include "Main.h"
#include "poissulleVerification.h"
#include "LatticStencil.h"
#include "BoundaryCondition.h"
#include "LBMalgorithm.h"
#include "Grid.h"
#include "postProccessing.h"



int main() {

	int N = 8;
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


	poissulleVerification poissulle;
	//poissulle.setVariablesFullWay(grid.getNx(), grid.getNy());
	poissulle.setVariablesHalfWay(grid.getNx(), grid.getNy());
	double rhoIn = poissulle.rhoIn;
	double rhout = poissulle.rhout;
	double tau = 1.0   /*2./(6* poissulle.nu+1)*/;
	double rho0 = poissulle.rho0;
	double cs = poissulle.cs;

	
	LBMalgorithm LBM = LBMalgorithm(grid, tau, cs, stencil);
	LBM.init();

	BoundaryCondition BC = BoundaryCondition(grid);


	postProccessing postProc;


	double velocityConvergence = 1e-8;
	bool convergence = false;
	double t = 0;
	int it = 0;
	double deltaTotalMoment = 0.0;
	int outputFrequency = 100;



	while (convergence == false && it<3000)
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

		BC.fullWayBouncBackTop();
		BC.fullWayBouncBackBot();




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
			//postProc.writeResults(grid);
			//std::cin >> it;
		}


		it++;
	}



	postProc.writeVectorFiedldVTK(grid);
	poissulle.erreurL2(grid.U);
	poissulle.writeNumericVelcity(grid.U);
	poissulle.writeAnalyticVelcity();
	postProc.writeScalarFieldVTK(grid);

	std::cin >> it;
	return 0;
}

