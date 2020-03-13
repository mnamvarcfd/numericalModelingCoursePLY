
#include <iostream> 
#include <cmath>
#include <time.h> 

// Local includes
#include "Main.h"
#include "F:\PhD Thesis\LBM_Code\AF_LBM_MN\newVersionCode\src\poissulleVerification.h"
#include "F:\PhD Thesis\LBM_Code\AF_LBM_MN\newVersionCode\src\couetteVerification.h"
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

	Grid gridold = Grid(nx_, ny_, nNode_, dx_, xMin_, xMax_, yMin_, yMax_);

	LatticStencil stencil;
	stencil.setStencil();


	//poissulleVerification poissulle;
	//poissulle.setVariablesFullWay(grid.getNx(), grid.getNy());
	//poissulle.setVariablesHalfWay(grid.getNx(), grid.getNy());
	//poissulle.setVariablesOnSite(grid.getNx(), grid.getNy()); 

	couetteVerification poissulle;
	poissulle.setVariables(grid.getNx(), grid.getNy());

	double rhoIn =  poissulle.rhoIn;
	double rhout =  poissulle.rhout;
	double tau = 2./(6* poissulle.nu+1);
	double rho0 = poissulle.rho0;
	double cs = poissulle.cs;
	std::cout << tau << "====================" << rho0 << std::endl;
	
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



	while (/*convergence == false &&*/ it<8000)
	{

		LBM.collision();
		//for (int j = 0; j < ny_; j++) {
		//	for (int i = 0; i < nx_; i++) {
		//		int iNode = j * nx_ + i;
		//		for (int k = 0; k < 9; k++)
		//			gridold.fVal[iNode][k] = grid.fVal[iNode][k];
		//	}
		//}


		LBM.streaming();
		BC.periodic();

		//BC.zouHeLeft(rhoIn);
		//BC.zouHeRigt(rhout);

		//BC.fullWayBouncBackTop();
		//BC.fullWayBouncBackBot();

		//BC.halfWayBouncBackTop(gridold);
		//BC.halfWayBouncBackBot(gridold);
		double uTop = 0.9;
		double uBot = -0.9;

		BC.ZouHeVelTop(uTop, 0.0);
		BC.zouHeVelTopLeftCorner(uTop, 0.0);
		BC.zouHeVelTopRigtCorner(uTop, 0.0);

		BC.ZouHeVelBot(uBot, 0.0);
		BC.zouHeVelBotLeftCorner(uBot, 0.0);
		BC.zouHeVelBotRigtCorner(uBot, 0.0);

		//BC.ZouHeVelLeft(0.0, 0.0);
		//BC.ZouHeVelRigt(0.0, 0.0);


		LBM.calMacroValue();


		deltaTotalMoment = postProc.calMom1Res(grid);


		//if ((deltaTotalMoment) < velocityConvergence)
		//{
		//	std::cout << "Convergence reached : " << deltaTotalMoment << std::endl;
		//	std::cout << "Times : " << it << std::endl;
		//	convergence = true;
		//}
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

