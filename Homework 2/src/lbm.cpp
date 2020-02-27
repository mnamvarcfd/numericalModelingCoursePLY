/*******************************************************************
*	lbm
*	Author : Bruno Blais
*
*	Student Name : Morteza Namvar
*
*	Role :  This is an LBM code to solve Poiseuille Flow in 2D with
*           Periodic boundary conditions
*
********************************************************************/

// C++ Standard includes
#include <iostream> 
#include <cmath>

#include <time.h> 
// Library includes
#include "tinyxml2.h"

// Local includes
#include "parser.h"
#include "domain.h"
#include "lattice.h"
#include "solutionExporterEulerian.h"


void calMacroValue(Domain &domain, double **xi, Lattice lat, double *g, double dt);
void bouncBackBC(Domain domain, Lattice lat);
void periodicBC(Domain domain, Lattice lat);
void streaming(double **Cxy, Domain domain, Lattice lat);
void writeVelocityProfile(std::string fileName, Lattice lat, Domain domain);
void writeResults(Domain domain, Lattice lat);
void collision(Domain domain, Lattice lat, double tau, double dt, double *omega, double **xi,
	double c1, double c2, double c3, double c4, double c5, double *g);






void setXi(double **Xi, double Xi_r)
{
    Xi[0][0] = 0.;
    Xi[1][0] = 0.;

    Xi[0][1] = Xi_r;
    Xi[1][1] = 0.;

    Xi[0][2] = 0.;
    Xi[1][2] = Xi_r;

    Xi[0][3] = -Xi_r;
    Xi[1][3] = 0.;

    Xi[0][4] = 0.;
    Xi[1][4] = -Xi_r;

    Xi[0][5] = Xi_r;
    Xi[1][5] = Xi_r;

    Xi[0][6] = -Xi_r;
    Xi[1][6] = Xi_r;

    Xi[0][7] = -Xi_r;
    Xi[1][7] = -Xi_r;

    Xi[0][8] = Xi_r;
    Xi[1][8] = -Xi_r;
}
void setWeights(double* omega)
{
    omega[0]= 4./9.;

    omega[1]= 1./9.;
    omega[2]= 1./9.;
    omega[3]= 1./9.;
    omega[4]= 1./9.;

    omega[5]= 1./36.;
    omega[6]= 1./36.;
    omega[7]= 1./36.;
    omega[8]= 1./36.;
}

int main() {
    Parser parser;
    parser.parse();

    Domain domain;

    // Domain interpretation

    // Mesh : get the values of nx and ny
    double dx = parser.getDx();

    // Give the values to the domain
    domain.setDx(dx);

    // Set the size of the physical domain
    double xmin = parser.getXmin();
    double xmax = parser.getXmax();
    double ymin = parser.getYmin();
    double ymax = parser.getYmax();
    domain.setXMin(xmin);
    domain.setXMax(xmax);
    domain.setYMin(ymin);
    domain.setYMax(ymax);

    // Finalize creation of the domain
    domain.setFinalize();

    // Once domain has been finalized you can get the value of dx and dy automatically
    int nx = domain.getNx();
    int ny = domain.getNy();
    int ntot = domain.getNTot();

    // Time stepping information from XML file
    double timeEnd = parser.getTimeEnd();
    int outputFrequency = parser.getOutputFrequency();

    // Physical property
    double tau = parser.getTau();
    double mu = parser.getMu();
    double rho0 = parser.getRho0();
    double tau_inv = 1./tau;
    
    // Calculate time step
    double dt= (tau-0.5)*rho0*dx*dx/(3.*mu);
    
    // Calculate Xi_r and cs
    double xi_r =dx/dt;
    double cs = 1./sqrt(3.)*xi_r;
    
    // Microscopic speed
    double **xi;
    xi = new double*[2];
    xi[0] = new double[9];
    xi[1] = new double[9];
    setXi(xi, xi_r);     // Weights 
    double *omega;
    omega = new double[9];
    setWeights(omega);
    
    // Creation of the lattices
    Lattice lat(domain.getNTot());
    
    // Creation of the gravity vector
    double g[2];
    g[0] = 0.;
    g[1] = -9.81;
    
    // Speed at the previous time step (This can be used to monitor convergence)
    double **u_old;
    u_old = new double*[ntot];
    for (int j=0 ; j<domain.getNTot() ; j++)
    {
        u_old[j] = new double[2];
        //u_old[j] = new double[2];  momment by M. namvar
        
        u_old[j][0] = 0.;
        u_old[j][1] = 0.;
    }
    
    // Analytical solution for speed
    double L = xmax-xmin;
    double vmax = (L)*(L)*rho0*g[1]/(8.*mu);
    for (int j=0 ; j<ntot ; j++)
    {
        int k=j%nx;
        double x=xmin +(k+0.5)*dx;
        
        lat.rho_[j] = rho0;
        lat.u_[j][0] = 0.;
        lat.u_[j][1] = vmax*(1.-4.*(x)*(x)/(L*L));
        
        for (int l=0 ; l<9 ; l++)
        {
            lat.f0_[j][l] = 0.;
        }
    }
    
    // Write the analytical solution
    writeLattice(domain,"AnalyticalSolution",0,lat);
	writeVelocityProfile("AnalyticalProfile.plt",lat, domain);

    // Initialization of lattice at rest
    for (int j=0 ; j<domain.getNTot() ; j++)
    {
        lat.rho_[j] = rho0;
        
        lat.u_[j][0] = 0.;
        lat.u_[j][1] = 0.;
        
        for (int k=0 ; k<9 ; k++)
        {
            lat.f0_[j][k] = omega[k]*rho0;
			lat.f_[j][k] = omega[k] * rho0; // m namvar
        }
    }
    
    
    // Pre-determined constants to avoid useless recalculations
    const double c1 = 1./(cs*cs);
    const double c2 = 1./(2.*cs*cs*cs*cs);
    const double c3 = -1./(2.*cs*cs);
    const double c4 = 1./(cs*cs*cs*cs);
    const double c5 = 1.-0.5/tau;

    double *moment1, *sol_analyt;
    moment1 = new double[2];
    sol_analyt = new double[ntot];

    double *fi, *S;
    fi = new double[9];
    S = new double[9];

    bool convergence=false;
    double velocityConvergence=1e-8;
    double t=0;
    int it=0;


	// Variables to use in the equilibrium calculation
	double deltaTotalMoment = 0.0;
	double totalMomentNew =  0.0;
	double totalMomentOld =  0.0;

	// Microscopic unit speed
	double **Cxy;
	Cxy = new double*[2];
	Cxy[0] = new double[9];
	Cxy[1] = new double[9];
	setXi(Cxy, 1.0); 

	// Time loop
	while (t < timeEnd && convergence == false)
	{

		collision(domain, lat, tau, dt, omega, xi, c1, c2, c3, c4, c5, g);

		streaming(Cxy, domain, lat);

		periodicBC(domain, lat);

		bouncBackBC(domain, lat);

		calMacroValue(domain, xi, lat, g, dt);

        // Convergence test for velocity
        totalMomentNew = 0.;
		for (int j = 0; j < ntot; j++)
			totalMomentNew += lat.u_[j][1];

		deltaTotalMoment = abs(totalMomentNew - totalMomentOld)/ abs(totalMomentNew);

	
        if ((deltaTotalMoment)<velocityConvergence)
        {
            std::cout << "Convergence reached : " << deltaTotalMoment <<std::endl;
            std::cout << "Times : "<< t  <<std::endl;
            convergence=true;
        }
		totalMomentOld = totalMomentNew;
		
		if (it%outputFrequency==0)
        {
            std::cout << "Time : " << t  << " - Convergence : "<< (deltaTotalMoment) <<std::endl;
            writeLattice(domain,"Lattice",it,lat);

			writeVelocityProfile("NumericalProfile.plt", lat, domain);
        }
        
        t+=dt;
        it++;
    }

    // Write the last result
    writeLattice(domain,"Lattice",it,lat);

    // Calculate infinity norm of the error
    double erreur=0.;
    for (int j=0 ; j<ntot ; j++)
    {
        int k=j%nx;
        double x=xmin +(k+0.5)*dx;
        double L = xmax-xmin;
        sol_analyt[j] = vmax*(1.-4.*(x)*(x)/(L*L));
        erreur = erreur + (sol_analyt[j]-lat.u_[j][1])*(sol_analyt[j]-lat.u_[j][1]);
    }
    double erreurL2 = sqrt(1./domain.getNTot()*erreur);
    std::cout << "Norm of the error :  "<< erreurL2 <<std::endl;


	writeVelocityProfile("NumericalProfile.plt",lat, domain);
	std::cin >> erreurL2;
    return 0;
}

