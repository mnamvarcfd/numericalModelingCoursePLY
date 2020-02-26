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
void FWBB(Domain domain, Lattice lat);
void periodic(Domain domain, Lattice lat);

int bottomStreamTo(int iNode, int k, double cx, double cy, Domain domain);
int topStreamTo(int iNode, int k, double cx, double cy, Domain domain);
int streamTo(int iNode, int k, double cx, double cy, Domain domain);
void streamingNew(int pop, Lattice lat, int k, double fi);

double equilibrum(double ux, double uy, double udotu, double ro,
	double *omega, int k, double cx, double cy, double c1, double c2, double c3);
double source(double ux, double uy, double *omega, int k, double cx, double cy,
	double c1, double c4, double c5, double *g);
double applyCollision(int iNode, Lattice lat, int k, double Feq, double Si, double tau, double dt);


void zouHeTop(Domain domain, Lattice lat, double rho);
void zouHeBottom(Domain domain, Lattice lat, double rho);


void writeVelocityProfile(std::string fileName, Lattice lat, Domain domain);







//void streaming(int nx, int ny, int j, double cx, double cy, Lattice &lat, double fi, int k);
//void streamingBottomPeriodicity(int nx, int ny, int j, double cx, double cy, Lattice &lat, double fi, int k);
//void streamingTopPeriodicity(int nx, int ny, int j, double cx, double cy, Lattice &lat, double fi, int k);
//void streamingLeftWall(int nx, int ny, int j, double cx, double cy, Lattice &lat, double fi, int k);
//void streamingRigtWall(int nx, int ny, int j, double cx, double cy, Lattice &lat, double fi, int k);
//
//void streamingBottomLeftWall(int nx, int ny, int j, double cx, double cy, Lattice &lat, double fi, int k);
//void streamingTopLeftWall(int nx, int ny, int j, double cx, double cy, Lattice &lat, double fi, int k);
//void streamingTopRigtWall(int nx, int ny, int j, double cx, double cy, Lattice &lat, double fi, int k);
//void streamingBottomRigtWall(int nx, int ny, int j, double cx, double cy, Lattice &lat, double fi, int k);
//double collision(int iNode, Lattice lat, double ux, double uy, double udotu, double ro,
//	double *omega, int k, double cx, double cy, double c1, double c2, double c3, 
//	double c4, double c5, double cs, double tau, double *g, double dt);
//double collisionWall(int iNode, Lattice lat, int k);
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

int main()
{
    Parser parser;
    parser.parse();

    Domain domain;

    // Domain interpretation

    // Mesh : get the values of nx and ny
    double dx = parser.getDx() /*1/9.*/ ;

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
    double xi_r =/*1.0*/  dx/dt; //momment by M.namvar
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

    double *Feq, *S;
    Feq = new double[9];
    S = new double[9];

    bool convergence=false;
    double velocityConvergence=1e-8;
    double t=0;
    int it=0;


	////clock_t t1, t2;
	////t1 = clock();

	// Variables to use in the equilibrium calculation
	double ro,ux,uy,udotu,fi,cx,cy,feq,si;

	// Time loop
    while(t<timeEnd && convergence==false /*&& it<8000*/)
    {
        // Space loop
        for (int j=0 ; j< ntot; j++)
        {
			ro = lat.rho_[j];
			ux = lat.u_[j][0];
			uy = lat.u_[j][1] /*- 0.1*/ ;
			udotu = ux * ux + uy * uy;

		  // Bottom left
			if (j == 0)
			{
				for (int k = 0; k < 9; k++) {
					cx = xi[0][k];
					cy = xi[1][k];

					feq = equilibrum(ux, uy, udotu, ro, omega, k, cx, cy, c1, c2, c3);
					si = source(ux, uy, omega, k, cx, cy, c1, c4, c5, g);
					fi = applyCollision(j, lat, k, feq, si, tau, dt);

					//fi = lat.f0_[j][k];  
					  
					if (k != 1 && k != 2 && k != 5) continue;

					int pop = streamTo(j, k, cx, cy, domain);
					streamingNew(pop, lat, k, fi);
				}
			}    
			// Bottom right
			else if (j == nx - 1)
			{
				for (int k = 0; k < 9; k++) {
					cx = xi[0][k];
					cy = xi[1][k];

					feq = equilibrum(ux, uy, udotu, ro, omega, k, cx, cy, c1, c2, c3);
					si = source(ux, uy, omega, k, cx, cy, c1, c4, c5, g);
					fi = applyCollision(j, lat, k, feq, si, tau, dt);

					//fi = lat.f0_[j][k];

					if (k != 3 && k != 6 && k != 2) continue;

					int pop = streamTo(j, k, cx, cy, domain);
					streamingNew(pop, lat, k, fi);
				}
			} 
			// Top left
			else if (j == ntot - nx)
			{
				
				for (int k = 0; k < 9; k++) {
					cx = xi[0][k];
					cy = xi[1][k];

					feq = equilibrum(ux, uy, udotu, ro, omega, k, cx, cy, c1, c2, c3);
					si = source(ux, uy, omega, k, cx, cy, c1, c4, c5, g);
					fi = applyCollision(j, lat, k, feq, si, tau, dt);

					//fi = lat.f0_[j][k];

					if (k != 1 && k != 4 && k != 8) continue;

					int pop = streamTo(j, k, cx, cy, domain);
					streamingNew(pop, lat, k, fi);
				}
			} 
			// Top right
			else if (j == ntot - 1) 
			{
				for (int k = 0; k < 9; k++) {
					cx = xi[0][k];
					cy = xi[1][k];

					feq = equilibrum(ux, uy, udotu, ro, omega, k, cx, cy, c1, c2, c3);
					si = source(ux, uy, omega, k, cx, cy, c1, c4, c5, g);
					fi = applyCollision(j, lat, k, feq, si, tau, dt);

					//fi = lat.f0_[j][k];

					if (k != 3 && k != 4 && k != 7) continue;

					int pop = streamTo(j, k, cx, cy, domain);
					streamingNew(pop, lat, k, fi);
				}
			}
			//Bottom Periodicity
			else if (j > 0 && j < (nx - 1)) 
			{
				for (int k = 0; k < 9; k++) {
					cx = xi[0][k];
					cy = xi[1][k];

					feq = equilibrum(ux, uy, udotu, ro, omega, k, cx, cy, c1, c2, c3);
					si = source(ux, uy, omega, k, cx, cy, c1, c4, c5, g);
					fi = applyCollision(j, lat, k, feq, si, tau, dt);

					//fi = lat.f0_[j][k];

					if (k == 7 || k == 4 || k == 8) continue;

					int pop = streamTo(j, k, cx, cy, domain);
					streamingNew(pop, lat, k, fi);
				}

				//periodicBottom(domain, lat, 0.9);
			}
			// Top Periodicity
			else if (j > (ntot - nx) && j < (ntot - 1)) 
			{
				for (int k = 0; k < 9; k++) {
					cx = xi[0][k];
					cy = xi[1][k];

					feq = equilibrum(ux, uy, udotu, ro, omega, k, cx, cy, c1, c2, c3);
					si = source(ux, uy, omega, k, cx, cy, c1, c4, c5, g);
					fi = applyCollision(j, lat, k, feq, si, tau, dt);

					//fi = lat.f0_[j][k];

					if (k == 6 || k == 2 || k == 5) continue;

					int pop = streamTo(j, k, cx, cy, domain);
					streamingNew(pop, lat, k, fi);

				}


				//periodicTop(domain, lat, 1.1);
			}
			// Right
			else if ((j + 1) % nx == 0) 
			{
				for (int k = 0; k < 9; k++) {
					cx = xi[0][k];
					cy = xi[1][k];

					feq = equilibrum(ux, uy, udotu, ro, omega, k, cx, cy, c1, c2, c3);
					si = source(ux, uy, omega, k, cx, cy, c1, c4, c5, g);
					fi = applyCollision(j, lat, k, feq, si, tau, dt);

					//fi = lat.f0_[j][k];

					if (k == 5 || k == 1 || k == 8) continue;
			
					int pop = streamTo(j, k, cx, cy, domain);
					streamingNew(pop, lat, k, fi);
				}
			}
			// Left
			else if (j%nx == 0) 
			{
				for (int k = 0; k < 9; k++) {
					cx = xi[0][k];
					cy = xi[1][k];

					feq = equilibrum(ux, uy, udotu, ro, omega, k, cx, cy, c1, c2, c3);
					si = source(ux, uy, omega, k, cx, cy, c1, c4, c5, g);
					fi = applyCollision(j, lat, k, feq, si, tau, dt);

					//fi = lat.f0_[j][k];

					if (k == 6 || k == 3 || k == 7) continue;

					int pop = streamTo(j, k, cx, cy, domain);
					streamingNew(pop, lat, k, fi);
				}
			}
			// non boundary
			else
			{
				for (int k = 0; k < 9; k++) {
					cx = xi[0][k];
					cy = xi[1][k];

					feq = equilibrum(ux, uy, udotu, ro, omega, k, cx, cy, c1, c2, c3);
					si = source(ux, uy, omega, k, cx, cy, c1, c4, c5, g);
					fi = applyCollision(j, lat, k, feq, si, tau, dt);

					int pop = streamTo(j, k, cx, cy, domain);
					streamingNew(pop, lat, k, fi);
				}
			}

        }

		lat.f_[ntot - nx][8] = lat.f_[1][8];
		lat.f_[0][5] = lat.f_[1][5];

		FWBB(domain, lat);
		periodic(domain, lat);
		//system("pause");

		//zouHeBottom(domain, lat, 0.9);
		//zouHeTop(domain, lat, 1.1);


		calMacroValue(domain, xi, lat, g, dt);


		//writeResults(domain, lat);
		//writeLattice(domain, "Lattice", it, lat);
		//std::cin >> fi;
		//break;

        // Convergence test for velocity
        double deltaUmax = 0.;
        for (int j=0 ; j<domain.getNTot() ; j++)
        {
            double deltaU = sqrt( (lat.u_[j][1] - u_old[j][1])*(lat.u_[j][1] - u_old[j][1]) );
            deltaUmax=std::max(deltaU,deltaUmax);
        }
        
        if ((deltaUmax/dt)<velocityConvergence)
        {
            std::cout << "Convergence reached : " << (deltaUmax/dt) <<std::endl;
            std::cout << "Times : "<< t  <<std::endl;
            convergence=true;
        }
        else if (it%outputFrequency==0)
        {
            std::cout << "Time : " << t  << " - Convergence : "<< (deltaUmax/dt) <<std::endl;
            writeLattice(domain,"Lattice",it,lat);

			writeVelocityProfile("NumericalProfile.plt", lat, domain);
        }

        // This copies the content of lattice n to lattice n-1
        // f0 = f
        // u0 = u
        for (int j=0 ; j< ntot; j++)
        {
            for (int n=0 ; n<9 ; n++)
            {
                lat.f0_[j][n] = lat.f_[j][n];
            }
            
            u_old[j][0] = lat.u_[j][0];
            u_old[j][1] = lat.u_[j][1];
        }
        
        t+=dt;
        it++;
    }

	////t2 = clock();
	////std::cout << " time is  " << (float)(t2 - t1) / CLOCKS_PER_SEC << std::endl;

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
