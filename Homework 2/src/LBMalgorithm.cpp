#include "LBMalgorithm.h"

#include "domain.h"
#include "lattice.h"



LBMalgorithm::LBMalgorithm()
{
}


LBMalgorithm::~LBMalgorithm()
{
}


//void streaming(Lattice lat, Domain domain){
//	int nx = domain.getNx();
//	int ny = domain.getNy();
//
//Do j=0,nx-1
//	Do i=imax-1,1,-1
//		f(1,i,j)=f(1,i-1,j)   ! Streaming From Right to Left
//	End Do
//	Do i=1,imax-2
//		f(7,i,j)=f(7,i+1,j)   ! Streaming From Left to Right
//	End Do
//End Do
//
//Do j=jmax-1,1,-1   ! Streaming From Top to Bottom
//
//	Do i=0,imax-1
//		f(2,i,j)=f(2,i,j-1)
//	End Do
//	Do i=imax-1,1,-1
//		f(3,i,j)=f(3,i-1,j-1) !????????????
//	End Do
//
//	Do i=0,imax-2
//		f(8,i,j)=f(8,i+1,j-1)
//	End Do
//
//End Do
//	 
//			 
//Do j=0,jmax-2   ! Streaming From Bottom to Top
//
//	Do i=0,imax-1
//		f(5,i,j)=f(5,i,j+1) 
//	End Do
//
//	Do i=0,imax-2
//		f(6,i,j)=f(6,i+1,j+1)
//	End Do
//
//	Do i=imax-1,1,-1
//		f(4,i,j)=f(4,i-1,j+1)
//	End Do
//
//End Do 
//				
//	
//}



//void collision(Lattice lat, Domain domain, double omega[], double xi[][]) {
//	int nx = domain.getNx();
//	int ny = domain.getNy();
//
//	double u2, uiei, uiei2;
//	double ro;
//	double ux;
//	double uy;
//	double udotu;
//	double feq;
//	double fi;
//
//
//	for (int j = 0; j < ny - 1; j++) {
//		for (int i = 0; i < nx - 1; i++) {
//
//			ro = lat.rho_[j];
//			ux = lat.u_[j][0];
//			uy = lat.u_[j][1];
//			udotu = lat.u_[j][0] * lat.u_[j][0] + lat.u_[j][1] * lat.u_[j][1];
//
//			int iNode = j * ny + i;
//
//			for (int k = 0; k < 9 - 1; k++) {
//				int cx = xi[k][0];
//				int cy = xi[k][1];
//
//				uiei = cx*ux + cy*uy;
//				uiei2 = uiei * uiei;
//				fi = lat.f0_[iNode][k];
//
//				feq = ro * omega[k] * (1 + uiei + uiei2 - udotu * c3);
//
//				//fi = lat.f0_[iNode][k] - (lat.f0_[iNode][k] - feq) / tau;
//
//				//S[n] = c5 * omega[n] * (((x[0][n] - ux) / cs + (x[0][n] * ux) / c4)*g[0] + ((x[1][n] - uy) / cs + (x[1][n] * uy) / c4)*g[1]);
//
//				int streamTo = (j + cy)*nx + i + cx;
//
//				lat.f_[iNode][k] = fi;
//			}
//
//		}
//	}
//}
//
//

