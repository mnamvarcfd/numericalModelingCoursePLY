/*****************************************
*	Class Domain
*	Author : Bruno Blais
*
*	Role : This class stores the information about the geometry and the domain
*
******************************************/


#ifndef DOMAIN_H
#define DOMAIN_H
#include <iostream>

class  Domain
{
 public:
    Domain();          // Public constructor
    ~Domain(); 	       // Public destructor
	int nx;
	int ny;
	int ntot;
    // Accessors
    int getNx()   {return nx_;};
    int getNy()   {return ny_;};
    int getNTot() {return nTot_;};
    int getLx()   {return Lx_;};
    int getLy()   {return Ly_;};

    double getXMin() {return xMin_;};
    double getXMax() {return xMax_;};
    double getYMin() {return yMin_;};
    double getYMax() {return yMax_;};
    double getDX()   {return dx_;};

    // Setters
    void setXMin(double l)
    {
        xMin_= l;
    }
    void setXMax(double l)
    {
        xMax_=l;
    }
    void setYMin(double l)
    {
        yMin_= l;
    }
    void setYMax(double l)
    {
        yMax_=l;
    }

    void setDx(double dx)
    {
        dx_=dx;
    }

    void setFinalize()
    {
        Lx_=xMax_-xMin_;
        Ly_=yMax_-yMin_;

        nx_= int ( (Lx_+1e-6*Lx_) / dx_);
        ny_= int ( (Ly_+1e-6*Ly_) / dx_);
		std::cout << "nx " << nx_ << std::endl;
        std::cout << "ny " << ny_ << std::endl;
        std::cout << "dx " << dx_ << std::endl;
        nTot_=nx_*ny_;
    }

 private:
        int nx_;
        int ny_;
        int nTot_;

        double xMin_;
        double xMax_;
        double yMin_;
        double yMax_;
        double Lx_;
        double Ly_;
        double dx_;
};
#endif
