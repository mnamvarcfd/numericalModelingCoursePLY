#include "lattice.h"

Lattice::Lattice(int n)
{
    f_ = new double*[n];
    f0_ = new double*[n];
    u_ = new double*[n];
    rho_ = new double[n];

    for (int i=0 ; i < n ; i ++)
    {
        f_[i] = new double[9];
        f0_[i] = new double[9];
        u_[i] = new double[2];
    }
}
