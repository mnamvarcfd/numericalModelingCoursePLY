#ifndef LATTICE_H
#define LATTICE_H


class Lattice
{
public:
    Lattice(int n);
    double **f_;
    double **f0_;
    double **u_;
    double *rho_;
};

#endif // LATTICE_H
