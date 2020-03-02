#include "LatticStencil.h"

LatticStencil::LatticStencil()
{
}

void LatticStencil::setStencil()
{
	cx[0] = 0.0;
	cx[1] = 0.0; cx[2] = 1.0; cx[3] = 1.0; cx[4] = 1.0;
	cx[5] = 0.0; cx[6] = -1.0; cx[7] = -1.0; cx[8] = -1.0;

	cy[0] = 0.0;
	cy[1] = 1.0; cy[2] = 1.0; cy[3] = 0.0; cy[4] = -1.0;
	cy[5] = -1.0; cy[6] = -1.0; cy[7] = 0.0; cy[8] = 1.0;

	wt[0] = 4. / 9;
	wt[1] = 1. / 9; wt[2] = 1. / 36; wt[3] = 1. / 9; wt[4] = 1. / 36;
	wt[5] = 1. / 9; wt[6] = 1. / 36; wt[7] = 1. / 9; wt[8] = 1. / 36;
}

LatticStencil::~LatticStencil()
{
}
