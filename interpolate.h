#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "memory.h"
#include "mkl_lapacke.h"
#include <tricubic.h>

using namespace std;

class Interpolate{
public:
  Interpolate(int, int, int, int, MKL_Complex16 **);
  ~Interpolate();

  void set_method();
  void execute(double *, MKL_Complex16 *);
  void reset_gamma();

  int UseGamma;

private:
  void tricubic_init();
  void tricubic(double *, MKL_Complex16 *);
  void trilinear(double *, MKL_Complex16 *);
  Memory *memory;

  int which;
  int Nx, Ny, Nz, Npt, ndim;
  int flag_reset_gamma, flag_allocated_dfs;

  MKL_Complex16 **data;
  MKL_Complex16 **Dfdx, **Dfdy, **Dfdz, **D2fdxdy, **D2fdxdz, **D2fdydz, **D3fdxdydz;
  double a[64], f[8], dfdx[8], dfdy[8], dfdz[8], d2fdxdy[8], d2fdxdz[8], d2fdydz[8], d3fdxdydz[8];
  int vidx[8];
};

#endif
