#ifndef _RICIAN_DECONV_H
#define _RICIAN_DECONV_H

#include "common.h"

#define DEBLUR_EPSILON 0.0001
#define DEBLUR_DT 1.0E-10
#define DEBLUR_ITERATIONS 10

#define DENOISE_EPSILON 5.0
#define DENOISE_DT 1.0E-20
#define DENOISE_ITERATIONS 500

#define U(a,b,c) (u[a+b*N+c*M*N])
#define G(a,b,c) (g[a+b*N+c*M*N])

#define U_CENTER U(i,j,k)
#define U_LEFT U(i,j-1,k)
#define U_RIGHT U(i,j+1,k)
#define U_UP U(i-1,j,k)
#define U_DOWN U(i+1,j,k)
#define U_IN U(i,j,k-1)
#define U_OUT U(i,j,k+1)

#define G_CENTER G(i,j,k)
#define G_LEFT G(i,j-1,k)
#define G_RIGHT G(i,j+1,k)
#define G_UP G(i-1,j,k)
#define G_DOWN G(i+1,j,k)
#define G_IN G(i,j,k-1)
#define G_OUT G(i,j,k+1)

void rician_deconv(double u[M * N * P], const double f[M * N * P],
		   double g[M * N * P], double conv[M * N * P], double Ksigma,
		   double sigma, double lambda, uint2_t deblur);

#endif
