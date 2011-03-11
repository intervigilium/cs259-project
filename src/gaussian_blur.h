#include "common.h"

#define GAUSSIAN_NUMSTEPS 3

#define U(a,b,c) (u[a+b*N+c*M*N])

void gaussian_blur(double u[M * N * P], double Ksigma);
