#include "common.h"

void rician_deconv(double u[M * N * P], const double f[M * N * P],
		   double g[M * N * P], double conv[M * N * P], double Ksigma,
		   double sigma, double lambda);
