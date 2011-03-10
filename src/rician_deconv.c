#include "rician_deconv.h"

void rician_deconv(double u[M * N * P], const double f[M * N * P],
		   double g[M * N * P], double conv[M * N * P],
		   double Ksigma, double sigma, double lambda)
{
#pragma AP interface ap_bus port=f pipeline
#pragma AP interface ap_bus port=u pipeline
#pragma AP interface ap_memory port=g pipeline
#pragma AP interface ap_memory port=conv pipeline

	double sigma2, gamma, r;
	double numer, denom;
	double u_stencil_up, u_stencil_center, u_stencil_down;
	double g_stencil_up, g_stencil_center, g_stencil_down;
	int i, j, k;
	int iteration;

	/* Initializations */
	sigma2 = SQR(sigma);
	gamma = lambda / sigma2;

	/* Main gradient descent loop */
	for (iteration = 1; iteration <= MAX_ITERATIONS; iteration++) {
		/* parallelize/pipeline this, no data deps */
		/* Approximate g = 1/|grad u| */
		for (k = 1; k < P - 1; k++) {
			for (j = 1; j < N - 1; j++) {
				u_stencil_center = U(0, j, k);
				u_stencil_down = U(1, j, k);
				for (i = 1; i < M - 1; i++) {
					u_stencil_up = u_stencil_center;
					u_stencil_center = u_stencil_down;
					u_stencil_down = U_DOWN;
					denom =
					    q3_sqrt(EPSILON +
						    SQR(u_stencil_center -
							U_RIGHT) +
						    SQR(u_stencil_center -
							U_LEFT) +
						    SQR(u_stencil_center -
							u_stencil_up) +
						    SQR(u_stencil_center -
							u_stencil_down) +
						    SQR(u_stencil_center -
							U_IN) +
						    SQR(u_stencil_center -
							U_OUT));
					G_CENTER = 1.0 / denom;
				}
			}
		}
		array_copy(u, conv);
		/* parallelize/pipeline this, no data deps */
		for (i = 0; i < M; i++) {
#pragma AP pipeline
			r = conv[i] * f[i] / sigma2;
			numer = r * 2.38944 + r * (0.950037 + r);
			denom =
			    4.65314 + r * (2.57541 +
					   r * (2.57541 + r * (1.48937 + r)));
			conv[i] -= f[i] * r;
		}

		gaussian_blur(conv, Ksigma);
		/* Update u by a semi-implict step */
		/* pipeline? data deps due to u[i][j][k] writeback */
		for (k = 1; k < P - 1; k++) {
			for (j = 1; j < N - 1; j++) {
				u_stencil_center = U(0, j, k);
				g_stencil_center = U(0, j, k);
				u_stencil_down = U(1, j, k);
				g_stencil_down = G(1, j, k);
				for (i = 1; i < M - 1; i++) {
					u_stencil_up = u_stencil_center;
					g_stencil_up = g_stencil_center;
					u_stencil_center = u_stencil_down;
					g_stencil_center = g_stencil_down;
					u_stencil_down = U_DOWN;
					g_stencil_down = G_DOWN;

					numer =
					    u_stencil_center +
					    DT * (U_RIGHT * G_RIGHT +
						  U_LEFT * G_LEFT +
						  U_RIGHT * G_RIGHT +
						  u_stencil_up *
						  g_stencil_up +
						  u_stencil_down *
						  g_stencil_down +
						  U_IN * G_IN +
						  U_OUT * G_OUT -
						  gamma * conv[i + j * M +
							       k * N * M]);
					denom =
					    1.0 + DT * (G_RIGHT + G_LEFT +
							g_stencil_down +
							g_stencil_up + G_IN +
							G_OUT);
					U_CENTER = numer / denom;
				}
			}
		}
	}
}
