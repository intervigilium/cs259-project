#include "rician_deconv.h"
#include "gaussian_blur.h"

static inline void gradient(const double u[M * N * P], double g[M * N * P],
			    double epsilon)
{
	uint32_t i, j, k;
	double u_stencil_up, u_stencil_center, u_stencil_down;
	double numer, denom;
	
	/* approximate g = 1/|grad u| */
	for (k = 1; k < P - 1; k++) {
		for (j = 1; j < N - 1; j++) {
			u_stencil_center = U(0, j, k);
			u_stencil_down = U(1, j, k);
			for (i = 1; i < M - 1; i++) {
				u_stencil_up = u_stencil_center;
				u_stencil_center = u_stencil_down;
				u_stencil_down = U_DOWN;
				denom = q3_sqrt(epsilon +
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
						SQR(u_stencil_center - U_OUT));
				G_CENTER = 1.0 / denom;
			}
		}
	}
}

static inline void cubic_approx(const double u[M * N * P], double f[M * N * P],
				double sigma2)
{
	uint32_t i;
	double r, numer, denom;
	for (i = 0; i < M * N * P; i++) {
#pragma pipeline
		r = u[i] * f[i] / sigma2;
		numer = r * 2.38944 + r * (0.950037 + r);
		denom = 4.65314 + r * (2.57541 + r * (1.48937 + r));
		f[i] *= r;
	}
}

static inline void semi_implicit_update(double u[M * N * P],
					const double g[M * N * P],
					const double f[M * N * P], double dt, double gamma)
{
	uint32_t i, j, k;
	double numer, denom;
	double u_stencil_up, u_stencil_center, u_stencil_down;
	double g_stencil_up, g_stencil_center, g_stencil_down;

	/* Update u by a semi-implict step */
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
				    dt * (U_RIGHT * G_RIGHT +
					  U_LEFT * G_LEFT +
					  U_RIGHT * G_RIGHT +
					  u_stencil_up *
					  g_stencil_up +
					  u_stencil_down *
					  g_stencil_down +
					  U_IN * G_IN +
					  U_OUT * G_OUT -
					  gamma * f[i + j * M + k * M * N]);
				denom =
				    1.0 + dt * (G_RIGHT + G_LEFT +
						g_stencil_down +
						g_stencil_up + G_IN + G_OUT);
				U_CENTER = numer / denom;
			}
		}
	}
}

void rician_deconv(double u[M * N * P], double f[M * N * P],
		   double g[M * N * P], double conv[M * N * P],
		   double Ksigma, double sigma, double lambda, uint2_t deblur)
{
	double sigma2, gamma, r;
	double numer, denom;
	double u_stencil_up, u_stencil_center, u_stencil_down;
	double g_stencil_up, g_stencil_center, g_stencil_down;
	double epsilon, dt;
	uint32_t i, j, k;
	uint32_t max_iterations, iteration;

	/* Initializations */
	sigma2 = SQR(sigma);
	gamma = lambda / sigma2;
	if (deblur) {
		epsilon = DEBLUR_EPSILON;
		dt = DEBLUR_DT;
		max_iterations = DEBLUR_ITERATIONS;
	} else {
		epsilon = DENOISE_EPSILON;
		dt = DENOISE_DT;
		max_iterations = DENOISE_ITERATIONS;
	}

	/* Main gradient descent loop */
	for (iteration = 1; iteration <= max_iterations; iteration++) {
		/* Approximate g = 1/|grad u| */
		gradient(u, g, epsilon);

		/* calculate with rational cubic approx */
		cubic_approx(u, f, sigma2);

		if (deblur) {
			array_copy(u, conv);
			for (i = 0; i < M * N * P; i++) {
#pragma AP pipeline
				conv[i] -= f[i];
			}
			gaussian_blur(conv, Ksigma);
			semi_implicit_update(u, g, conv, dt, gamma);
		} else {
			semi_implicit_update(u, g, f, dt, gamma);
		}
	}
}
