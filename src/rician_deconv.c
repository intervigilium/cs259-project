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
#pragma AP pipeline
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

static inline double cubic_approx(const double u, double f, double sigma2)
{
	double r, numer, denom;

	r = u * f;
	r = r / sigma2;
	numer = r * 2.38944 + r * (0.950037 + r);
	denom = 4.65314 + r * (2.57541 + r * (1.48937 + r));
	r = numer / denom;
	return f * r;

}

static inline uint2_t semi_implicit_convergence(double u[M * N * P],
						const double g[M * N * P],
						const double f[M * N * P],
						double dt, double gamma,
						double sigma2)
{
	uint32_t i, j, k;
	double numer, denom;
	double u_last, r;
	double u_stencil_up, u_stencil_center, u_stencil_down;
	double g_stencil_up, g_stencil_center, g_stencil_down;
	double g_left_cache[M], g_right_cache[M], g_in_cache[M], g_out_cache[M];
	double left_mul_cache[M], right_mul_cache[M], in_mul_cache[M],
	    out_mul_cache[M];
	double u_res_cache[M];
	double f_center;

	/* Update u by a semi-implict step */
	for (k = 1; k < P - 1; k++) {
		for (j = 1; j < N - 1; j++) {
			for (i = 0; i < M; i++) {
#pragma AP pipeline
				g_left_cache[i] = G(1, j - 1, k);
				g_right_cache[i] = G(1, j + 1, k);
				g_in_cache[i] = G(1, j, k - 1);
				g_out_cache[i] = G(1, j, k + 1);
				left_mul_cache[i] =
				    g_left_cache[i] * U(1, j - 1, k);
				right_mul_cache[i] =
				    g_right_cache[i] * U(1, j + 1, k);
				in_mul_cache[i] =
				    g_in_cache[i] * U(1, j, k - 1);
				out_mul_cache[i] =
				    g_in_cache[i] * U(1, j, k + 1);
			}
			u_stencil_center = U(0, j, k);
			g_stencil_center = U(0, j, k);
			u_stencil_down = U(1, j, k);
			g_stencil_down = G(1, j, k);
			for (i = 1; i < M - 1; i++) {
				u_last = u_stencil_center;
				u_stencil_up = u_stencil_center;
				g_stencil_up = g_stencil_center;
				u_stencil_center = u_stencil_down;
				g_stencil_center = g_stencil_down;
				u_stencil_down = U_DOWN;
				g_stencil_down = G_DOWN;

				numer =
				    u_stencil_center +
				    dt * (right_mul_cache[i] +
					  left_mul_cache[i] +
					  u_stencil_up * g_stencil_up +
					  u_stencil_down * g_stencil_down +
					  in_mul_cache[i] +
					  out_mul_cache[i] - gamma * cubic_approx_cache[i]);
				denom =
				    1.0 + dt * (g_right_cache[i] +
						g_left_cache[i] +
						g_stencil_down + g_stencil_up +
						g_in_cache[i] + g_out_cache[i]);
				u_stencil_center = numer / denom;

//                              if (fast_fabs(u_last - u_stencil_center) <=
//                                  DENOISE_TOLERANCE) {
//                                      return 1;
//                              }
//                              U_CENTER = u_stencil_center;
				u_res_cache[i] = u_stencil_center;
			}
			for (i = 1; i < M - 1; i++) {
#pragma AP pipeline
				U_CENTER = u_res_cache[i];
			}
		}
	}
}

static inline void semi_implicit_update(double u[M * N * P],
					const double g[M * N * P],
					const double f[M * N * P], double dt,
					double gamma)
{
	uint32_t i, j, k;
	double numer, denom;
	double u_stencil_up, u_stencil_center, u_stencil_down;
	double g_stencil_up, g_stencil_center, g_stencil_down;
	double g_left_cache[M], g_right_cache[M], g_in_cache[M], g_out_cache[M];

	/* Update u by a semi-implict step */
	for (k = 1; k < P - 1; k++) {
		for (j = 1; j < N - 1; j++) {
			for (i = 0; i < M; i++) {
#pragma AP pipeline
				g_left_cache[i] = G(1, j - 1, k);
				g_right_cache[i] = G(1, j + 1, k);
				g_in_cache[i] = G(1, j, k - 1);
				g_out_cache[i] = G(1, j, k + 1);
			}
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
				    dt * (U_RIGHT * g_right_cache[i] +
					  U_LEFT * g_left_cache[i] +
					  u_stencil_up *
					  g_stencil_up +
					  u_stencil_down *
					  g_stencil_down +
					  U_IN * g_in_cache[i] +
					  U_OUT * g_out_cache[i] - gamma * F(i,
									     j,
									     k));
				denom =
				    1.0 + dt * (g_right_cache[i] +
						g_left_cache[i] +
						g_stencil_down + g_stencil_up +
						g_in_cache[i] + g_out_cache[i]);
				u_stencil_center = numer / denom;
				U_CENTER = u_stencil_center;
			}
		}
	}
}

void rician_deconv_deblur(double u[M * N * P], double f[M * N * P],
			  double g[M * N * P], double conv[M * N * P],
			  double Ksigma, double sigma, double lambda)
{
	uint32_t i;
	double sigma2, gamma;
	uint32_t iteration;

	/* Initializations */
	sigma2 = SQR(sigma);
	gamma = lambda / sigma2;

	/* Main gradient descent loop */
	for (iteration = 1; iteration <= DEBLUR_ITERATIONS; iteration++) {
		/* Approximate g = 1/|grad u| */
		gradient(u, g, DEBLUR_EPSILON);
		array_copy(u, conv);

		/* calculate with rational cubic approx */
		for (i = 0; i < M * N * P; i++) {
#pragma AP pipeline
			f[i] = cubic_approx(u[i], f[i], sigma2);
			conv[i] -= f[i];
		}
		gaussian_blur(conv, Ksigma);
		semi_implicit_update(u, g, conv, DEBLUR_DT, gamma);
	}
}

void rician_deconv_denoise(double u[M * N * P], double f[M * N * P],
			   double g[M * N * P], double sigma, double lambda)
{
	double sigma2, gamma;
	uint32_t iteration;
	uint2_t converged;

	/* Initializations */
	sigma2 = SQR(sigma);
	gamma = lambda / sigma2;

	/* Main gradient descent loop */
	for (iteration = 1; iteration <= DENOISE_ITERATIONS; iteration++) {
		/* Approximate g = 1/|grad u| */
		gradient(u, g, DENOISE_EPSILON);
		converged = 1;
		converged =
		    semi_implicit_convergence(u, g, f, DENOISE_DT, gamma,
					      sigma);
//              if (converged) {
//                      return;
//              }
	}
}
