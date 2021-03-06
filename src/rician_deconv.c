#include "rician_deconv.h"
#include "gaussian_blur.h"

static inline void gradient(const double u[M][N][P], double g[M][N][P],
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

static inline double cubic_approx(double u, double f, double sigma2)
{
	double r, numer, denom;

	r = u * f;
	r = r / sigma2;
	numer = r * 2.38944 + r * (0.950037 + r);
	denom = 4.65314 + r * (2.57541 + r * (1.48937 + r));
	r = numer / denom;
	return f * r;

}

static inline uint2_t semi_implicit_convergence(double u[M][N][P],
						const double g[M][N][P],
						const double f[M][N][P],
						double dt, double gamma,
						double sigma2)
{
	uint32_t i, j, k;
	double u_last;
	double u_stencil_center;
	double g_in, g_out;
	double g_right_cache[M], g_center_cache[M], g_left_cache[M];
	double u_right_cache[M], u_center_cache[M], u_left_cache[M];
	double u_stencil_cache[M], g_stencil_cache[M];
	double up_cache[M];
	/* caches j and i direction accesses */

	/* Update u by a semi-implict step */
	for (k = 1; k < P - 1; k++) {
		for (i = 0; i < M; i++) {
#pragma AP pipeline
			/* load up initial j+1 caches */
			u_center_cache[i] = U(i, 0, k);
			u_right_cache[i] = U(i, 1, k);

			g_center_cache[i] = G(i, 0, k);
			g_right_cache[i] = G(i, 1, k);
		}
		for (j = 1; j < N - 1; j++) {
			for (i = 1; i < M - 1; i++) {
#pragma AP pipeline
				/* load up next set of j+1 caches */
				u_left_cache[i] = u_center_cache[i];
				u_center_cache[i] = u_right_cache[i];
				u_right_cache[i] = U(i, j + 1, k);

				g_left_cache[i] = g_center_cache[i];
				g_center_cache[i] = g_right_cache[i];
				g_right_cache[i] = G(i, j + 1, k);

				g_in = G(i, j, k - 1);
				g_out = G(i, j, k + 1);

				/* precompute stencil at i,j,k for u */
				u_stencil_cache[i] =
				    g_left_cache[i] * u_left_cache[i];
				u_stencil_cache[i] +=
				    g_right_cache[i] * u_right_cache[i];
				u_stencil_cache[i] += g_in * U(i, j, k - 1);
				u_stencil_cache[i] += g_out * U(i, j, k + 1);
				/* stencil for u_down * g_down */
				u_stencil_cache[i] +=
				    g_center_cache[i + 1] * u_center_cache[i +
									   1];
				/* add cubic approx to stencil data */
				u_stencil_cache[i] +=
				    cubic_approx(u_center_cache[i], F(i, j, k),
						 sigma2) * gamma;
				u_stencil_cache[i] *= dt;
				u_stencil_cache[i] += u_center_cache[i];

				/* precompute stencil at i,j,k for g */
				g_stencil_cache[i] = g_left_cache[i];
				g_stencil_cache[i] += g_right_cache[i];
				g_stencil_cache[i] += g_in;
				g_stencil_cache[i] += g_out;
				/* stencil for g_up, g_down */
				g_stencil_cache[i] += g_center_cache[i - 1];
				g_stencil_cache[i] += g_center_cache[i + 1];
				g_stencil_cache[i] += gamma;
				g_stencil_cache[i] *= dt;
				g_stencil_cache[i] += 1.0;
			}
			for (i = 1; i < M - 1; i++) {
#pragma AP pipeline
				/* no guarantee u_center_cache[i-1] will get the right value if inlined into previous loop */
				up_cache[i] =
				    u_center_cache[i - 1] * g_center_cache[i -
									   1] *
				    dt;
			}
			for (i = 1; i < M - 1; i++) {
#pragma AP pipeline
				u_stencil_center =
				    (u_stencil_cache[i] +
				     up_cache[i]) / g_stencil_cache[i];

//                              if (fast_fabs(u_last - u_stencil_center) <=
//                                  DENOISE_TOLERANCE) {
//                                      return 1;
//                              }

				/* save result of this i,j,k */
				u_last = u_stencil_center;
				u_center_cache[i] = u_stencil_center;
				U_CENTER = u_stencil_center;
			}
		}
	}
	return 0;
}

static inline void semi_implicit_update(double u[M][N][P],
					const double g[M][N][P],
					const double f[M][N][P], double dt,
					double gamma)
{
	uint32_t i, j, k;
	double u_stencil_center;
	double g_in, g_out;
	double g_right_cache[M], g_center_cache[M], g_left_cache[M];
	double u_right_cache[M], u_center_cache[M], u_left_cache[M];
	double u_stencil_cache[M], g_stencil_cache[M];
	double up_cache[M];

	/* Update u by a semi-implict step */
	for (k = 1; k < P - 1; k++) {
		for (i = 0; i < M; i++) {
#pragma AP pipeline
			/* load up initial j+1 caches */
			u_center_cache[i] = U(i, 0, k);
			u_right_cache[i] = U(i, 1, k);

			g_center_cache[i] = G(i, 0, k);
			g_right_cache[i] = G(i, 1, k);
		}
		for (j = 1; j < N - 1; j++) {
			for (i = 0; i < M; i++) {
#pragma AP pipeline
				/* load up next set of j+1 caches */
				u_left_cache[i] = u_center_cache[i];
				u_center_cache[i] = u_right_cache[i];
				u_right_cache[i] = U(i, j + 1, k);

				g_left_cache[i] = g_center_cache[i];
				g_center_cache[i] = g_right_cache[i];
				g_right_cache[i] = G(i, j + 1, k);

				g_in = G(i, j, k - 1);
				g_out = G(i, j, k + 1);

				/* precompute stencil at i,j,k for u */
				u_stencil_cache[i] =
				    g_left_cache[i] * u_left_cache[i];
				u_stencil_cache[i] +=
				    g_right_cache[i] * u_right_cache[i];
				u_stencil_cache[i] += g_in * U(i, j, k - 1);
				u_stencil_cache[i] += g_out * U(i, j, k + 1);
				/* stencil for u_down * g_down */
				u_stencil_cache[i] +=
				    g_center_cache[i + 1] * u_center_cache[i +
									   1];

				/* subtract conv data to stencil data */
				u_stencil_cache[i] -= F(i, j, k) * gamma;
				u_stencil_cache[i] *= dt;
				u_stencil_cache[i] += u_center_cache[i];

				/* precompute stencil at i,j,k for g */
				g_stencil_cache[i] = g_left_cache[i];
				g_stencil_cache[i] += g_right_cache[i];
				g_stencil_cache[i] += g_in;
				g_stencil_cache[i] += g_out;
				/* stencil for g_up, g_down */
				g_stencil_cache[i] += g_center_cache[i - 1];
				g_stencil_cache[i] += g_center_cache[i + 1];
				g_stencil_cache[i] *= dt;
				g_stencil_cache[i] += 1.0;
			}
			for (i = 1; i < M - 1; i++) {
#pragma AP pipeline
				up_cache[i] =
				    u_center_cache[i - 1] * g_center_cache[i -
									   1] *
				    dt;
			}
			for (i = 1; i < M - 1; i++) {
#pragma AP pipeline
				u_stencil_center =
				    (u_stencil_cache[i] +
				     up_cache[i]) / g_stencil_cache[i];

				u_center_cache[i] = u_stencil_center;
				U_CENTER = u_stencil_center;
			}
		}
	}
}

void rician_deconv_deblur(double u[M][N][P], double f[M][N][P],
			  double g[M][N][P], double conv[M][N][P],
			  double Ksigma, double sigma, double lambda)
{
	uint32_t i, j, k;
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
		gaussian_blur(conv, Ksigma);

		/* calculate with rational cubic approx */
		for (k = 0; k < P; k++) {
			for (j = 0; j < N; j++) {
				for (i = 0; i < M; i++) {
#pragma AP pipeline
					conv[i][j][k] -=
					    cubic_approx(U(i, j, k), F(i, j, k),
							 sigma2);
				}
			}
		}
		gaussian_blur(conv, Ksigma);
		semi_implicit_update(u, g, conv, DEBLUR_DT, gamma);
	}
}

void rician_deconv_denoise(double u[M][N][P], double f[M][N][P],
			   double g[M][N][P], double sigma, double lambda)
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
