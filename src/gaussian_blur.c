#include "gaussian_blur.h"

void gaussian_blur(double u[M * N * P], double Ksigma)
{
	double lambda = (Ksigma * Ksigma) / (2.0 * GAUSSIAN_NUMSTEPS);
	double nu =
	    (1.0 + 2.0 * lambda - q3_sqrt(1.0 + 4.0 * lambda)) / (2.0 * lambda);
	double BoundaryScale = 1.0 / (1.0 - nu);
	double PostScale = 1;
	double m_cache[M], n_cache[N], p_cache[P];
	double r;
	uint32_t steps, i, j, k, idx;
	uint32_t plane, col, row;

	for (steps = 0; steps < 3 * GAUSSIAN_NUMSTEPS; steps++) {
#pragma AP unroll
		/* PostScale = (nu / lambda) ^ (3*GAUSSIAN_NUMSTEPS) */
		PostScale *= nu / lambda;
	}

	for (steps = 0; steps < GAUSSIAN_NUMSTEPS; steps++) {
		/* all of these loops have data dependencies
		 * due to u[i][j][k] writebacks */
		/* move up by one plane, ie k++ */

		/* downward */
		for (k = 0; k < P * N; k++) {
#pragma AP pipeline
			plane = k / N;
			col = k % N;
			U(0, col, plane) *= BoundaryScale;
		}
		for (k = 0; k < P; k++) {
			for (j = 0; j < N; j++) {
#pragma AP pipeline
				n_cache[j] = U(0, j, k);
			}
			for (i = 1; i < M; i++) {
				for (j = 0; j < N; j++) {
#pragma AP pipeline
					r = U(i, j, k) + nu * n_cache[j];
					U(i, j, k) = r;
					n_cache[j] = r;
				}
			}
		}

		/* upward */
		for (k = 0; k < P * N; k++) {
#pragma AP pipeline
			plane = k / N;
			col = k % N;
			U(M - 1, col, plane) *= BoundaryScale;
		}
		for (k = 0; k < P; k++) {
			for (j = 0; j < N; j++) {
#pragma AP pipeline
				n_cache[j] = U(M - 1, j, k);
			}
			for (i = 0; i < M - 1; i++) {
				/* retarded autopilot cant count backwards */
				idx = M - 2 - i;
				for (j = 0; j < N; j++) {
#pragma AP pipeline
					r = U(idx, j, k) + nu * n_cache[j];
					U(idx, j, k) = r;
					n_cache[j] = r;
				}
			}
		}

		/* right */
		for (k = 0; k < P * M; k++) {
#pragma AP pipeline
			plane = k / M;
			row = k % M;
			U(row, 0, plane) *= BoundaryScale;
		}
		for (i = 0; i < M; i++) {
			for (k = 0; k < P; k++) {
#pragma AP pipeline
				p_cache[k] = U(i, 0, k);
			}
			for (j = 1; j < N; j++) {
				for (k = 0; k < P; k++) {
#pragma AP pipeline
					r = U(i, j, k) + nu * p_cache[k];
					U(i, j, k) = r;
					p_cache[k] = r;
				}
			}
		}

		/* left */
		for (k = 0; k < P * M; k++) {
#pragma AP pipeline
			plane = k / M;
			row = k % M;
			U(row, N - 1, plane) *= BoundaryScale;
		}
		for (i = 0; i < M; i++) {
			for (k = 0; k < P; k++) {
#pragma AP pipeline
				p_cache[k] = U(i, N - 1, k);
			}
			for (j = 0; j < N - 1; j++) {
				idx = N - 2 - j;
				for (k = 0; k < P; k++) {
#pragma AP pipeline
					r = U(i, idx, k) + nu * p_cache[k];
					U(i, idx, k) = r;
					p_cache[k] = r;
				}
			}
		}

		/* out */
		for (j = 0; j < N * M; j++) {
#pragma AP pipeline
			col = j / M;
			row = j % M;
			U(row, col, 0) *= BoundaryScale;
		}
		for (j = 0; j < N; j++) {
			for (i = 0; i < M; i++) {
#pragma AP pipeline
				m_cache[i] = U(i, j, 0);
			}
			for (k = 1; k < P; k++) {
				for (i = 0; i < M; i++) {
#pragma AP pipeline
					r = U(i, j, k) + nu * m_cache[i];
					U(i, j, k) = r;
					m_cache[i] = r;
				}
			}
		}

		/* in */
		for (j = 0; j < N * M; j++) {
#pragma AP pipeline
			col = j / M;
			row = j % M;
			U(row, col, P - 1) *= BoundaryScale;
		}
		for (j = 0; j < N; j++) {
			for (i = 0; i < M; i++) {
#pragma AP pipeline
				m_cache[i] = U(i, j, P - 1);
			}
			for (k = 0; k < P - 1; k++) {
				idx = P - 2 - k;
				for (i = 0; i < M; i++) {
#pragma AP pipeline
					r = U(i, j, idx) + nu * m_cache[i];
					U(i, j, idx) = r;
					m_cache[i] = r;
				}
			}
		}
	}

	for (i = 0; i < M * N * P; i++) {
#pragma AP pipeline
		u[i] *= PostScale;
	}
}
