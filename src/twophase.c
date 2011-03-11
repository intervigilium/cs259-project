#include "twophase.h"

static inline void neumann_bc(double curvature_motion_part[M * N * P])
{
	uint32_t i, j, k;
	for (j = 0; j < N; j++) {
		for (k = 0; k < P; k++) {
			CMP(0, j, k) = CMP(1, j, k);
			CMP(M - 1, j, k) = CMP(M - 2, j, k);
		}
	}

	for (i = 0; i < M; i++) {
		for (k = 0; k < P; k++) {
			CMP(i, 0, k) = CMP(i, 1, k);
			CMP(i, N - 1, k) = CMP(i, N - 2, k);
		}
	}

	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++) {
			CMP(i, j, 0) = CMP(i, j, 1);
			CMP(i, j, P - 1) = CMP(i, j, P - 2);
		}
	}

	CMP(0, 0, 0) = CMP(1, 1, 1);
	CMP(M - 1, 0, 0) = CMP(M - 2, 1, 1);
	CMP(0, N - 1, 0) = CMP(1, N - 2, 1);
	CMP(0, 0, P - 1) = CMP(1, 1, P - 2);
	CMP(M - 1, N - 1, 0) = CMP(M - 2, N - 2, 1);
	CMP(M - 1, 0, P - 1) = CMP(M - 2, 1, P - 2);
	CMP(0, N - 1, P - 1) = CMP(1, N - 2, P - 2);
	CMP(M - 1, N - 1, P - 1) = CMP(M - 2, N - 2, P - 2);
}

void two_phase_3d_op_explicit(double phi[M * N * P],
			      const double u0[M * N * P],
			      double curvature_motion_part[M * N * P],
			      double dt, double c1, double c2)
{
	double mu = TP_MU;
	double nu = TP_NU;
	double lambda1 = TP_LAMBDA1;
	double lambda2 = TP_LAMBDA2;

	double dx = 1.0;
	double dy = 1.0;
	double dz = 1.0;

	double dx2 = dx * 2.0;
	double dy2 = dy * 2.0;
	double dz2 = dz * 2.0;

	double Dx_p, Dx_m;
	double Dy_p, Dy_m;
	double Dz_p, Dz_m;
	double Dx_0, Dy_0, Dz_0;

	double Dxx, Dyy, Dzz;
	double Dxy, Dxz, Dyz;

	double Grad, K;

	double stencil[3][3][3];
	double numer, denom;

	uint32_t i, j, k, l;

	for (i = 1; i < M - 1; i++) {
		for (j = 1; j < N - 1; j++) {
			for (k = 1; k < P - 1; k++) {

				/* stencil code */
				stencil[0][0][0] = stencil[0][0][1];
				stencil[0][1][0] = stencil[0][1][1];
				stencil[0][2][0] = stencil[0][2][1];

				stencil[0][0][1] = stencil[0][0][2];
				stencil[0][1][1] = stencil[0][1][2];
				stencil[0][2][1] = stencil[0][2][2];

				stencil[0][0][2] = PHI(i - 1, j - 1, k + 1);
				stencil[0][1][2] = PHI(i - 1, j, k + 1);
				stencil[0][2][2] = PHI(i - 1, j + 1, k + 1);

				stencil[1][0][0] = stencil[1][0][1];
				stencil[1][1][0] = stencil[1][2][1];
				stencil[1][2][0] = stencil[1][2][1];

				stencil[1][0][1] = stencil[1][0][2];
				stencil[1][1][1] = stencil[1][1][2];
				stencil[1][2][1] = stencil[1][2][2];

				stencil[1][0][2] = PHI(i, j - 1, k + 1);
				stencil[1][1][2] = PHI(i, j, k + 1);
				stencil[1][2][2] = PHI(i, j + 1, k + 1);

				stencil[2][0][0] = stencil[2][0][1];
				stencil[2][1][0] = stencil[2][1][1];
				stencil[2][2][0] = stencil[2][2][1];

				stencil[2][0][1] = stencil[2][0][2];
				stencil[2][1][1] = stencil[2][1][2];
				stencil[2][2][1] = stencil[2][2][2];

				stencil[2][0][2] = PHI(i + 1, j - 1, k + 1);
				stencil[2][1][2] = PHI(i + 1, j, k + 1);
				stencil[2][2][2] = PHI(i + 1, j + 1, k + 1);

				/* regular calculation here */
				Dx_p =
				    (stencil[2][1][1] - stencil[1][1][1]) / dx;
				Dx_m =
				    (stencil[1][1][1] - stencil[0][1][1]) / dx;
				Dy_p =
				    (stencil[1][2][1] - stencil[1][1][1]) / dy;
				Dy_m =
				    (stencil[1][1][1] - stencil[1][0][1]) / dy;
				Dz_p =
				    (stencil[1][1][2] - stencil[1][1][1]) / dz;
				Dz_m =
				    (stencil[1][1][1] - stencil[1][1][0]) / dz;

				Dx_0 =
				    (stencil[2][1][1] - stencil[0][1][1]) / dx2;
				Dy_0 =
				    (stencil[1][2][1] - stencil[1][0][1]) / dy2;
				Dz_0 =
				    (stencil[1][1][2] - stencil[1][1][0]) / dz2;

				Dxx = (Dx_p - Dx_m) / dx;
				Dyy = (Dy_p - Dy_m) / dy;
				Dzz = (Dz_p - Dz_m) / dz;

				Dxy =
				    (stencil[2][2][1] - stencil[2][0][1] -
				     stencil[0][2][1] -
				     stencil[0][0][1]) / (4 * dx * dy);
				Dxz =
				    (stencil[2][1][2] - stencil[2][1][0] -
				     stencil[0][1][2] +
				     stencil[0][1][0]) / (4 * dx * dz);
				Dyz =
				    (stencil[1][2][2] - stencil[1][2][0] -
				     stencil[1][0][2] +
				     stencil[1][0][0]) / (4 * dy * dz);

				Grad = (SQR(Dx_0) + SQR(Dy_0) + SQR(Dz_0));
				denom = Grad;

				/* denom = denom^1.5 */
				for (l = 0; l < 3; l++) {
					denom *= denom;
				}
				q3_sqrt(denom);

				numer = (Dx_0 * Dx_0 * Dyy -
					 2.0 * Dx_0 * Dy_0 * Dxy +
					 Dy_0 * Dy_0 * Dxx + Dx_0 * Dx_0 * Dzz -
					 2.0 * Dx_0 * Dz_0 * Dxz +
					 Dz_0 * Dz_0 * Dxx + Dy_0 * Dy_0 * Dzz -
					 2.0 * Dy_0 * Dz_0 * Dyz +
					 Dz_0 * Dz_0 * Dyy);

				K = numer / denom;

				CMP(i, j, k) =
				    Grad * (mu * K +
					    lambda1 * (U(i, j, k) -
						       c1) * (U(i, j,
								k) - c1) -
					    lambda2 * (U(i, j, k) -
						       c2) * (U(i, j, k) - c2));
			}
		}
	}

	/* NeumannBC functionality inlined, pipeline this */
	neumann_bc(curvature_motion_part);

	/* pipeline this */
	for (i = 0; i < M * N * P; i++) {
		phi[i] += curvature_motion_part[i] * dt;
	}
}
