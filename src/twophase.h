#ifndef _TWOPHASE_H
#define _TWOPHASE_H

#include "common.h"

#define TP_MU 11704.5
#define TP_NU 0.0
#define TP_LAMBDA1 1.0
#define TP_LAMBDA2 1.0
#define TP_EPSILON 10e-10

#define CMP(a,b,c) (curvature_motion_part[a+b*N+c*N*M])
#define PHI(a,b,c) (phi[a+b*N+c*N*M])
#define U(a,b,c) (u0[a+b*N+c*N*M])

void two_phase_3d_op_explicit(double phi[M * N * P], const double u0[M * N * P],
			      double curvature_motion_part[M * N * P],
			      double dt, double c1, double c2);

#endif
