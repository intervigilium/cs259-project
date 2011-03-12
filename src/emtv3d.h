#ifndef _EMTV3D_H
#define _EMTV3D_H

#include "common.h"

#define LAMBDA_X(i, x_s, x_d, ray_len) (ray_len*((double)i-x_s)/(x_d-x_s))
#define LAMBDA_Y(j, y_s, y_d, ray_len) (ray_len*((double)j-y_s)/(y_d-y_s))
#define LAMBDA_Z(k, z_s, z_d, ray_len) (ray_len*((double)k-z_s)/(z_d-z_s))

#define ARR_IDX(arr, x, y, z) (arr[(x)+(y)*M+(z)*M*N])

void raytracer_projection_forward(double img[M * N * P], double D, uint32_t q,
				  double ss, double d, double dtheta,
				  double ssz, uint32_t qz,
				  double sino[M * N * P])
#endif
