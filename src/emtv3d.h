#ifndef _EMTV3D_H
#define _EMTV3D_H

#include "common.h"

#define RT_PROJ_D 125
#define RT_PROJ_D1 64
#define RT_PROJ_Q 150
#define RT_PROJ_QZ 128
#define RT_PROJ_SS 0.5
#define RT_PROJ_SSZ 1.0
#define RT_PROJ_DTHETA 10.0

#define LAMBDA_X(i, x_s, x_d, ray_len) (ray_len*((double)i-x_s)/(x_d-x_s))
#define LAMBDA_Y(j, y_s, y_d, ray_len) (ray_len*((double)j-y_s)/(y_d-y_s))
#define LAMBDA_Z(k, z_s, z_d, ray_len) (ray_len*((double)k-z_s)/(z_d-z_s))

#define ARR_IDX(arr, x, y, z) (arr[(x)+(y)*M+(z)*M*N])

void raytracer_projection_forward(double img[M * N * P],
				  double sino[M * N * P]);

#endif
