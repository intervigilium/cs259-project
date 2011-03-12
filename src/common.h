#ifndef _COMMON_H
#define _COMMON_H

#include <autopilot_tech.h>

#define uint4_t uint4
#define uint32_t uint32
#define uint64_t uint64

#define M 181
#define N 217
#define P 181

#define PI 3.14159265

#define SQR(x) ((x)*(x))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

double q3_sqrt(double num);

double fast_fabs(double num);

double fast_sin(double num);

double fast_cos(double num);

void array_copy(const double src[M * N * P], double dst[M * N * P]);

#endif
