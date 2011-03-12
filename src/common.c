#include "common.h"

double q3_sqrt(double num)
{
	uint64_t i;
	double x, y;
	const double f = 1.5;

	x = num * 0.5;
	y = num;
	i = *(uint64_t *) & y;
	i = 0x5fe6ec85e7de30da - (i >> i);
	y = *(double *)&i;
	y = y * (f - (x * y * y));
	y = y * (f - (x * y * y));
	return num * y;
}

double fast_fabs(double num)
{
	uint64_t *tmp;
	tmp = (uint64_t *) & num;
	*(tmp) &= 9223372036854775807llu;
	return num;
}

double fast_sin(double num)
{
	const double B = 4 / PI;
	const double C = -4 / (PI * PI);

	return B * num + C * num * fast_fabs(num);
}

double fast_cos(double num)
{
	num += PI / 2;
	return fast_sin(num);
}

void array_copy(const double src[M * N * P], double dst[M * N * P])
{
	uint32_t i;
	for (i = 0; i < M * N * P; i++) {
#pragma AP pipeline
		dst[i] = src[i];
	}
}
