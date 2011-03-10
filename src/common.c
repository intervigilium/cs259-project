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

void array_copy(const double src[M*N*P], double dst[M*N*P])
{
	int i;
	for (i = 0; i < M*N*P; i++) {
#pragma AP pipeline
		dst[i] = src[i];
	}
}
