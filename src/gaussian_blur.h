#ifndef _GAUSSIAN_BLUR_H
#define _GAUSSIAN_BLUR_H

#include "common.h"

#define GAUSSIAN_NUMSTEPS 3

#define U(a,b,c) (u[(a)][(b)][(c)])

void gaussian_blur(double u[M][N][P], double Ksigma);

#endif
