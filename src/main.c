#include "gaussian_blur.h"
#include "rician_deconv.h"

void main(double u1[M][N][P], double u2[M][N][P], double u3[M][N][P],
	  double u4[M][N][P], double Ksigma, double sigma, double lambda)
{
	/* for deblur */
	/*
	   rician_deconv(u1, u2, u3, u4, Ksigma, sigma, lambda, 1);
	 */

	/* for denoise */
	/*
	   rician_deconv(u1, u2, u3, u4, Ksigma, sigma, lambda, 0);
	 */

	/* for registration */
	/*
	   gaussian_blur(u1, Ksigma);
	   gaussian_blur(u2, Ksigma);
	   gaussian_blur(u3, Ksigma);
	 */
}
