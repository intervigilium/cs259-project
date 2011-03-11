#include "common.h"
#include "gaussian_blur.h"
#include "rician_deconv.h"
#include "twophase.h"

int main(double u1[M * N * P], double u2[M * N * P], double u3[M * N * P],
	  double u4[M * N * P], double in1, double in2, double in3)
{
	/* for deblur */
	/*
	   rician_deconv(u1, u2, u3, u4, in1, in2, in3, 1);
	 */

	/* for denoise */
	/*
	   rician_deconv(u1, u2, u3, u4, in1, in2, in3, 0);
	 */

	/* for registration */
	/*
	   gaussian_blur(u1, in1);
	   gaussian_blur(u2, in1);
	   gaussian_blur(u3, in1);
	 */

    /* for segmentation */
    /*
    two_phase_3d_op_explicit(u1, u2, u3, in1, in2, in3);
    */
	
	return 0;
}
