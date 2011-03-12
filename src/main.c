#include "common.h"
#include "gaussian_blur.h"
#include "rician_deconv.h"
#include "twophase.h"
#include "emtv3d.h"

int main(double u1[M * N * P], double u2[M * N * P], double u3[M * N * P],
	 double u4[M * N * P], double in1, double in2, double in3,
	 uint4_t select)
{
#pragma AP interface ap_bus port=u1 pipeline
#pragma AP interface ap_bus port=u2 pipeline
#pragma AP interface ap_bus port=u3 pipeline
#pragma AP interface ap_bus port=u4 pipeline
	if (select == 0) {
		/* for deblur */
		rician_deconv_deblur(u1, u2, u3, u4, in1, in2, in3);
	} else if (select == 1) {
		/* for denoise */
		rician_deconv_denoise(u1, u2, u3, in1, in2);
	} else if (select == 2) {
		/* for registration */
		gaussian_blur(u1, in1);
		gaussian_blur(u2, in1);
		gaussian_blur(u3, in1);
	} else if (select == 3) {
		/* for segmentation */
		two_phase_3d_op_explicit(u1, u2, u3, in1, in2, in3);
	} else if (select == 4) {
		/* for emtv3d */
//              raytracer_projection_forward(u1, u2);
	}

	return 0;
}
