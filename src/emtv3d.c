#include "emtv3d.h"

static inline double raytracer_forward(double src_x, double src_y, double src_z,
				       double det_x, double det_y, double det_z,
				       double img[M * N * P])
{
	double sino = 0.0;

	double ray_dir_x, ray_dir_y, ray_dir_z;;
	double ray_len;

	double len_x, len_y, len_z;
	double absvalue_x, absvalue_y, absvalue_z;
	double lambda_min, lambda_max;
	int signx, signy, signz;

	int v_x, v_y, v_z;
	int u_x, u_y, u_z;

	int index;
	double x_min, y_min, z_min;

	double lambda_0;
	double lambda_x, lambda_y, lambda_z;

	double Tx, Ty, Tz;
	double length;

	/* get direction */
	ray_dir_x = det_x - src_x;
	ray_dir_y = det_y - src_y;
	ray_dir_z = det_z - src_z;
	/* get ray length */
	ray_len = SQR(ray_dir_x) + SQR(ray_dir_y) + SQR(ray_dir_z);
	ray_len = q3_sqrt(ray_len);

	signx = (ray_dir_x > 0) ? 1 : -1;
	signy = (ray_dir_y > 0) ? 1 : -1;
	signz = (ray_dir_z > 0) ? 1 : -1;

	absvalue_x = fast_fabs(ray_dir_x);
	absvalue_y = fast_fabs(ray_dir_y);
	absvalue_z = fast_fabs(ray_dir_z);

	/* get x=1 Lx Ly Lz */
	len_x = (double)((absvalue_x > 1.e-4) ? (ray_len / absvalue_x) : 1.e6);
	len_y = (double)((absvalue_y > 1.e-4) ? (ray_len / absvalue_y) : 1.e6);
	len_z = (double)((absvalue_z > 1.e-4) ? (ray_len / absvalue_z) : 1.e6);

	/* get the entry and exit point between Ray & Image
	 * distance between source and entry point */
	lambda_min = MAX(MAX(0.0,
			     MIN(LAMBDA_X(0, src_x, det_x, ray_len),
				 LAMBDA_X(M, src_x, det_x, ray_len))
			 ),
			 MAX(MIN
			     (LAMBDA_Y(0, src_y, det_y, ray_len),
			      LAMBDA_Y(N, src_y, det_y, ray_len)),
			     MIN(LAMBDA_Z(0, src_z, det_z, ray_len),
				 LAMBDA_Z(P, src_z, det_z, ray_len))
			 )
	    );

	/* distance between source and exit point */
	lambda_max = MIN(MIN(ray_len,
			     MAX(LAMBDA_X(0, src_x, det_x, ray_len),
				 LAMBDA_X(M, src_x, det_x, ray_len))
			 ),
			 MIN(MAX
			     (LAMBDA_Y(0, src_y, det_y, ray_len),
			      LAMBDA_Y(N, src_y, det_y, ray_len)),
			     MAX(LAMBDA_Z(0, src_z, det_z, ray_len),
				 LAMBDA_Z(P, src_z, det_z, ray_len))
			 )
	    );

	/* no intersection between Ray & Image */
	if (lambda_min >= lambda_max)
		return sino;

	/* get the position of entry point */
	x_min = src_x + lambda_min * ray_dir_x / ray_len;
	y_min = src_y + lambda_min * ray_dir_y / ray_len;
	z_min = src_z + lambda_min * ray_dir_z / ray_len;

	x_min = (x_min) > 0 ? x_min : 0;
	y_min = (y_min) > 0 ? y_min : 0;
	z_min = (z_min) > 0 ? z_min : 0;

	/* v current pixel, u next pixel */
	if (fast_fabs(z_min - (int)(z_min + 0.5)) < 1.e-4) {
		/* integer pixel */
		v_z = (int)(z_min + 0.5);
		u_z = v_z + signz;
		index = 3;
		if (ray_dir_z < 0.0) {
			v_z -= 1;
		}
	} else {
		/* frac pixel */
		v_z = (int)(z_min);
		if (ray_dir_z < 0.0)
			u_z = (int)(z_min);
		else
			u_z = (int)(z_min + 1);
	}

	if (fast_fabs(y_min - (int)(y_min + 0.5)) < 1.e-4) {
		v_y = (int)(y_min + 0.5);
		u_y = v_y + signy;
		index = 2;
		if (ray_dir_y < 0.0)
			v_y -= 1;
	} else {
		v_y = (int)(y_min);
		if (ray_dir_y < 0.0)
			u_y = (int)(y_min);
		else
			u_y = (int)(y_min + 1);
	}

	if (fast_fabs(x_min - (int)(x_min + 0.5)) < 1.e-4) {
		v_x = (int)(x_min + 0.5);
		u_x = v_x + signx;
		index = 1;
		if (ray_dir_x < 0.0)
			v_x -= 1;
	} else {
		v_x = (int)(x_min);
		if (ray_dir_x < 0.0)
			u_x = (int)(x_min);
		else
			u_x = (int)(x_min + 1);
	}

	/* set the beginning pixel
	 * lambda_x .y .z are the distance between source point and next ray point */
	lambda_0 = lambda_min;
	lambda_x =
	    (absvalue_x < 1e-4) ? 1e6 : LAMBDA_X(u_x, src_x, det_x, ray_len);
	lambda_y =
	    (absvalue_y < 1e-4) ? 1e6 : LAMBDA_Y(u_y, src_y, det_y, ray_len);
	lambda_z =
	    (absvalue_z < 1e-4) ? 1e6 : LAMBDA_Z(u_z, src_z, det_z, ray_len);

	/* the main loop
	 * large distance on X axis */
	if (fast_fabs(ray_dir_x) > fast_fabs(ray_dir_y)) {
		/* intersection with x */
		if (index == 1) {
			/* align the x axis without intersection of y and z */
			if (MIN(lambda_z, lambda_y) > lambda_max) {
				while (1) {
					sino += ARR_IDX(img, v_x, v_y, v_z);

					v_x += signx;

					if (v_x >= M || v_x < 0)
						return sino * len_x;
				}
			}
			/* cross z firstly */
			if (lambda_z < lambda_y) {
				Ty = lambda_y - lambda_z;
				Tz = lambda_z - lambda_0;
				Tx = lambda_x - lambda_0;

				/* intersection Z firstly */
				if (Tz < Tx) {
					length = Tz;
					sino +=
					    ARR_IDX(img, v_x, v_y,
						    v_z) * length;

					Tx -= Tz;
					Tz = len_z;

					/* update v_z */
					v_z += signz;

					if (v_z >= P || v_z < 0)
						return sino;
				} else {
					length = Tx;
					sino +=
					    ARR_IDX(img, v_x, v_y,
						    v_z) * length;

					Tz -= Tx;
					Tx = len_x;
					/* update v_x */
					v_x += signx;
					/* exit img */
					if (v_x >= M || v_x < 0)
						return sino;

					/* intersection with z */
					while (Tz >= Tx) {
						length = Tx;
						sino +=
						    ARR_IDX(img, v_x, v_y,
							    v_z) * length;

						Tz -= Tx;
						v_x += signx;
						if (v_x >= M || v_x < 0)
							return sino;
					}

					length = Tz;
					sino +=
					    ARR_IDX(img, v_x, v_y,
						    v_z) * length;

					Tx -= Tz;
					Tz = len_z;
					v_z += signz;
					if (v_z >= P || v_z < 0)
						return sino;
				}
			} else {
				Tz = lambda_z - lambda_y;
				Ty = lambda_y - lambda_0;
				Tx = lambda_x - lambda_0;

				if (Ty < Tx) {
					length = Ty;
					sino +=
					    ARR_IDX(img, v_x, v_y,
						    v_z) * length;

					Tx -= Ty;
					Ty = len_y;
					v_y += signy;
					if (v_y >= N || v_y < 0)
						return sino;
				} else {
					length = Tx;
					sino +=
					    ARR_IDX(img, v_x, v_y,
						    v_z) * length;

					Ty -= Tx;
					Tx = len_x;
					v_x += signx;
					if (v_x >= M || v_x < 0)
						return sino;

					while (Ty >= Tx) {
						length = Tx;
						sino +=
						    ARR_IDX(img, v_x, v_y,
							    v_z) * length;

						Ty -= Tx;
						v_x += signx;
						if (v_x >= M || v_x < 0)
							return sino;
					}

					length = Ty;
					sino +=
					    ARR_IDX(img, v_x, v_y,
						    v_z) * length;

					Tx -= Ty;
					Ty = len_y;
					v_y += signy;
					if (v_y >= N || v_y < 0)
						return sino;
				}
			}
		} else if (index == 2) {
			/* intersection with y */
			Ty = len_y;
			Tz = lambda_z - lambda_0;
			Tx = lambda_x - lambda_0;
		} else {
			/* intersection with z */
			Tz = len_z;
			Ty = lambda_y - lambda_0;
			Tx = lambda_x - lambda_0;
		}

		/* intersection begin with Y, Z */
		while (1) {
			if (Tz < Ty) {
				/* intersection with Z firstly */
				Ty -= Tz;
				if (Tz < Tx) {
					length = Tz;
					sino +=
					    ARR_IDX(img, v_x, v_y,
						    v_z) * length;

					Tx -= Tz;
					Tz = len_z;
					v_z += signz;
					if (v_z >= P || v_z < 0)
						return sino;
				} else {
					length = Tx;
					sino +=
					    ARR_IDX(img, v_x, v_y,
						    v_z) * length;

					Tz -= Tx;
					Tx = len_x;
					v_x += signx;
					if (v_x >= M || v_x < 0)
						return sino;

					while (Tz >= Tx) {
						length = Tx;
						sino +=
						    ARR_IDX(img, v_x, v_y,
							    v_z) * length;

						Tz -= Tx;
						v_x += signx;
						if (v_x >= M || v_x < 0)
							return sino;
					}
					length = Tz;
					sino +=
					    ARR_IDX(img, v_x, v_y,
						    v_z) * length;

					Tx -= Tz;
					Tz = len_z;
					v_z += signz;
					if (v_z >= P || v_z < 0)
						return sino;
				}
			} else {
				/* intersection with Y firstly */
				Tz -= Ty;
				if (Ty < Tx) {
					length = Ty;
					sino +=
					    ARR_IDX(img, v_x, v_y,
						    v_z) * length;

					Tx -= Ty;
					Ty = len_y;
					v_y += signy;
					if (v_y >= N || v_y < 0)
						return sino;
				} else {
					length = Tx;
					sino +=
					    ARR_IDX(img, v_x, v_y,
						    v_z) * length;

					Ty -= Tx;
					Tx = len_x;
					v_x += signx;
					if (v_x >= M || v_x < 0)
						return sino;

					while (Ty >= Tx) {
						length = Tx;
						sino +=
						    ARR_IDX(img, v_x, v_y,
							    v_z) * length;

						Ty -= Tx;
						v_x += signx;
						if (v_x >= M || v_x < 0)
							return sino;
					}

					length = Ty;
					sino +=
					    ARR_IDX(img, v_x, v_y,
						    v_z) * length;

					Tx -= Ty;
					Ty = len_y;
					v_y += signy;
					if (v_y >= N || v_y < 0)
						return sino;
				}
			}
		}
	} else {
		/* large distance on Y axis */
		/* intersection with y integer pixel */
		if (index == 2) {
			if (MIN(lambda_z, lambda_x) > lambda_max) {
				while (1) {
					sino += ARR_IDX(img, v_x, v_y, v_z);

					v_y += signy;
					if (v_y >= N || v_y < 0)
						return sino * len_y;
				}
			}
			if (lambda_z < lambda_x) {
				Tx = lambda_x - lambda_z;
				Tz = lambda_z - lambda_0;
				Ty = lambda_y - lambda_0;
				if (Tz < Ty) {
					length = Tz;
					sino +=
					    ARR_IDX(img, v_x, v_y,
						    v_z) * length;

					Ty -= Tz;
					Tz = len_z;
					v_z += signz;
					if (v_z >= P || v_z < 0)
						return sino;
				} else {
					length = Ty;
					sino +=
					    ARR_IDX(img, v_x, v_y,
						    v_z) * length;

					Tz -= Ty;
					Ty = len_y;
					v_y += signy;
					if (v_y >= N || v_y < 0)
						return sino;
					while (Tz >= Ty) {
						length = Ty;
						sino +=
						    ARR_IDX(img, v_x, v_y,
							    v_z) * length;

						Tz -= Ty;
						v_y += signy;
						if (v_y >= N || v_y < 0)
							return sino;
					}
					length = Tz;
					sino +=
					    ARR_IDX(img, v_x, v_y,
						    v_z) * length;

					Ty -= Tz;
					Tz = len_z;
					v_z += signz;
					if (v_z >= P || v_z < 0)
						return sino;
				}
			} else {
				Tz = lambda_z - lambda_x;
				Tx = lambda_x - lambda_0;
				Ty = lambda_y - lambda_0;
				if (Tx < Ty) {
					length = Tx;
					sino +=
					    ARR_IDX(img, v_x, v_y,
						    v_z) * length;

					Ty -= Tx;
					Tx = len_x;
					v_x += signx;
					if (v_x >= M || v_x < 0)
						return sino;
				} else {
					length = Ty;
					sino +=
					    ARR_IDX(img, v_x, v_y,
						    v_z) * length;

					Tx -= Ty;
					Ty = len_y;
					v_y += signy;
					if (v_y >= N || v_y < 0)
						return sino;
					while (Tx >= Ty) {
						length = Ty;
						sino +=
						    ARR_IDX(img, v_x, v_y,
							    v_z) * length;

						Tx -= Ty;
						v_y += signy;
						if (v_y >= N || v_y < 0)
							return sino;
					}
					length = Tx;
					sino +=
					    ARR_IDX(img, v_x, v_y,
						    v_z) * length;

					Ty -= Tx;
					Tx = len_x;
					v_x += signx;
					if (v_x >= M || v_x < 0)
						return sino;
				}
			}
		} else if (index == 1) {
			Tx = len_x;
			Tz = lambda_z - lambda_0;
			Ty = lambda_y - lambda_0;
		} else {
			Tz = len_z;
			Tx = lambda_x - lambda_0;
			Ty = lambda_y - lambda_0;
		}

		while (1) {
			if (Tz < Tx) {
				Tx -= Tz;
				if (Tz < Ty) {
					length = Tz;
					sino +=
					    ARR_IDX(img, v_x, v_y,
						    v_z) * length;

					Ty -= Tz;
					Tz = len_z;
					v_z += signz;
					if (v_z >= P || v_z < 0)
						return sino;
				} else {
					length = Ty;
					sino +=
					    ARR_IDX(img, v_x, v_y,
						    v_z) * length;

					Tz -= Ty;
					Ty = len_y;
					v_y += signy;
					if (v_y >= N || v_y < 0)
						return sino;
					while (Tz >= Ty) {
						length = Ty;
						sino +=
						    ARR_IDX(img, v_x, v_y,
							    v_z) * length;

						Tz -= Ty;
						v_y += signy;
						if (v_y >= N || v_y < 0)
							return sino;
					}
					length = Tz;
					sino +=
					    ARR_IDX(img, v_x, v_y,
						    v_z) * length;

					Ty -= Tz;
					Tz = len_z;
					v_z += signz;
					if (v_z >= P || v_z < 0)
						return sino;
				}
			} else {
				Tz -= Tx;
				if (Tx < Ty) {
					length = Tx;
					sino +=
					    ARR_IDX(img, v_x, v_y,
						    v_z) * length;

					Ty -= Tx;
					Tx = len_x;
					v_x += signx;
					if (v_x >= M || v_x < 0)
						return sino;
				} else {
					length = Ty;
					sino +=
					    ARR_IDX(img, v_x, v_y,
						    v_z) * length;

					Tx -= Ty;
					Ty = len_y;
					v_y += signy;
					if (v_y >= N || v_y < 0)
						return sino;

					while (Tx >= Ty) {
						length = Ty;
						sino +=
						    ARR_IDX(img, v_x, v_y,
							    v_z) * length;

						Tx -= Ty;
						v_y += signy;
						if (v_y >= N || v_y < 0)
							return sino;
					}

					length = Tx;
					sino +=
					    ARR_IDX(img, v_x, v_y,
						    v_z) * length;

					Ty -= Tx;
					Tx = len_x;
					v_x += signx;
					if (v_x >= M || v_x < 0)
						return sino;
				}
			}
		}
	}
}

void raytracer_projection_forward(double img[M * N * P], double sino[M * N * P])
{
	/* formerly parameters, actually constants from CPU_Routine.c */
	double D = RT_PROJ_D;
	double d = RT_PROJ_D1;
	double ss = RT_PROJ_SS;
	double ssz = RT_PROJ_SSZ;

	uint32_t q = RT_PROJ_Q;
	uint32_t qz = RT_PROJ_QZ;

	double src_x, src_y, src_z;
	double det_x, det_y, det_z;

	uint32_t NS = (uint32_t) (359.0 / RT_PROJ_DTHETA);
	double thetaS = 0.;
	double dtheta2 = (double)(RT_PROJ_DTHETA * PI / 180.);
	double ss2 = (double)(ss * PI / 180.);
	double aa;
	uint32_t i, j, k;

	for (i = 0; i <= NS; i++) {
		src_x = D * fast_cos(thetaS) + d;
		src_y = D * fast_sin(thetaS) + d;
		src_z = d;
		for (j = -q; j < q + 1; j++) {
			det_x = D * fast_cos(j * ss2 + PI + thetaS) + d;
			det_y = D * fast_sin(j * ss2 + PI + thetaS) + d;
			for (k = -qz; k < qz + 1; k++) {
				det_z = k * ssz + d;
				aa = raytracer_forward(src_x, src_y, src_z,
						       det_x, det_y, det_z,
						       img);
				ARR_IDX(sino, (j + q), (k + qz), i) = aa;
			}
		}
		thetaS += dtheta2;
	}
}
