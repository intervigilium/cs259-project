#include <autopilot_tech.h>

#define uint2_t uint2
#define uint32_t uint32
#define uint64_t uint64

#define M 60
#define N 60
#define P 60

#define SQR(x) ((x)*(x))

double q3_sqrt(double num);

double fast_fabs(double num);

void array_copy(const double src[M*N*P], double dst[M*N*P]);
