#include <autopilot_tech.h>

#define uint2_t uint2
#define uint32_t uint32
#define uint64_t uint64

#define M 60
#define N 60
#define P 60

#define GAUSSIAN_NUMSTEPS 3

#define MAX_ITERATIONS 10

#define SQR(x) ((x)*(x))
#define U(a,b,c) (u[a+b*N+c*M*N])
#define G(a,b,c) (g[a+b*N+c*M*N])

#define U_CENTER U(i,j,k)
#define U_LEFT U(i,j-1,k)
#define U_RIGHT U(i,j+1,k)
#define U_UP U(i-1,j,k)
#define U_DOWN U(i+1,j,k)
#define U_IN U(i,j,k-1)
#define U_OUT U(i,j,k+1)

#define G_CENTER G(i,j,k)
#define G_LEFT G(i,j-1,k)
#define G_RIGHT G(i,j+1,k)
#define G_UP G(i-1,j,k)
#define G_DOWN G(i+1,j,k)
#define G_IN G(i,j,k-1)
#define G_OUT G(i,j,k+1)

double q3_sqrt(double num);

double fast_fabs(double num);

void array_copy(const double src[M*N*P], double dst[M*N*P]);
