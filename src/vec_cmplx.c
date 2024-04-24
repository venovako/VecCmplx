#include "vec_cmplx.h"

#ifdef VEC_CMPLX_TEST
int main(/* int argc, char *argv[] */)
{
  (void)printf("libvec_cmplx built on %s with %s for %s on %s ", __DATE__, VEC_CMPLX_COMPILER, VEC_CMPLX_OS, VEC_CMPLX_ARCH);
#ifdef NDEBUG
  (void)printf("with optimization level %d ", NDEBUG);
#else /* !NDEBUG */
  (void)printf("for debugging ");
#endif /* ?NDEBUG */
  (void)printf("and with OpenMP %d\n", _OPENMP);
  return (IS_STD_MXCSR ? EXIT_SUCCESS : EXIT_FAILURE);
}
#else /* !VEC_CMPLX_TEST */
static_assert(CHAR_BIT == 8, "CHAR_BIT != 8");
static_assert(sizeof(float) == 4, "sizeof(float) != 4");
static_assert(sizeof(double) == 8, "sizeof(double) != 8");
static_assert(sizeof(long double) >= 8, "sizeof(long double) < 8");
static_assert(sizeof(__float128) == 16, "sizeof(__float128) != 16");
#endif /* ?VEC_CMPLX_TEST */
