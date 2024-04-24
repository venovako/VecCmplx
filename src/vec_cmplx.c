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

void VSprintf(const int f, const char *const h, const VS v)
{
#ifdef _OPENMP
#pragma omp critical
  {
#endif /* _OPENMP */
  alignas(PVN_VECLEN) float s[VSL];

  if ((h ? dprintf(f, "\nL: %s\n", h) : 0) < 0)
    perror("dprintf 0");

  _mm512_store_ps(s, v);

  char a[17] = { '\0' };
  for (unsigned i = 0u; i < VSL; ++i)
    if (20 != dprintf(f, "%X: %s\n", i, pvn_stoa(a, s[i])))
      perror("dprintf 1");
#ifdef _OPENMP
  }
#endif /* _OPENMP */
}

void VDprintf(const int f, const char *const h, const VD v)
{
#ifdef _OPENMP
#pragma omp critical
  {
#endif /* _OPENMP */
  alignas(PVN_VECLEN) double d[VDL];

  if ((h ? dprintf(f, "\nL: %s\n", h) : 0) < 0)
    perror("dprintf 0");

  _mm512_store_pd(d, v);

  char a[26] = { '\0' };
  for (unsigned i = 0u; i < VDL; ++i)
    if (29 != dprintf(f, "%u: %s\n", i, pvn_dtoa(a, d[i])))
      perror("dprintf 1");
#ifdef _OPENMP
  }
#endif /* _OPENMP */
}

void MSprintf(const int f, const char *const h, const MS m)
{
#ifdef _OPENMP
#pragma omp critical
  {
#endif /* _OPENMP */
  if ((h ? dprintf(f, "\n%s: ", h) : 0) < 0)
    perror("dprintf 0");

  const unsigned u = MS2U(m);
  for (unsigned i = 0u, o = (1u << (VSL - 1u)); i < VSL; ++i) {
    if (1 != dprintf(f, "%c", ((u & o) ? '1' : '0')))
      perror("dprintf 1");
    o >>= 1u;
  }

  if (8 != dprintf(f, " (%04X)\n", u))
    perror("dprintf 2");
#ifdef _OPENMP
  }
#endif /* _OPENMP */
}

void MDprintf(const int f, const char *const h, const MD m)
{
#ifdef _OPENMP
#pragma omp critical
  {
#endif /* _OPENMP */
  if ((h ? dprintf(f, "\n%s: ", h) : 0) < 0)
    perror("dprintf 0");

  const unsigned u = MD2U(m);
  for (unsigned i = 0u, o = (1u << (VDL - 1u)); i < VDL; ++i) {
    if (1 != dprintf(f, "%c", ((u & o) ? '1' : '0')))
      perror("dprintf 1");
    o >>= 1u;
  }

  if (6 != dprintf(f, " (%02X)\n", u))
    perror("dprintf 2");
#ifdef _OPENMP
  }
#endif /* _OPENMP */
}
#endif /* ?VEC_CMPLX_TEST */
