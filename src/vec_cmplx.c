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
  alignas(PVN_VECLEN) const float x[VSL] = { 15.0f, -14.0f, 13.0f, -12.0f, 11.0f, -10.0f, 9.0f, -8.0f, 7.0f, -6.0f, 5.0f, -4.0f, 3.0f, -2.0f, 1.0f, -0.0f };
  alignas(PVN_VECLEN) const float y[VSL] = { 0.0f, -1.0f, 2.0f, -3.0f, 4.0f, -5.0f, 6.0f, -7.0f, 8.0f, -9.0f, 10.0f, -11.0f, 12.0f, -13.0f, 14.0f, -15.0f };
  alignas(PVN_VECLEN) float z[VSL] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
  const ssize_t n = (ssize_t)(VSL / 2u);
  int info = 0;
  vec_cmul0f_(&n, x, y, z, &info);
  (void)printf("vec_cmul0f_=%d\n", info);
  return (IS_STD_MXCSR ? EXIT_SUCCESS : EXIT_FAILURE);
}
#else /* !VEC_CMPLX_TEST */
/*#include <sleefquad.h>*/

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

void vec_cmul0f_(const ssize_t *const n, const float *const x, const float *const y, float *const z, int *const info)
{
  PVN_ASSERT(n);
  PVN_ASSERT(x);
  PVN_ASSERT(y);
  PVN_ASSERT(z);
  PVN_ASSERT(info);
  if ((*info < 0) || (*info > 3))
    *info = -5;
  if (!PVN_IS_VECALIGNED(z))
    *info = -4;
  if (!PVN_IS_VECALIGNED(y))
    *info = -3;
  if (!PVN_IS_VECALIGNED(x))
    *info = -2;
  if (*n < 0)
    *info = -1;
  if (!*n || (*info < 0))
    return;
  const size_t m = ((size_t)*n << 1u);
  register const VI si = _mm512_set_epi32(15, 13, 11, 9, 7, 5, 3, 1, 14, 12, 10, 8, 6, 4, 2, 0);
  register const VI mi = _mm512_set_epi32(15, 7, 14, 6, 13, 5, 12, 4, 11, 3, 10, 2, 9, 1, 8, 0);
  register VS xp = ((*info & 1) ? _mm512_set_ps(x[1], x[1], x[1], x[1], x[1], x[1], x[1], x[1], x[0], x[0], x[0], x[0], x[0], x[0], x[0], x[0]) : _mm512_setzero_ps());
  register VS yp = ((*info & 2) ? _mm512_set_ps(y[1], y[1], y[1], y[1], y[1], y[1], y[1], y[1], y[0], y[0], y[0], y[0], y[0], y[0], y[0], y[0]) : _mm512_setzero_ps());
  for (size_t i = 0u, rem = m; i < m; (i += VSL), (rem -= VSL)) {
    if (rem < VSL) {
      alignas(PVN_VECLEN) float tail[VSL] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
      if (!(*info & 1)) {
        for (size_t j = 0u; j < rem; ++j)
          tail[j] = x[i + j];
        register const VS px = _mm512_load_ps(tail); VSP(px);
        xp = _mm512_permutexvar_ps(si, px);
      }
      VSP(xp);
      if (!(*info & 2)) {
        for (size_t j = 0u; j < rem; ++j)
          tail[j] = y[i + j];
        register const VS py = _mm512_load_ps(tail); VSP(py);
        yp = _mm512_permutexvar_ps(si, py);
      }
      VSP(yp);
    }
    else {
      if (!(*info & 1)) {
        register const VS px = _mm512_load_ps(x + i); VSP(px);
        xp = _mm512_permutexvar_ps(si, px);
      }
      VSP(xp);
      if (!(*info & 2)) {
        register const VS py = _mm512_load_ps(y + i); VSP(py);
        yp = _mm512_permutexvar_ps(si, py);
      }
      VSP(yp);
    }
#ifdef __AVX512DQ__
    register const VD xr = _mm512_cvtps_pd(_mm512_extractf32x8_ps(xp, 0)); VDP(xr);
    register const VD xi = _mm512_cvtps_pd(_mm512_extractf32x8_ps(xp, 1)); VDP(xi);
    register const VD yr = _mm512_cvtps_pd(_mm512_extractf32x8_ps(yp, 0)); VDP(yr);
    register const VD yi = _mm512_cvtps_pd(_mm512_extractf32x8_ps(yp, 1)); VDP(yi);
#else /* !__AVX512DQ__ */
    register const VD xr = _mm512_cvtps_pd(_mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(xp), 0))); VDP(xr);
    register const VD xi = _mm512_cvtps_pd(_mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(xp), 1))); VDP(xi);
    register const VD yr = _mm512_cvtps_pd(_mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(yp), 0))); VDP(yr);
    register const VD yi = _mm512_cvtps_pd(_mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(yp), 1))); VDP(yi);
#endif /* ?__AVX512DQ__ */
    register const VD zr = _mm512_fmsub_pd(xr, yr, _mm512_mul_pd(xi, yi)); VDP(zr);
    register const VD zi = _mm512_fmadd_pd(xr, yi, _mm512_mul_pd(xi, yr)); VDP(zi);
#ifdef __AVX512DQ__
    register const VS pz = _mm512_permutexvar_ps(mi, _mm512_insertf32x8(_mm512_zextps256_ps512(_mm512_cvtpd_ps(zr)), _mm512_cvtpd_ps(zi), 1)); VSP(pz);
#else /* !__AVX512DQ__ */
    register const VS pz = _mm512_permutexvar_ps(mi, _mm512_castpd_ps(_mm512_insertf64x4(_mm512_castps_pd(_mm512_zextps256_ps512(_mm512_cvtpd_ps(zr))), _mm256_castps_pd(_mm512_cvtpd_ps(zi)), 1))); VSP(pz);
#endif /* ?__AVX512DQ__ */
    if (rem < VSL) {
      alignas(PVN_VECLEN) float tail[VSL] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
      _mm512_store_ps(tail, pz);
      for (size_t j = 0u; j < rem; ++j)
        z[i + j] = tail[j];
    }
    else
      _mm512_store_ps((z + i), pz);
  }
}

void vec_cmul0_(const ssize_t *const n, const double *const x, const double *const y, double *const z, int *const info)
{
  PVN_ASSERT(n);
  PVN_ASSERT(x);
  PVN_ASSERT(y);
  PVN_ASSERT(z);
  PVN_ASSERT(info);
  if ((*info < 0) || (*info > 3))
    *info = -5;
  if (!PVN_IS_VECALIGNED(z))
    *info = -4;
  if (!PVN_IS_VECALIGNED(y))
    *info = -3;
  if (!PVN_IS_VECALIGNED(x))
    *info = -2;
  if (*n < 0)
    *info = -1;
  if (!*n || (*info < 0))
    return;
  const size_t m = ((size_t)*n << 1u);
  /* TODO */
}
#endif /* ?VEC_CMPLX_TEST */
