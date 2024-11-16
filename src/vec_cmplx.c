#include "vec_cmplx.h"

#ifdef VEC_CMPLX_TEST
int main(/* int argc, char *argv[] */)
{
  (void)fprintf(stdout, "libvec_cmplx built on %s with %s for %s on %s ", __DATE__, VEC_CMPLX_COMPILER, VEC_CMPLX_OS, VEC_CMPLX_ARCH);
#ifdef NDEBUG
  (void)fprintf(stdout, "with optimization level %d ", NDEBUG);
#else /* !NDEBUG */
  (void)fprintf(stdout, "for debugging ");
#endif /* ?NDEBUG */
  (void)fprintf(stdout, "and with OpenMP %d\n", _OPENMP);
  (void)fflush(stdout);
  alignas(PVN_VECLEN) const float x[VSL] = { 15.0f, -14.0f, 13.0f, -12.0f, 11.0f, -10.0f, 9.0f, -8.0f, 7.0f, -6.0f, 5.0f, -4.0f, 3.0f, -2.0f, 1.0f, -0.0f };
  alignas(PVN_VECLEN) const float y[VSL] = { 0.0f, -1.0f, 2.0f, -3.0f, 4.0f, -5.0f, 6.0f, -7.0f, 8.0f, -9.0f, 10.0f, -11.0f, 12.0f, -13.0f, 14.0f, -15.0f };
  alignas(PVN_VECLEN) float z[VSL] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
  const ssize_t n = (ssize_t)(VSL >> 1u);
  int info = 0;
  vec_cmul0_(&n, x, y, z, &info);
  (void)fflush(stdout);
  (void)fprintf(stdout, "vec_cmul0_=%d\n", info);
#ifdef __AVX512VL__
  const ssize_t inc = 2;
  vec_cmul1_(&n, x, (x + 1), &inc, y, (y + 1), &inc, z, (z + 1), &inc, &info);
  (void)fflush(stdout);
  (void)fprintf(stdout, "vec_cmul1_=%d\n", info);
#endif /* __AVX512VL__ */
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

void vec_cmul0_(const ssize_t *const n, const float *const x, const float *const y, float *const z, int *const info)
{
  PVN_ASSERT(n);
  PVN_ASSERT(x);
  PVN_ASSERT(y);
  PVN_ASSERT(z);
  PVN_ASSERT(info);
  if (*n < 0)
    *info = -1;
  if (!*n || (*info < 0))
    return;
  if (PVN_IS_VECALIGNED(x))
    *info |= 4;
  if (PVN_IS_VECALIGNED(y))
    *info |= 8;
  if (PVN_IS_VECALIGNED(z))
    *info |= 16;
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
        register const VS px = ((*info & 4) ? _mm512_load_ps(x + i) : _mm512_loadu_ps(x + i)); VSP(px);
        xp = _mm512_permutexvar_ps(si, px);
      }
      VSP(xp);
      if (!(*info & 2)) {
        register const VS py = ((*info & 8) ? _mm512_load_ps(y + i) : _mm512_loadu_ps(y + i)); VSP(py);
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
    else if (*info & 16)
      _mm512_store_ps((z + i), pz);
    else
      _mm512_storeu_ps((z + i), pz);
  }
}

#ifdef __AVX512VL__
void vec_cmul1_(const ssize_t *const n, const float *const rx, const float *const ix, const ssize_t *const incx, const float *const ry, const float *const iy, const ssize_t *const incy, float *const rz, float *const iz, const ssize_t *const incz, int *const info)
{
  PVN_ASSERT(n);
  PVN_ASSERT(rx);
  PVN_ASSERT(ix);
  PVN_ASSERT(incx);
  PVN_ASSERT(ry);
  PVN_ASSERT(iy);
  PVN_ASSERT(incy);
  PVN_ASSERT(rz);
  PVN_ASSERT(iz);
  PVN_ASSERT(incz);
  PVN_ASSERT(info);
  *info = ((*n < 0) ? -1 : 0);
  if (!*n || (*info < 0))
    return;
  if (PVN_IS_VECALIGNED(rx))
    *info |= 1;
  if (PVN_IS_VECALIGNED(ix))
    *info |= 2;
  if (PVN_IS_VECALIGNED(ry))
    *info |= 4;
  if (PVN_IS_VECALIGNED(iy))
    *info |= 8;
  if (PVN_IS_VECALIGNED(rz))
    *info |= 16;
  if (PVN_IS_VECALIGNED(iz))
    *info |= 32;
  const size_t m = (size_t)*n;
  register const __m256i gx = (*incx ? _mm256_set_epi32(*incx * 7, *incx * 6, *incx * 5, *incx * 4, *incx * 3, *incx * 2, *incx, 0) : _mm256_setzero_si256());
  register const __m256i gy = (*incy ? _mm256_set_epi32(*incy * 7, *incy * 6, *incy * 5, *incy * 4, *incy * 3, *incy * 2, *incy, 0) : _mm256_setzero_si256());
  register const __m256i sz = (*incz ? _mm256_set_epi32(*incz * 7, *incz * 6, *incz * 5, *incz * 4, *incz * 3, *incz * 2, *incz, 0) : _mm256_setzero_si256());
  for (size_t i = 0u, rem = m; i < m; (i += 8u), (rem -= 8u)) {
    register __m256 xr_ /*= _mm256_setzero_ps()*/;
    register __m256 xi_ /*= _mm256_setzero_ps()*/;
    register __m256 yr_ /*= _mm256_setzero_ps()*/;
    register __m256 yi_ /*= _mm256_setsero_ps()*/;
    if (rem < 8u) {
      alignas(32u) float tail[8u] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
      for (size_t j = 0u; j < rem; ++j)
        tail[j] = rx[*incx * (i + j)];
      xr_ = _mm256_load_ps(tail); VSP(xr_);
      for (size_t j = 0u; j < rem; ++j)
        tail[j] = ix[*incx * (i + j)];
      xi_ = _mm256_load_ps(tail); VSP(xi_);
      for (size_t j = 0u; j < rem; ++j)
        tail[j] = ry[*incy * (i + j)];
      yr_ = _mm256_load_ps(tail); VSP(yr_);
      for (size_t j = 0u; j < rem; ++j)
        tail[j] = iy[*incy * (i + j)];
      yi_ = _mm256_load_ps(tail); VSP(yi_);
    }
    else {
      xr_ = ((*incx == 1) ? ((*info & 1) ? _mm256_load_ps(rx + i) : _mm256_loadu_ps(rx + i)) : _mm256_i32gather_ps((rx + (*incx * i)), gx, 4)); VSP(xr_);
      xi_ = ((*incx == 1) ? ((*info & 2) ? _mm256_load_ps(ix + i) : _mm256_loadu_ps(ix + i)) : _mm256_i32gather_ps((ix + (*incx * i)), gx, 4)); VSP(xi_);
      yr_ = ((*incy == 1) ? ((*info & 4) ? _mm256_load_ps(ry + i) : _mm256_loadu_ps(ry + i)) : _mm256_i32gather_ps((ry + (*incy * i)), gy, 4)); VSP(yr_);
      yi_ = ((*incy == 1) ? ((*info & 8) ? _mm256_load_ps(iy + i) : _mm256_loadu_ps(iy + i)) : _mm256_i32gather_ps((iy + (*incy * i)), gy, 4)); VSP(yi_);
    }
    register const VD xr = _mm512_cvtps_pd(xr_); VDP(xr);
    register const VD xi = _mm512_cvtps_pd(xi_); VDP(xi);
    register const VD yr = _mm512_cvtps_pd(yr_); VDP(yr);
    register const VD yi = _mm512_cvtps_pd(yi_); VDP(yi);
    register const VD zr = _mm512_fmsub_pd(xr, yr, _mm512_mul_pd(xi, yi)); VDP(zr);
    register const VD zi = _mm512_fmadd_pd(xr, yi, _mm512_mul_pd(xi, yr)); VDP(zi);
    register const __m256 zr_ = _mm512_cvtpd_ps(zr); VSP(zr_);
    register const __m256 zi_ = _mm512_cvtpd_ps(zi); VSP(zi_);
    if (*incz) {
      if (rem < 8u) {
        alignas(32u) float tail[8u] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
        _mm256_store_ps(tail, zr_);
        for (size_t j = 0u; j < rem; ++j)
          rz[*incz * (i + j)] = tail[j];
        _mm256_store_ps(tail, zi_);
        for (size_t j = 0u; j < rem; ++j)
          iz[*incz * (i + j)] = tail[j];
      }
      else if (*incz == 1) {
        if (*info & 16)
          _mm256_store_ps((rz + i), zr_);
        else
          _mm256_storeu_ps((rz + i), zr_);
        if (*info & 32)
          _mm256_store_ps((iz + i), zi_);
        else
          _mm256_storeu_ps((iz + i), zi_);
      }
      else {
        _mm256_i32scatter_ps((rz + (*incz * i)), sz, zr_, 4);
        _mm256_i32scatter_ps((iz + (*incz * i)), sz, zi_, 4);
      }
    }
  }
}
#endif /* __AVX512VL__ */
#endif /* ?VEC_CMPLX_TEST */
