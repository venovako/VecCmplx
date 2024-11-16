#ifndef VEC_CMPLX_H
#define VEC_CMPLX_H

#include "pvn.h"

#ifdef __AVX512F__
#include <immintrin.h>
#else /* !__AVX512F__ */
#error AVX-512 instructions not available
#endif /* ?__AVX512F__ */

/* default MXCSR */

#ifdef STD_MXCSR
#error STD_MXCSR already defined
#else /* !STD_MXCSR */
#define STD_MXCSR 0x1F80u
#endif /* ?STD_MXCSR */

#ifdef IS_STD_MXCSR
#error IS_STD_MXCSR already defined
#else /* !IS_STD_MXCSR */
#define IS_STD_MXCSR ((_mm_getcsr() & 0xFFC0u) == (STD_MXCSR))
#endif /* ?IS_STD_MXCSR */

/* vector types */

/* vector type containing integers */
#ifdef VI
#error VI already defined
#else /* !VI */
#define VI __m512i
#endif /* ?VI */

/* vector type containing floats */
#ifdef VS
#error VS already defined
#else /* !VS */
#define VS __m512
#endif /* ?VS */

/* vector type containing doubles */
#ifdef VD
#error VD already defined
#else /* !VD */
#define VD __m512d
#endif /* ?VD */

/* mask types */

/* mask type for float lanes */
#ifdef MS
#error MS already defined
#else /* !MS */
#define MS __mmask16
#endif /* ?MS */

/* mask type for double lanes */
#ifdef MD
#error MD already defined
#else /* !MD */
#define MD __mmask8
#endif /* ?MD */

/* vector length in 32-bit lanes */
#ifdef VSL
#error VSL already defined
#else /* !VSL */
#define VSL 16u
#endif /* ?VSL */

/* vector length in 64-bit lanes */
#ifdef VDL
#error VDL already defined
#else /* !VDL */
#define VDL 8u
#endif /* ?VDL */

/* mask operations */

#ifdef MS2U
#error MS2U already defined
#else /* !MS2U */
#define MS2U(m) _cvtmask16_u32(m)
#endif /* ?MS2U */
#ifdef MD2U
#error MD2U already defined
#else /* !MD2U */
#ifdef __AVX512DQ__
#define MD2U(m) _cvtmask8_u32(m)
#else /* !__AVX512DQ__ */
#define MD2U(m) _cvtmask16_u32(m)
#endif /* ?__AVX512DQ__ */
#endif /* ?MD2U */

/* printout */

PVN_EXTERN_C void VSprintf(const int f, const char *const h, const VS v);
PVN_EXTERN_C void VDprintf(const int f, const char *const h, const VD v);

PVN_EXTERN_C void MSprintf(const int f, const char *const h, const MS m);
PVN_EXTERN_C void MDprintf(const int f, const char *const h, const MD m);

#ifdef VSP
#error VSP already defined
#else /* !VSP */
#ifdef PVN_PRINTOUT
#define VSP(v) VSprintf(PVN_PRINTOUT, #v, (v))
#else /* !PVN_PRINTOUT */
#define VSP(v)
#endif /* ?PVN_PRINTOUT */
#endif /* ?VSP */

#ifdef VDP
#error VDP already defined
#else /* !VDP */
#ifdef PVN_PRINTOUT
#define VDP(v) VDprintf(PVN_PRINTOUT, #v, (v))
#else /* !PVN_PRINTOUT */
#define VDP(v)
#endif /* ?PVN_PRINTOUT */
#endif /* ?VDP */

#ifdef MSP
#error MSP already defined
#else /* !MSP */
#ifdef PVN_PRINTOUT
#define MSP(m) MSprintf(PVN_PRINTOUT, #m, (m))
#else /* !PVN_PRINTOUT */
#define MSP(m)
#endif /* ?PVN_PRINTOUT */
#endif /* ?MSP */

#ifdef MDP
#error MDP already defined
#else /* !MDP */
#ifdef PVN_PRINTOUT
#define MDP(m) MDprintf(PVN_PRINTOUT, #m, (m))
#else /* !PVN_PRINTOUT */
#define MDP(m)
#endif /* ?PVN_PRINTOUT */
#endif /* ?MDP */

static inline size_t n2VS(const size_t n)
{
  return ((n + (VSL - 1u)) / VSL);
}

static inline size_t n2VD(const size_t n)
{
  return ((n + (VDL - 1u)) / VDL);
}

static inline size_t n2NS(const size_t n)
{
  return (n2VS(n) * VSL);
}

static inline size_t n2ND(const size_t n)
{
  return (n2VD(n) * VDL);
}

/* (Re0,Im0, Re1,Im1, ..., Re7,Im7) */

static inline VS Sxmuly(const VS x, const VS y)
{
  register const VI si = _mm512_set_epi32(15, 13, 11, 9, 7, 5, 3, 1, 14, 12, 10, 8, 6, 4, 2, 0);
  register const VI mi = _mm512_set_epi32(15, 7, 14, 6, 13, 5, 12, 4, 11, 3, 10, 2, 9, 1, 8, 0);
  register const VS xp = _mm512_permutexvar_ps(si, x); VSP(xp);
  register const VS yp = _mm512_permutexvar_ps(si, y); VSP(yp);
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
  return _mm512_permutexvar_ps(mi, _mm512_insertf32x8(_mm512_zextps256_ps512(_mm512_cvtpd_ps(zr)), _mm512_cvtpd_ps(zi), 1));
#else /* !__AVX512DQ__ */
  return _mm512_permutexvar_ps(mi, _mm512_castpd_ps(_mm512_insertf64x4(_mm512_castps_pd(_mm512_zextps256_ps512(_mm512_cvtpd_ps(zr))), _mm256_castps_pd(_mm512_cvtpd_ps(zi)), 1)));
#endif /* ?__AVX512DQ__ */
}

static inline VD Dxmuly(const VD x, const VD y)
{
  /* TODO: use Sleef quad */
  return _mm512_setzero_pd();
}
#endif /* !VEC_CMPLX_H */
