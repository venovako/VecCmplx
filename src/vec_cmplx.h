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

/* half-length vector type containing integers */
#ifdef VI_2
#error VI_2 already defined
#else /* !VI_2 */
#define VI_2 __m256i
#endif /* ?VI_2 */

/* vector type containing floats */
#ifdef VS
#error VS already defined
#else /* !VS */
#define VS __m512
#endif /* ?VS */

/* half-length vector type containing floats */
#ifdef VS_2
#error VS_2 already defined
#else /* !VS_2 */
#define VS_2 __m256
#endif /* ?VS_2 */

/* vector type containing doubles */
#ifdef VD
#error VD already defined
#else /* !VD */
#define VD __m512d
#endif /* ?VD */

/* half-length vector type containing doubles */
#ifdef VD_2
#error VD_2 already defined
#else /* !VD_2 */
#define VD_2 __m256d
#endif /* ?VD_2 */

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

/* half the maximal vector length */
#ifdef PVN_VECLEN_2
#error PVN_VECLEN_2 already defined
#else /* !PVN_VECLEN_2 */
#define PVN_VECLEN_2 (PVN_VECLEN >> 1u)
#endif /* ?PVN_VECLEN_2 */

/* vector length in 32-bit lanes */
#ifdef VSL
#error VSL already defined
#else /* !VSL */
#define VSL 16u
#endif /* ?VSL */

/* half the vector length in 32-bit lanes */
#ifdef VSL_2
#error VSL_2 already defined
#else /* !VSL_2 */
#define VSL_2 8u
#endif /* ?VSL_2 */

/* vector length in 64-bit lanes */
#ifdef VDL
#error VDL already defined
#else /* !VDL */
#define VDL 8u
#endif /* ?VDL */

/* half the vector length in 64-bit lanes */
#ifdef VDL_2
#error VDL_2 already defined
#else /* !VDL_2 */
#define VDL_2 4u
#endif /* ?VDL_2 */

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
PVN_EXTERN_C void vec_cmul0_(const ssize_t *const n, const float *const x, const float *const y, float *const z, int *const info);
#ifdef __AVX512VL__
/* (Re0, Re1, ...), (Im0, Im1, ...) */
PVN_EXTERN_C void vec_cmul1_(const ssize_t *const n, const float *const rx, const float *const ix, const ssize_t *const incx, const float *const ry, const float *const iy, const ssize_t *const incy, float *const rz, float *const iz, const ssize_t *const incz, int *const info);
#endif /* __AVX512VL__ */
#endif /* !VEC_CMPLX_H */
