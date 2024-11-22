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

/* printout */

PVN_EXTERN_C void VSprintf(const int f, const char *const h, const VS v);
PVN_EXTERN_C void VDprintf(const int f, const char *const h, const VD v);

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

/* (Re0,Im0, Re1,Im1, ..., ReN,ImN) */
PVN_EXTERN_C void vec_cmul0_(const ssize_t *const n, const float *const x, const float *const y, float *const z, int *const info);
PVN_EXTERN_C void vec_zmul0_(const ssize_t *const n, const double *const x, const double *const y, double *const z, int *const info);
/* (Re0,Re1, ..., ReN); (Im0,Im1, ..., ImN) */
PVN_EXTERN_C void vec_cmul1_(const ssize_t *const n, const float *const rx, const float *const ix, const ssize_t *const incx, const float *const ry, const float *const iy, const ssize_t *const incy, float *const rz, float *const iz, const ssize_t *const incz, int *const info);
PVN_EXTERN_C void vec_zmul1_(const ssize_t *const n, const double *const rx, const double *const ix, const ssize_t *const incx, const double *const ry, const double *const iy, const ssize_t *const incy, double *const rz, double *const iz, const ssize_t *const incz, int *const info);
#endif /* !VEC_CMPLX_H */
