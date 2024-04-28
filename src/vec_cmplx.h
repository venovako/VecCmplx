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

/* vector bitwise operations */

#ifdef VSAND
#error VSAND already defined
#endif /* VSAND */
#ifdef VDAND
#error VDAND already defined
#endif /* VDAND */

#ifdef VSANDNOT
#error VSANDNOT already defined
#endif /* VSANDNOT */
#ifdef VDANDNOT
#error VDANDNOT already defined
#endif /* VDANDNOT */

#ifdef VSOR
#error VSOR already defined
#endif /* VSOR */
#ifdef VDOR
#error VDOR already defined
#endif /* VDOR */

#ifdef VSXOR
#error VSXOR already defined
#endif /* VSXOR */
#ifdef VDXOR
#error VDXOR already defined
#endif /* VDXOR */

#ifdef __AVX512DQ__
#define VSAND(x,y) _mm512_and_ps((x),(y))
#define VSANDNOT(x,y) _mm512_andnot_ps((x),(y))
#define VSOR(x,y) _mm512_or_ps((x),(y))
#define VSXOR(x,y) _mm512_xor_ps((x),(y))

#define VDAND(x,y) _mm512_and_pd((x),(y))
#define VDANDNOT(x,y) _mm512_andnot_pd((x),(y))
#define VDOR(x,y) _mm512_or_pd((x),(y))
#define VDXOR(x,y) _mm512_xor_pd((x),(y))
#else /* !__AVX512DQ__ */
#define VSAND(x,y) _mm512_castsi512_ps(_mm512_and_epi32(_mm512_castps_si512(x),_mm512_castps_si512(y)))
#define VSANDNOT(x,y) _mm512_castsi512_ps(_mm512_andnot_epi32(_mm512_castps_si512(x),_mm512_castps_si512(y)))
#define VSOR(x,y) _mm512_castsi512_ps(_mm512_or_epi32(_mm512_castps_si512(x),_mm512_castps_si512(y)))
#define VSXOR(x,y) _mm512_castsi512_ps(_mm512_xor_epi32(_mm512_castps_si512(x),_mm512_castps_si512(y)))

#define VDAND(x,y) _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(x),_mm512_castpd_si512(y)))
#define VDANDNOT(x,y) _mm512_castsi512_pd(_mm512_andnot_epi64(_mm512_castpd_si512(x),_mm512_castpd_si512(y)))
#define VDOR(x,y) _mm512_castsi512_pd(_mm512_or_epi64(_mm512_castpd_si512(x),_mm512_castpd_si512(y)))
#define VDXOR(x,y) _mm512_castsi512_pd(_mm512_xor_epi64(_mm512_castpd_si512(x),_mm512_castpd_si512(y)))
#endif /* ?__AVX512DQ__ */

/* unary arithmetic operations */

#ifdef VSABS
#error VSABS already defined
#else /* !VSABS */
#define VSABS(x) VSANDNOT(_mm512_set1_ps(-0.0f),(x))
#endif /* ?VSABS */
#ifdef VSNEG
#error VSNEG already defined
#else /* !VSNEG */
#define VSNEG(x) VSXOR((x),_mm512_set1_ps(-0.0f))
#endif /* ?VSNEG */
#ifdef VSSGN
#error VSSGN already defined
#else /* !VSSGN */
#define VSSGN(x) VSAND((x),_mm512_set1_ps(-0.0f))
#endif /* ?VSSGN */

#ifdef VDABS
#error VDABS already defined
#else /* !VDABS */
#define VDABS(x) VDANDNOT(_mm512_set1_pd(-0.0),(x))
#endif /* ?VDABS */
#ifdef VDNEG
#error VDNEG already defined
#else /* !VDNEG */
#define VDNEG(x) VDXOR((x),_mm512_set1_pd(-0.0))
#endif /* ?VDNEG */
#ifdef VDSGN
#error VDSGN already defined
#else /* !VDSGN */
#define VDSGN(x) VDAND((x),_mm512_set1_pd(-0.0))
#endif /* ?VDSGN */

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

#ifdef MSOR
#error MSOR already defined
#else /* !MSOR */
#define MSOR(a,b) _kor_mask16((a),(b))
#endif /* ?MSOR */
#ifdef MDOR
#error MDOR already defined
#else /* !MDOR */
#ifdef __AVX512DQ__
#define MDOR(a,b) _kor_mask8((a),(b))
#else /* !__AVX512DQ__ */
#define MDOR(a,b) (__mmask8)_kor_mask16((a),(b))
#endif /* ?__AVX512DQ__ */
#endif /* ?MDOR */

#ifdef MSXOR
#error MSXOR already defined
#else /* !MSXOR */
#define MSXOR(a,b) _kxor_mask16((a),(b))
#endif /* ?MSXOR */
#ifdef MDXOR
#error MDXOR already defined
#else /* !MDXOR */
#ifdef __AVX512DQ__
#define MDXOR(a,b) _kxor_mask8((a),(b))
#else /* !__AVX512DQ__ */
#define MDXOR(a,b) (__mmask8)_kxor_mask16((a),(b))
#endif /* ?__AVX512DQ__ */
#endif /* ?MDXOR */

#ifdef MSAND
#error MSAND already defined
#else /* !MSAND */
#define MSAND(a,b) _kand_mask16((a),(b))
#endif /* ?MSAND */
#ifdef MDAND
#error MDAND already defined
#else /* !MDAND */
#ifdef __AVX512DQ__
#define MDAND(a,b) _kand_mask8((a),(b))
#else /* !__AVX512DQ__ */
#define MDAND(a,b) (__mmask8)_kand_mask16((a),(b))
#endif /* ?__AVX512DQ__ */
#endif /* ?MDAND */

#ifdef MSANDN
#error MSANDN already defined
#else /* !MSANDN */
#define MSANDN(a,b) _kandn_mask8((a),(b))
#endif /* ?MSANDN */
#ifdef MDANDN
#error MDANDN already defined
#else /* !MDANDN */
#ifdef __AVX512DQ__
#define MDANDN(a,b) _kandn_mask8((a),(b))
#else /* !__AVX512DQ__ */
#define MDANDN(a,b) (__mmask8)_kandn_mask16((a),(b))
#endif /* ?__AVX512DQ__ */
#endif /* ?MDANDN */

/* printout */

extern void VSprintf(const int f, const char *const h, const VS v);
extern void VDprintf(const int f, const char *const h, const VD v);

extern void MSprintf(const int f, const char *const h, const MS m);
extern void MDprintf(const int f, const char *const h, const MD m);

#ifdef VSP
#error VSP already defined
#else /* !VSP */
#ifdef PRINTOUT
#define VSP(v) VSprintf((PRINTOUT), #v, (v))
#else /* !PRINTOUT */
#define VSP(v)
#endif /* ?PRINTOUT */
#endif /* ?VSP */

#ifdef VDP
#error VDP already defined
#else /* !VDP */
#ifdef PRINTOUT
#define VDP(v) VDprintf((PRINTOUT), #v, (v))
#else /* !PRINTOUT */
#define VDP(v)
#endif /* ?PRINTOUT */
#endif /* ?VDP */

#ifdef MSP
#error MSP already defined
#else /* !MSP */
#ifdef PRINTOUT
#define MSP(m) MSprintf((PRINTOUT), #m, (m))
#else /* !PRINTOUT */
#define MSP(m)
#endif /* ?PRINTOUT */
#endif /* ?MSP */

#ifdef MDP
#error MDP already defined
#else /* !MDP */
#ifdef PRINTOUT
#define MDP(m) MDprintf((PRINTOUT), #m, (m))
#else /* !PRINTOUT */
#define MDP(m)
#endif /* ?PRINTOUT */
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
#endif /* !VEC_CMPLX_H */
