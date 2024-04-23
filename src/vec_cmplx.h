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

#endif /* !VEC_CMPLX_H */
