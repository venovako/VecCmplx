# MARCH=common-avx512 for KNLs
AR=ar
ARFLAGS=rsv
include ../../libpvn/src/pvn.mk
CC=$(PVN_CC)
CFLAGS=$(PVN_CFLAGS) -D_GNU_SOURCE -D_LARGEFILE64_SOURCE $(PVN_CPPFLAGS)
LDFLAGS=$(PVN_LDFLAGS) -L. -lvec_cmplx $(PVN_LIBS)
