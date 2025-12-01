# MARCH=knl for KNLs
AR=ar
ARFLAGS=rsv
include ../../libpvn/src/pvn.mk
CC=$(PVN_CC)
CFLAGS=$(PVN_CFLAGS) -D_LARGEFILE64_SOURCE
ifeq ($(OS),Linux)
CFLAGS += -D_GNU_SOURCE
else # !Linux
CFLAGS += -m64
endif # ?Linux
CFLAGS += $(PVN_CPPFLAGS)
LDFLAGS=$(PVN_LDFLAGS) -L. -lvec_cmplx $(PVN_LIBS)
