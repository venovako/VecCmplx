OS=$(shell uname)
ARCH=$(shell uname -m)
RM=rm -rfv
AR=ar
ARFLAGS=rsv
ifndef LIBPVN
LIBPVN=../../libpvn
endif # !LIBPVN
include $(LIBPVN)/src/pvn.mk
CC=$(PVN_CC)
CFLAGS=$(PVN_CFLAGS) -DVEC_CMPLX_OS="\"$(OS)\"" -DVEC_CMPLX_ARCH="\"$(ARCH)\"" -DVEC_CMPLX_COMPILER="\"$(CC)\"" $(PVN_CPPFLAGS)
LDFLAGS=$(PVN_LDFLAGS) -L. -lvec_cmplx $(PVN_LIBS)

SRCS=vec_cmplx.c
OBJS=$(SRCS:.c=.o)
LIBS=libvec_cmplx.a
EXES=$(SRCS:.c=.exe)

.PHONY: all help clean

all: $(LIBS) $(EXES)

help:
	@echo $(MAKE) "[all|clean|help]"

libvec_cmplx.a: $(OBJS)
	$(AR) $(ARFLAGS) $@ $(OBJS)

%.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@

%.exe : %.c
	$(CC) $(CFLAGS) -DVEC_CMPLX_TEST="\"dynamic_exe\"" $< -o $@ $(LDFLAGS)

clean:
	-$(RM) $(EXES)
	-$(RM) $(LIBS)
	-$(RM) $(OBJS)
	-$(RM) *.optrpt
	-$(RM) *.dSYM
