OS=$(shell uname)
ARCH=$(shell uname -m)
RM=rm -rfv
ifndef COMPILER
ifeq ($(OS),Linux)
COMPILER=icx
else # !Linux
COMPILER=gcc
endif # ?Linux
endif # !COMPILER
include $(COMPILER).mk
CFLAGS += -DVEC_CMPLX_OS="\"$(OS)\"" -DVEC_CMPLX_ARCH="\"$(ARCH)\"" -DVEC_CMPLX_COMPILER="\"$(COMPILER)\"" -DPVN_SLEEF=0

SRCS=vec_cmplx.c
OBJS=$(SRCS:.c=.o)
LIBS=libvec_cmplx.a
EXES=$(SRCS:.c=.exe)

.PHONY: all help clean

all: $(LIBS) $(EXES)

help:
	@echo $(MAKE) "[COMPILER=icx|gcc] [all|clean|help]"

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
