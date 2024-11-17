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
CFLAGS += -DVEC_CMPLX_OS="\"$(OS)\"" -DVEC_CMPLX_ARCH="\"$(ARCH)\"" -DVEC_CMPLX_COMPILER="\"$(COMPILER)\""
ifndef LIBPVN
LIBPVN=../../libpvn/src
endif # !LIBPVN
CFLAGS += -I. -I$(LIBPVN)
LDFLAGS += -L. -lvec_cmplx -L$(LIBPVN) -lpvn
ifndef SLEEF
SLEEF=/usr/local
endif # !SLEEF
CFLAGS += -I$(SLEEF)/include
ifneq ($(realpath $(SLEEF)/lib64),)
LDFLAGS += -L$(SLEEF)/lib64
else # lib
LDFLAGS += -L$(SLEEF)/lib
endif # ?lib64
LDFLAGS += -lsleefquad -lsleef
ifdef QUADMATH
CFLAGS += -DPVN_QUADMATH="\"$(QUADMATH)\""
LDFLAGS += $(QUADMATH)
endif # QUADMATH
ifeq ($(findstring BSD,$(OS)),BSD)
LDFLAGS += -lexecinfo
else # !BSD
LDFLAGS += -ldl
endif # ?BSD
LDFLAGS += -lm
ifdef PRINTOUT
CFLAGS += -DPVN_PRINTOUT=STD$(PRINTOUT)_FILENO
endif # PRINTOUT
ifdef PROFILE
CFLAGS += -DPVN_PROFILE=$(PROFILE)u -fno-inline -finstrument-functions
endif # PROFILE

SRCS=vec_cmplx.c
OBJS=$(SRCS:.c=.o)
LIBS=libvec_cmplx.a
EXES=$(SRCS:.c=.exe)

.PHONY: all help clean

all: $(LIBS) $(EXES)

help:
	@echo $(MAKE) "[COMPILER=icx|gcc] [COMPILER_PREFIX=...] [COMPILER_SUFFIX=...] [MARCH=...] [NDEBUG=0|1|2|3|...] [SLEEF=...] [LIBPVN=...] [PRINTOUT=OUT|ERR] [PROFILE=...] [all|clean|help]"

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
	-$(RM) *.opt.yaml
	-$(RM) *.dSYM
