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
CFLAGS += -DVEC_CMPLX_OS="\"$(OS)\"" -DVEC_CMPLX_ARCH="\"$(ARCH)\"" -DVEC_CMPLX_COMPILER="\"$(COMPILER)\"" -I../../libpvn/src
ifdef PRINTOUT
CFLAGS += -DPRINTOUT=STD$(PRINTOUT)_FILENO
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
	@echo $(MAKE) "[COMPILER=icx|gcc] [COMPILER_PREFIX=...] [COMPILER_SUFFIX=...] [NDEBUG=0|1|2|3|...] [PRINTOUT=OUT|ERR] [PROFILE=...] [all|clean|help]"

libvec_cmplx.a: $(OBJS)
	$(AR) $(ARFLAGS) $@ $(OBJS)

%.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@

%.exe : %.c
	$(CC) $(CFLAGS) -DVEC_CMPLX_TEST="\"dynamic_exe\"" $< -o $@ -L. -lvec_cmplx -L../../libpvn/src -lpvn $(LDFLAGS)

clean:
	-$(RM) $(EXES)
	-$(RM) $(LIBS)
	-$(RM) $(OBJS)
	-$(RM) *.opt.yaml
	-$(RM) *.dSYM
