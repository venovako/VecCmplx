AR=xiar
ARFLAGS=-qnoipo -lib rsv
CC=$(COMPILER_PREFIX)icx$(COMPILER_SUFFIX)
ifdef NDEBUG
CFLAGS=-DNDEBUG=$(NDEBUG) -O$(NDEBUG) -fno-math-errno -qopt-report=3
ifndef PROFILE
CFLAGS += -inline-level=2
endif # !PROFILE
else # DEBUG
CFLAGS=-O0 -g -debug extended -debug inline-debug-info -debug pubnames -debug parallel -ftrapv
endif # ?NDEBUG
ifndef CPU
CPU=Host
# common-avx512 for KNLs
endif # !CPU
CFLAGS += -D_GNU_SOURCE -D_LARGEFILE64_SOURCE -std=gnu18 -fPIC -fexceptions -fasynchronous-unwind-tables -fp-model=precise -fp-speculation=safe -fimf-precision=high -fprotect-parens -fma -no-ftz -qopenmp -fno-omit-frame-pointer -mprefer-vector-width=512 -traceback -vec-threshold0 -x$(CPU)
LDFLAGS=-rdynamic -static-libgcc -ldl
ifndef QUADMATH
QUADMATH=$(shell gcc -print-file-name=libquadmath.a)
ifeq ($(QUADMATH),libquadmath.a)
QUADMATH=-lquadmath
endif # ?QUADMATH
endif # !QUADMATH
CFLAGS += -DPVN_QUADMATH="\"$(QUADMATH) -limf\""
LDFLAGS += $(QUADMATH) -lm
