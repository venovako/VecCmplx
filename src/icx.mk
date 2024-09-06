AR=xiar
ARFLAGS=-qnoipo -lib rsv
CC=$(COMPILER_PREFIX)icx$(COMPILER_SUFFIX)
ifdef NDEBUG
CFLAGS=-DNDEBUG=$(NDEBUG) -O$(NDEBUG) -fno-math-errno -inline-level=2 -qopt-report=3
else # DEBUG
CFLAGS=-O0 -g -debug extended -debug inline-debug-info -debug pubnames -debug parallel -ftrapv
endif # ?NDEBUG
CFLAGS += -D_GNU_SOURCE -D_LARGEFILE64_SOURCE -std=gnu18 -fPIC -fexceptions -fasynchronous-unwind-tables -fp-model=precise -fp-speculation=safe -fimf-precision=high -fprotect-parens -fma -no-ftz -qopenmp -fno-omit-frame-pointer -mprefer-vector-width=512 -traceback -vec-threshold0 -xHost
LDFLAGS=-rdynamic -static-libgcc -ldl
ifndef QUADMATH
QUADMATH=$(shell gcc -print-file-name=libquadmath.a)
ifeq ($(QUADMATH),libquadmath.a)
QUADMATH=-lquadmath
endif # ?QUADMATH
endif # !QUADMATH
CFLAGS += -DPVN_QUADMATH="\"$(QUADMATH) -limf\""
LDFLAGS += $(QUADMATH) -lm
