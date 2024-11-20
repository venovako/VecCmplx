AR=ar
ARFLAGS=rsv
CC=$(COMPILER_PREFIX)icx$(COMPILER_SUFFIX)
CFLAGS=-fno-math-errno
ifdef NDEBUG
CFLAGS += -DNDEBUG=$(NDEBUG) -O$(NDEBUG) -qopt-report=3
ifndef PROFILE
CFLAGS += -inline-level=2
endif # !PROFILE
else # DEBUG
CFLAGS += -O0 -g -debug extended -debug inline-debug-info -debug pubnames -ftrapv #-debug parallel
endif # ?NDEBUG
ifndef MARCH
MARCH=Host
# common-avx512 for KNLs
endif # !MARCH
CFLAGS += -D_GNU_SOURCE -D_LARGEFILE64_SOURCE -std=gnu18 -fPIC -fexceptions -fasynchronous-unwind-tables -fp-model=precise -fp-speculation=safe -fimf-precision=high -fprotect-parens -fma -no-ftz -fno-omit-frame-pointer -mprefer-vector-width=512 -traceback -vec-threshold0 -x$(MARCH) #-qopenmp
LDFLAGS=-rdynamic -static-libgcc
ifndef QUADMATH
QUADMATH=$(shell gcc -print-file-name=libquadmath.a)
ifeq ($(QUADMATH),libquadmath.a)
QUADMATH=-lquadmath
endif # ?QUADMATH
endif # !QUADMATH
