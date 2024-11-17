ARFLAGS=rsv
CC=$(COMPILER_PREFIX)gcc$(COMPILER_SUFFIX)
CFLAGS=-fno-math-errno
ifdef NDEBUG
CFLAGS += -DNDEBUG=$(NDEBUG) -O$(NDEBUG)
else # DEBUG
CFLAGS += -Og -ggdb3 -ftrapv
endif # ?NDEBUG
CFLAGS += -D_LARGEFILE64_SOURCE -Wall -Wextra
ifeq ($(OS),Linux)
CFLAGS += -D_GNU_SOURCE
else # !Linux
CFLAGS += -m64
endif # ?Linux
ifndef MARCH
MARCH=native
# knl for KNLs
endif # !MARCH
CFLAGS += -std=gnu$(shell if [ `$(CC) -dumpversion | cut -f1 -d.` -ge 14 ]; then echo 23; else echo 18; fi) -fPIC -fexceptions -fasynchronous-unwind-tables -ffp-contract=fast -fopenmp -fno-omit-frame-pointer -fvect-cost-model=unlimited -march=$(MARCH)
LDFLAGS=-rdynamic -static-libgcc
ifndef QUADMATH
QUADMATH=$(shell $(CC) -print-file-name=libquadmath.a)
ifeq ($(QUADMATH),libquadmath.a)
QUADMATH=-lquadmath
endif # ?QUADMATH
endif # !QUADMATH
