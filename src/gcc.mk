ARFLAGS=rsv
CC=$(COMPILER_PREFIX)gcc$(COMPILER_SUFFIX)
ifdef NDEBUG
CFLAGS=-DNDEBUG=$(NDEBUG) -O$(NDEBUG) -fno-math-errno
else # DEBUG
CFLAGS=-Og -ggdb3 -ftrapv
endif # ?NDEBUG
CFLAGS += -D_LARGEFILE64_SOURCE -Wall -Wextra
ifeq ($(OS),Linux)
CFLAGS += -D_GNU_SOURCE
else # !Linux
CFLAGS += -m64
endif # ?Linux
ifndef CPU
CPU=native
# e.g., knl for KNL
endif # !CPU
CFLAGS += -std=gnu$(shell if [ `$(CC) -dumpversion | cut -f1 -d.` -ge 14 ]; then echo 23; else echo 18; fi) -fPIC -fexceptions -fasynchronous-unwind-tables -ffp-contract=fast -fopenmp -fno-omit-frame-pointer -fvect-cost-model=unlimited -march=$(CPU)
LDFLAGS=-rdynamic -static-libgcc
ifeq ($(findstring BSD,$(OS)),BSD)
LDFLAGS += -lexecinfo
else # !BSD
LDFLAGS += -ldl
endif # ?BSD
ifndef QUADMATH
QUADMATH=$(shell $(CC) -print-file-name=libquadmath.a)
ifeq ($(QUADMATH),libquadmath.a)
QUADMATH=-lquadmath
endif # ?QUADMATH
endif # !QUADMATH
CFLAGS += -DPVN_QUADMATH="\"$(QUADMATH)\""
LDFLAGS += $(QUADMATH) -lm
