ARFLAGS=rsv
CC=$(COMPILER_PREFIX)gcc$(COMPILER_SUFFIX)
ifdef NDEBUG
CFLAGS=-DNDEBUG=$(NDEBUG) -O$(NDEBUG) -fno-math-errno -fvect-cost-model=unlimited
else # DEBUG
CFLAGS=-Og -ggdb3 -ftrapv
endif # ?NDEBUG
CFLAGS += -D_LARGEFILE64_SOURCE -Wall -Wextra
ifeq ($(OS),Linux)
CFLAGS += -D_GNU_SOURCE
else # !Linux
CFLAGS += -m64
endif # ?Linux
CFLAGS += -std=gnu18 -fPIC -fexceptions -fasynchronous-unwind-tables -ffp-contract=fast -fno-omit-frame-pointer -fopenmp
ifeq ($(ARCH),ppc64le)
CFLAGS += -mcpu=native -mtraceback=full
else # !ppc64le
CFLAGS += -march=native
endif # ?ppc64le
LDFLAGS=-rdynamic -static-libgcc
ifeq ($(findstring BSD,$(OS)),BSD)
LDFLAGS += -lexecinfo
else # !BSD
LDFLAGS += -ldl
endif # ?BSD
ifeq ($(findstring 86,$(ARCH)),86)
ifndef QUADMATH
QUADMATH=$(shell $(CC) -print-file-name=libquadmath.a)
ifeq ($(QUADMATH),libquadmath.a)
QUADMATH=-lquadmath
endif # ?QUADMATH
endif # !QUADMATH
CFLAGS += -DPVN_QUADMATH="\"$(QUADMATH)\""
LDFLAGS += $(QUADMATH)
endif # ?86
LDFLAGS += -lm
