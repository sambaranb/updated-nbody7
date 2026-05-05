ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Obtaining the code
#
# Two modes: local-rsync (for private development, currently active) and
# github-tarball (for public release, currently commented out). Switch by
# uncommenting the tarball rules and removing the sync-src rule when going
# public.
VERSION = 9b09a1c4efe8d5cbf348cb4232fd1dac0baef095

# --- Local-rsync mode (private development) ---------------------------------
# Sync from the developer's working tree on every build so local edits on the
# mpi-parallel branch are always picked up. `sync-src' is .PHONY so it always
# runs; the order-only dependency below ensures the sync happens before the
# inner build but does not cause unnecessary rebuilds based on timestamps.
#LOCAL_NBODY7_TREE = /Users/sambaran/updated-nbody7
#
#.PHONY: sync-src
#sync-src:
#	@mkdir -p src/updated-nbody7-$(VERSION)
#	rsync -a --delete \
#	    --exclude='.git' \
#	    --exclude='GPU2/lib/.claude' \
#	    --exclude='*.o' --exclude='*.mod' --exclude='*.a' \
#	    --exclude='*.air' --exclude='*.metallib' \
#	    --exclude='default.profraw' \
#	    $(LOCAL_NBODY7_TREE)/ src/updated-nbody7-$(VERSION)/
#
#src/updated-nbody7-$(VERSION): sync-src

# --- GitHub-tarball mode (public release) -----------------------------------
# Uncomment the rules below and remove the sync-src rule above when publishing.

nbody7.tar.gz:
	$(DOWNLOAD) https://github.com/sambaranb/updated-nbody7/archive/$(VERSION).tar.gz >$@

src/updated-nbody7-$(VERSION): nbody7.tar.gz
	tar zxf $<
	touch updated-nbody7-$(VERSION)
	mv updated-nbody7-$(VERSION) src


# Building the code into a static library
#
# v0 scope: CPU-only NBODY7 behind an MPI-IPC AMUSE worker. Only the flags
# actually emitted by support/config.mk are listed. CUDA is left as a
# commented placeholder for a planned CUDA-variant worker. The other flags
# referenced by the AMUSE Makefile template (HDF5, GSL, BLAS/LAPACK,
# SAPPORO_LIGHT, G6LIB, SIMPLE_HASH, STOPCONDMPI, AMUSE_MPI, NETCDF4, MPFR,
# GMP, FFTW, OPENMP_CFLAGS/CXXFLAGS) are not used by NBODY7 and have been
# removed as dead references. Re-add them here (and in support/configure.ac)
# if a future variant needs them.

DEPFLAGS += $(STOPCOND_CFLAGS)
DEPFLAGS += $(FORSOCKETS_CFLAGS)
DEPFLAGS += $(OPENMP_FCFLAGS)

# Where the NBODY7 source headers (common6.h, amuse.h) live after the
# upstream tree has been synced/extracted into src/. interface.f90 uses
# bare `include 'common6.h'` / `include 'amuse.h'`, resolved via this -I.
DEPFLAGS += -Isrc/updated-nbody7-$(VERSION)/Block

# CUDA (for a future nbody7_cuda_worker variant): also uncomment
# AMUSE_CUDA() in support/configure.ac so NVCC, CUDA_FLAGS, CUDA_LDFLAGS get
# emitted into config.mk.
#DEPFLAGS += $(CUDA_FLAGS)

CFLAGS  += $(DEPFLAGS)
FCFLAGS += $(DEPFLAGS)

#LDFLAGS += $(CUDA_LDFLAGS)

LDLIBS += -lm $(STOPCOND_LIBS)
LDLIBS += $(FORSOCKETS_LIBS)

# libnbody7.a contains C++ object files (gpunb.cpp, gpupot.cpp,
# cnbint6.cpp, cxvpred.cpp, wtime.cpp, gpuirr.cpp) and OpenMP-annotated
# Fortran (intgrt.omp.f). The Fortran-driven worker link therefore
# needs the C++ standard library and the OpenMP runtime explicitly.
# -lc++ is the clang/libc++ name on macOS (Apple and conda-clang).
# -lstdc++ is the GCC/libstdc++ name on Linux (Ubuntu, Debian, Fedora, etc.).
ifeq ($(shell uname -s),Darwin)
LDLIBS += -lc++
else
LDLIBS += -lstdc++
endif
LDLIBS += -fopenmp


CODELIB = src/libnbody7.a


##### Modify build command as needed #####
.PHONY: $(CODELIB)
$(CODELIB): | src/updated-nbody7-$(VERSION)
	$(MAKE) -C src -j $(CPU_COUNT) VERSION=$(VERSION) all


# Building the workers
nbody7_worker.f90: interface.py
	amusifier --type=f90 interface.py Nbody7Interface -o $@

nbody7_worker.o: nbody7_worker.f90
	$(MPIFC) -c -o $@ $(FCFLAGS) $<

# interface.f90 includes common6.h and amuse.h from the synced upstream
# tree, so the sync (or tarball-extract) must complete first. Use an
# order-only prerequisite so the dep doesn't cause needless rebuilds.
interface.o: interface.f90 | src/updated-nbody7-$(VERSION)
	$(MPIFC) -ffixed-form -ffixed-line-length-none -o $@ -c $(FCFLAGS) $<

nbody7_worker: interface.o nbody7_worker.o $(CODELIB)
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS)


%.o: %.f90
	$(MPIFC) -o $@ -c $(FCFLAGS) $<


# Which packages contain which workers?
##### Add and/or remove as needed #####
amuse-nbody7_contains: nbody7_worker

amuse-nbody7-SUFFIX_contains: nbody7_SUFFIX_worker


# Building and installing packages
develop-%: %_contains
	support/shared/uninstall.sh $*
	python -m pip install -e packages/$*

install-%: %_contains
	support/shared/uninstall.sh $*
	python -m pip install packages/$*

package-%: %_contains
	python3 -m pip install -vv --no-cache-dir --no-deps --no-build-isolation --prefix ${PREFIX} packages/$*

test-%:
	cd packages/$* && pytest


# Cleaning up
.PHONY: clean
clean:
	$(MAKE) -C src clean
	rm -rf *.o *.mod *worker*

.PHONY: distclean
distclean: clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache
