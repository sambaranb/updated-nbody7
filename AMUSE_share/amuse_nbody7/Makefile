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
VERSION = fb9214bb78d39d69f08b0717f299074d0c365f9a

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
# v1 scope: NBODY7 behind an MPI-IPC AMUSE worker, with cpu/sse/avx/gpu
# variants built side-by-side as separate per-variant wheels. The flags
# below are variant-agnostic and apply to interface.f90 / nbody7_worker.f90
# compilation; per-variant differences (CUDA libs for gpu, SIMD compiler
# flags for sse/avx) are scoped via target-specific variables on the
# worker rules below. The other flags referenced by the AMUSE Makefile
# template (HDF5, GSL, BLAS/LAPACK, SAPPORO_LIGHT, G6LIB, SIMPLE_HASH,
# STOPCONDMPI, AMUSE_MPI, NETCDF4, MPFR, GMP, FFTW, OPENMP_CFLAGS/CXXFLAGS)
# are not used by NBODY7 and have been removed as dead references.

DEPFLAGS += $(STOPCOND_CFLAGS)
DEPFLAGS += $(FORSOCKETS_CFLAGS)
DEPFLAGS += $(OPENMP_FCFLAGS)

# Where the NBODY7 source headers (common6.h, amuse.h) live after the
# upstream tree has been synced/extracted into src/. interface.f90 uses
# bare `include 'common6.h'` / `include 'amuse.h'`, resolved via this -I.
DEPFLAGS += -Isrc/updated-nbody7-$(VERSION)/Block

CFLAGS  += $(DEPFLAGS)
FCFLAGS += $(DEPFLAGS)

LDLIBS += -lm $(STOPCOND_LIBS)
LDLIBS += $(FORSOCKETS_LIBS)

# Each libnbody7_<variant>.a contains C++ object files (gpunb*.cpp,
# gpupot*.cpp, cnbint*.cpp, cxvpred.cpp, wtime.cpp, gpuirr*.cpp) and
# OpenMP-annotated Fortran (intgrt.omp.f). The Fortran-driven worker
# link therefore needs the C++ standard library and the OpenMP runtime
# explicitly. -lc++ is the clang/libc++ name on macOS (Apple and
# conda-clang). -lstdc++ is the GCC/libstdc++ name on Linux.
ifeq ($(shell uname -s),Darwin)
LDLIBS += -lc++
else
LDLIBS += -lstdc++
endif
LDLIBS += -fopenmp


# Per-variant static libraries built by src/Makefile. Each worker links
# against exactly one of these. .PHONY so that the inner make is always
# re-entered (it owns its own dependency tracking against the upstream
# tree, and the wildcard archive recipe must re-evaluate each time).
#
# These are explicit rules rather than a single src/libnbody7_%.a:
# pattern rule because GNU Make disables pattern-rule matching for any
# target declared .PHONY. Using a pattern rule with .PHONY would silently
# skip the codelib build and leave the worker link to fail on a missing
# archive.
.PHONY: src/libnbody7_cpu.a src/libnbody7_sse.a src/libnbody7_avx.a src/libnbody7_gpu.a

# The explicit `all` goal is required: src/Makefile (and ../support/config.mk
# below it) only includes config.mk when MAKECMDGOALS contains something
# other than clean/distclean. Omitting `all` leaves MAKECMDGOALS empty in
# the sub-make, which silently skips the include and produces empty
# MPIFC / FC / CXX / AR / RANLIB and a broken upstream invocation.
src/libnbody7_cpu.a: | src/updated-nbody7-$(VERSION)
	$(MAKE) -C src -j $(CPU_COUNT) VERSION=$(VERSION) VARIANT=cpu all

src/libnbody7_sse.a: | src/updated-nbody7-$(VERSION)
	$(MAKE) -C src -j $(CPU_COUNT) VERSION=$(VERSION) VARIANT=sse all

src/libnbody7_avx.a: | src/updated-nbody7-$(VERSION)
	$(MAKE) -C src -j $(CPU_COUNT) VERSION=$(VERSION) VARIANT=avx all

src/libnbody7_gpu.a: | src/updated-nbody7-$(VERSION)
	$(MAKE) -C src -j $(CPU_COUNT) VERSION=$(VERSION) VARIANT=gpu all


# Building the workers
#
# nbody7_worker.f90 and interface.f90 are variant-agnostic: the same
# Fortran glue and the same Nbody7Interface class signatures work for
# every variant. The variant difference is purely which compiled NBODY7
# routines (cpu/sse/avx/gpu) the worker contains, set at link time by
# linking against a different libnbody7_<variant>.a.
nbody7_worker.f90: interface.py
	amusifier --type=f90 interface.py Nbody7Interface -o $@

nbody7_worker.o: nbody7_worker.f90
	$(MPIFC) -c -o $@ $(FCFLAGS) $<

# interface.f90 includes common6.h and amuse.h from the synced upstream
# tree, so the sync (or tarball-extract) must complete first. Use an
# order-only prerequisite so the dep doesn't cause needless rebuilds.
interface.o: interface.f90 | src/updated-nbody7-$(VERSION)
	$(MPIFC) -ffixed-form -ffixed-line-length-none -o $@ -c $(FCFLAGS) $<

# Per-variant linked workers. Each spawns a separate Python-level
# selection via Nbody7Interface(mode=MODE_<variant>) which dispatches
# name_of_the_worker to the matching binary at runtime.
nbody7_worker:     interface.o nbody7_worker.o src/libnbody7_cpu.a
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS)

nbody7_sse_worker: interface.o nbody7_worker.o src/libnbody7_sse.a
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS)

nbody7_avx_worker: interface.o nbody7_worker.o src/libnbody7_avx.a
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS)

# GPU worker additionally needs the CUDA runtime; CUDA_LDFLAGS is emitted
# by support/configure.ac when AMUSE_CUDA() is active.
nbody7_gpu_worker: LDLIBS += $(CUDA_LDFLAGS)
nbody7_gpu_worker: interface.o nbody7_worker.o src/libnbody7_gpu.a
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS)


%.o: %.f90
	$(MPIFC) -o $@ -c $(FCFLAGS) $<


# Which packages contain which workers?
amuse-nbody7_contains:     nbody7_worker
amuse-nbody7-sse_contains: nbody7_sse_worker
amuse-nbody7-avx_contains: nbody7_avx_worker
amuse-nbody7-gpu_contains: nbody7_gpu_worker


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
