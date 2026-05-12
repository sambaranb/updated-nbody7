#avx       = enable
chainmode = archain
ksblock   = on

ifeq ($(avx), enable)
MKOPTIONS += avx=enable
endif
ifeq ($(chainmode), archain)
MKOPTIONS += chainmode=archain
endif
ifeq ($(ksblock), on)
MKOPTIONS += ksblock=on
endif

ifeq ($(chainmode), archain)
 ifeq ($(ksblock), on)
  BINNAME = nbody7b
 else
  BINNAME = nbody7
 endif
else
 ifeq ($(ksblock), on)
  BINNAME = nbody6b
 else
  BINNAME = nbody6
 endif
endif

FC       = gfortran
FFLAGS   = -O3 -Wall -march=native
CXX      = g++
CXXFLAGS = -O3 -Wall -march=native
ifeq ($(avx), enable)
FFLAGS   += -mavx
CXXFLAGS += -mavx
endif

# CUDA toolkit location.
# Default order:
#   1. CUDA_PATH=... on the make command line (or env)        -- highest priority
#   2. $(CONDA_PREFIX)  (set by `conda activate <env>`; e.g. Amuse-env-gpu)
#   3. /usr/local/cuda  (typical system install)
# AMUSE worker builds pass CUDA_PATH=$(CUDA_TK) explicitly via src/Makefile,
# where CUDA_TK is detected by support/configure (AMUSE_CUDA macro).
ifeq ($(origin CUDA_PATH), undefined)
ifneq ($(CONDA_PREFIX),)
CUDA_PATH := $(CONDA_PREFIX)
else
CUDA_PATH := /usr/local/cuda
endif
endif

# SDK_PATH points at an NVIDIA cuda-samples checkout (public repo:
# https://github.com/NVIDIA/cuda-samples). Required when the .cu files
# are compiled with -DWITH_CUDA5 (our default), because lib/gpunb.velocity.cu
# and lib/gpupot.gpu.cu then `#include <helper_cuda.h>` from cuda-samples'
# Common/ directory. The conda CUDA toolkit and /usr/local/cuda do NOT
# ship these samples -- clone the repo separately and set SDK_PATH.
# The -I line below is added only when SDK_PATH is non-empty, so cpu/sse/avx
# builds remain unaffected.
#SDK_PATH = /path/to/cuda-samples

# Compute capability auto-detection (first attached GPU). nvidia-smi
# returns "8.6" etc.; strip the dot to feed -arch=sm_86. Falls back to
# 86 (Ampere/A40) when nvidia-smi is unavailable or the flag is too old.
# Override with NVCC_ARCH=80 (etc.) on the make command line.
ifeq ($(origin NVCC_ARCH), undefined)
NVCC_ARCH := $(shell nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -1 | tr -d '. ')
endif
ifeq ($(NVCC_ARCH),)
NVCC_ARCH := 86
endif

NVCC = $(CUDA_PATH)/bin/nvcc
NVCC += -arch=sm_$(NVCC_ARCH) -Xptxas -v
ifneq ($(SDK_PATH),)
NVCC += -I $(SDK_PATH)/Common
endif
ifeq ($(avx), enable)
NVCC += -Xcompiler "-O3 -Wall -fopenmp -mavx"
else
NVCC += -Xcompiler "-O3 -Wall -fopenmp"
endif
NVCC += -DWITH_CUDA5
# NTHREAD=128 (the WITH_MAX_THREADS branch in lib/gpunb.velocity.cu and
# lib/gpupot.gpu.cu) needs enough per-SM headroom to host 128 threads per
# block at those kernels' register-pressure footprint. Empirically, sm_86
# (Ampere A40, RTX 30-series) is fine; sm_75 (Turing, RTX 2080) spills
# and trips cudaErrorIllegalAddress at the next device->host memcpy.
# Gate the define on detected compute capability; below sm_86 we fall
# back to the upstream-safe NTHREAD=64. Loosen the threshold to 80 only
# after confirming on a specific Ampere data-center card (e.g. A100).
ifeq ($(shell test $(NVCC_ARCH) -ge 86 && echo y),y)
NVCC += -DWITH_MAX_THREADS
endif

RUNDIR = ./run_versions

# BUILD can be overridden so per-variant invocations emit into separate
# directories (Build_cpu, Build_sse, Build_avx, Build_gpu). Used by the
# AMUSE worker build to archive multiple variants from the same source
# tree without object-file collisions.
BUILD = Build

ifeq ($(avx), enable)
all: gpu avx
else
all: gpu sse
endif

$(BUILD):
	mkdir -p $(BUILD)

$(RUNDIR):
	mkdir -p $(RUNDIR)

gpu: gpunb.gpu.o gpupot.gpu.o intgrt.o gpucor.o nbintp.o cnbint.o nbint.o start.o adjust.o energy2.o phicor.o cmfirr.o cmfirr2.o checkl2.o kspert.o swap.o scale.o wtime.o cxvpred.o gpuirr.o jpred.o jpred2.o fpcorr2.o repair.o sweep2.o ksres3.o check3.o bhplot.o fpoly0.o kspinit.o kspreg.o | $(BUILD) $(RUNDIR)
	cp -f gpunb.gpu.o gpunb.o
	cp -f gpupot.gpu.o gpupot.o
	cp -f *.o ./$(BUILD)
	make gpu -C ./$(BUILD) -f $(shell pwd)/Makefile_cpu.build -j4 $(MKOPTIONS) FC=$(FC) CUDA_PATH=$(CUDA_PATH)
	mv -f ./$(BUILD)/nbody6 $(RUNDIR)/$(BINNAME).gpu

avx: gpunb.avx.o gpupot.avx.o intgrt.o gpucor.o nbintp.o cnbint.o nbint.o start.o adjust.o energy2.o phicor.o cmfirr.o cmfirr2.o checkl2.o kspert.o swap.o scale.o wtime.o cxvpred.o gpuirr.o jpred.o jpred2.o fpcorr2.o repair.o sweep2.o ksres3.o check3.o bhplot.o fpoly0.o kspinit.o kspreg.o | $(BUILD) $(RUNDIR)
	cp -f gpunb.avx.o gpunb.o
	cp -f gpupot.avx.o gpupot.o
	cp -f *.o ./$(BUILD)
	make sse -C ./$(BUILD) -f $(shell pwd)/Makefile_cpu.build -j4 $(MKOPTIONS) FC=$(FC)
	mv -f ./$(BUILD)/nbody6 $(RUNDIR)/$(BINNAME).avx

sse: gpunb.sse.o gpupot.sse.o intgrt.o gpucor.o nbintp.o cnbint.o nbint.o start.o adjust.o energy2.o phicor.o cmfirr.o cmfirr2.o checkl2.o kspert.o swap.o scale.o wtime.o cxvpred.o gpuirr.o jpred.o jpred2.o fpcorr2.o repair.o sweep2.o ksres3.o check3.o bhplot.o fpoly0.o kspinit.o kspreg.o | $(BUILD) $(RUNDIR)
	cp -f gpunb.sse.o gpunb.o
	cp -f gpupot.sse.o gpupot.o
	cp -f *.o ./$(BUILD)
	make sse -C ./$(BUILD) -f $(shell pwd)/Makefile_cpu.build -j4 $(MKOPTIONS) FC=$(FC)
	mv -f ./$(BUILD)/nbody6 $(RUNDIR)/$(BINNAME).sse

cpu: gpunb.cpu.o gpupot.cpu.o intgrt.o gpucor.o nbintp.o cnbint6.o nbint.o start.o adjust.o energy2.o phicor.o cmfirr.o cmfirr2.o checkl2.o kspert.o swap.o scale.o wtime.o cxvpred.o gpuirr.cpu.o jpred.o jpred2.o fpcorr2.o repair.o sweep2.o ksres3.o check3.o bhplot.o fpoly0.o kspinit.o kspreg.o | $(BUILD) $(RUNDIR)
	cp -f gpunb.cpu.o gpunb.o
	cp -f gpupot.cpu.o gpupot.o
	cp -p cnbint6.o cnbint.o
	cp -p gpuirr.cpu.o gpuirr.o
	cp -f *.o ./$(BUILD)
	make cpu -C ./$(BUILD) -f $(shell pwd)/Makefile_cpu.build -j4 $(MKOPTIONS) FC=$(FC)
	mv -f ./$(BUILD)/nbody6 $(RUNDIR)/$(BINNAME).cpu

# Same Stage-1 + Stage-2 object compilation as `cpu`, but stops before
# the throwaway link of `nbody6`. Used by the AMUSE-nbody7 worker build
# which archives the resulting Build/*.o into libnbody7.a and supplies
# its own main. Keeps the AMUSE path independent of platform-specific
# C++ stdlib / system-SDK linker quirks (e.g. -lstdc++ vs -lc++).
cpu-objects: gpunb.cpu.o gpupot.cpu.o intgrt.o gpucor.o nbintp.o cnbint6.o nbint.o start.o adjust.o energy2.o phicor.o cmfirr.o cmfirr2.o checkl2.o kspert.o swap.o scale.o wtime.o cxvpred.o gpuirr.cpu.o jpred.o jpred2.o fpcorr2.o repair.o sweep2.o ksres3.o check3.o bhplot.o fpoly0.o kspinit.o kspreg.o | $(BUILD)
	cp -f gpunb.cpu.o gpunb.o
	cp -f gpupot.cpu.o gpupot.o
	cp -p cnbint6.o cnbint.o
	cp -p gpuirr.cpu.o gpuirr.o
	cp -f *.o ./$(BUILD)
	make cpu-objects -C ./$(BUILD) -f $(shell pwd)/Makefile_cpu.build -j4 $(MKOPTIONS) FC=$(FC)

# Same Stage-1 + Stage-2 object compilation as `sse`/`avx`/`gpu`, but
# stops before the throwaway link of `nbody6`. Used by the AMUSE-nbody7
# variant worker builds (amuse-nbody7-sse / -avx / -gpu) which archive
# the resulting Build_<variant>/*.o into libnbody7_<variant>.a and
# supply their own main via interface.f90. Stage-2 reuses the existing
# cpu-objects target (variant-agnostic: same OBJECTS list, no link).
sse-objects: gpunb.sse.o gpupot.sse.o intgrt.o gpucor.o nbintp.o cnbint.o nbint.o start.o adjust.o energy2.o phicor.o cmfirr.o cmfirr2.o checkl2.o kspert.o swap.o scale.o wtime.o cxvpred.o gpuirr.o jpred.o jpred2.o fpcorr2.o repair.o sweep2.o ksres3.o check3.o bhplot.o fpoly0.o kspinit.o kspreg.o | $(BUILD)
	cp -f gpunb.sse.o gpunb.o
	cp -f gpupot.sse.o gpupot.o
	cp -f *.o ./$(BUILD)
	make cpu-objects -C ./$(BUILD) -f $(shell pwd)/Makefile_cpu.build -j4 $(MKOPTIONS) FC=$(FC)

avx-objects: gpunb.avx.o gpupot.avx.o intgrt.o gpucor.o nbintp.o cnbint.o nbint.o start.o adjust.o energy2.o phicor.o cmfirr.o cmfirr2.o checkl2.o kspert.o swap.o scale.o wtime.o cxvpred.o gpuirr.o jpred.o jpred2.o fpcorr2.o repair.o sweep2.o ksres3.o check3.o bhplot.o fpoly0.o kspinit.o kspreg.o | $(BUILD)
	cp -f gpunb.avx.o gpunb.o
	cp -f gpupot.avx.o gpupot.o
	cp -f *.o ./$(BUILD)
	make cpu-objects -C ./$(BUILD) -f $(shell pwd)/Makefile_cpu.build -j4 $(MKOPTIONS) FC=$(FC)

gpu-objects: gpunb.gpu.o gpupot.gpu.o intgrt.o gpucor.o nbintp.o cnbint.o nbint.o start.o adjust.o energy2.o phicor.o cmfirr.o cmfirr2.o checkl2.o kspert.o swap.o scale.o wtime.o cxvpred.o gpuirr.o jpred.o jpred2.o fpcorr2.o repair.o sweep2.o ksres3.o check3.o bhplot.o fpoly0.o kspinit.o kspreg.o | $(BUILD)
	cp -f gpunb.gpu.o gpunb.o
	cp -f gpupot.gpu.o gpupot.o
	cp -f *.o ./$(BUILD)
	make cpu-objects -C ./$(BUILD) -f $(shell pwd)/Makefile_cpu.build -j4 $(MKOPTIONS) FC=$(FC) CUDA_PATH=$(CUDA_PATH)

metal: gpunb.metallib gpunb.metal.o gpupot.metallib gpupot.metal.o intgrt.o gpucor.o nbintp.o cnbint6.o nbint.o start.o adjust.o energy2.o phicor.o cmfirr.o cmfirr2.o checkl2.o kspert.o swap.o scale.o wtime.o cxvpred.o gpuirr.cpu.o jpred.o jpred2.o fpcorr2.o repair.o sweep2.o ksres3.o check3.o bhplot.o fpoly0.o kspinit.o kspreg.o | $(BUILD) $(RUNDIR)
	cp -f gpunb.metal.o gpunb.o
	cp -f gpupot.metal.o gpupot.o
	cp -p cnbint6.o cnbint.o
	cp -p gpuirr.cpu.o gpuirr.o
	cp -f *.o ./$(BUILD)
	cp -f *.metallib ./$(BUILD)
	make metal -C ./$(BUILD) -f $(shell pwd)/Makefile_cpu.build -j4 $(MKOPTIONS) FC=$(FC)
	mv -f ./$(BUILD)/nbody6 $(RUNDIR)/$(BINNAME).metal

cnbint.o: ./lib/cnbint.cpp
	$(CXX) $(CXXFLAGS) -c $^

cnbint6.o: ./lib/cnbint6.cpp
	$(CXX) $(CXXFLAGS) -c $^

wtime.o: ./lib/wtime.cpp
	$(CXX) $(CXXFLAGS) -c $^

cxvpred.o: ./lib/cxvpred.cpp
	$(CXX) $(CXXFLAGS) -fopenmp -c $^

ifeq ($(avx), enable)
gpuirr.o: gpuirr.s
	$(CXX) -c $^ -o $@

gpuirr.s: gpuirr.avx.s
	grep -E -v 'v?movaps\s(%xmm[0-9]+),\s?\1$$' < $^ > $@ #Thanks, Sakuraba

gpuirr.avx.s: irrlib/gpuirr.avx.cpp
	$(CXX) $(CXXFLAGS) -mavx -fopenmp -S $^ -o $@

else
gpuirr.o: ./irrlib/gpuirr.sse.cpp
	$(CXX) $(CXXFLAGS) -fopenmp -c $^ -o $@

endif
intgrt.o: intgrt.omp.f
	$(FC) $(FFLAGS) -fopenmp $^ -c -o $@ 

intgrt.s: intgrt.omp.f
	$(FC) $(FFLAGS) -fopenmp $^ -S -o $@

start.o: start.f
	$(FC) $(FFLAGS) -fopenmp $^ -c -o $@

gpunb.gpu.o: lib/gpunb.velocity.cu
	$(NVCC) $^ -c -o $@

gpunb.velocity.su: gpunb.velocity.cubin
	cuobjdump -sass $< | c++filt > $@

gpunb.velocity.cubin: lib/gpunb.velocity.cu
	$(NVCC) -cubin --ptxas-options=-v $<

gpupot.gpu.o: lib/gpupot.gpu.cu
	$(NVCC) $^ -c -o $@

gpunb.avx.o: lib/gpunb.avx.cpp
	$(CXX) $(CXXFLAGS) -fopenmp $^ -c -o $@

gpupot.avx.o: lib/gpupot.avx.cpp
	$(CXX) $(CXXFLAGS) -fopenmp $^ -c -o $@

gpunb.sse.o: lib/gpunb.sse.cpp
	$(CXX) $(CXXFLAGS) -fopenmp $^ -c -o $@

gpupot.sse.o: lib/gpupot.sse.cpp
	$(CXX) $(CXXFLAGS) -fopenmp $^ -c -o $@

gpuirr.cpu.o: ./irrlib/gpuirr.cpp
	$(CXX) $(CXXFLAGS) -fopenmp -c $^ -o $@

gpunb.cpu.o: lib/gpunb.cpp
	$(CXX) $(CXXFLAGS) -DWITH_MAX_CPU_THREADS -fopenmp $^ -c -o $@

gpupot.cpu.o: lib/gpupot.cpp
	$(CXX) $(CXXFLAGS) -fopenmp $^ -c -o $@

# Metal shader and host compilation (requires Apple clang)
METALCC = xcrun -sdk macosx metal
OBJCXX  = clang++

gpupot.metallib: lib/gpupot.metal
	$(METALCC) -c $^ -o gpupot.air
	xcrun -sdk macosx metallib gpupot.air -o $@

gpupot.metal.o: lib/gpupot.metal.mm
	$(OBJCXX) -O3 -Wall -fobjc-arc -c $^ -o $@

gpunb.metallib: lib/gpunb.metal
	$(METALCC) -c $^ -o gpunb.air
	xcrun -sdk macosx metallib gpunb.air -o $@

gpunb.metal.o: lib/gpunb.metal.mm
	$(OBJCXX) -O3 -Wall -fobjc-arc -c $^ -o $@

clean:
	rm -f *.o *.s *.air *.metallib ./$(BUILD)/*.o

