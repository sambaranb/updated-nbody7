# The classic standalone NBODY7 executables and the AMUSE_NBODY7 workers are
# derived from the same public updated-nbody7 source code: https://github.com/sambaranb/updated-nbody7
# This ensures that any update to the source code consistently reaches both pathways.

The AMUSE_NBODY7 interface is under testing. While applying it in various
astrophysical studies, reporting issues, and requesting features is highly
encouraged, the standalone NBODY7 executable (i.e., the amusein=0 variant) is,
for now, the most recommended choice for long-running N-body simulations.
The subfolders contain input and output files of short example standalone runs. 


# Good news for Apple Silicon users

The hard-coded X86 SIMD instructions in some of the
C++ routines (e.g., gpunb, gpupot, gpuirr, cnbint) are no more an obstacle!
Pure CPU counterparts are now implemented. Furthermore, Apple Metal implementations (beta) of gpunb
and gpupot have been created for Apple Silicon platforms, as a substitute for NVIDIA/CUDA.

The CPU/Metal variant has been reasonably tested (by Banerjee so far) against the standard SSE/CUDA variant.
See the example folder test_metal_version/ for such a comparison plot. 
Upon request, Banerjee is happy to provide a Jupyter Notebook for more comparisons.
(The data volume would make cloning inconvenient.)
You are encouraged to do more and more tests and find potential issues.
For now, the Metal implementation should be considered as beta.


# Download and compile standalone NBODY7

While the AMUSE installer downloads the same source code along with the multi-target
Makefile, for a standalone use, it is perhaps cleaner to clone:

git clone -b mpi-parallel --single-branch https://github.com/sambaranb/updated-nbody7
(You are cloning with the mpi-parallel branch, which is currently the most updated
branch. It will probably be merged with the main branch in the near future.)

cd updated-nbody7
(This is the NBODY7 checkout directory)

cd GPU2/

# All platforms: the 'pure' CPU variant
make cpu

# X86_only: X86 SIMD variants
make sse
make avx
(Remember to uncomment the avx=enable line in the beginning of GPU2/Makefile
for the AVX compilation; alternatively, supply avx=enable from the command line.)

# X86 + NVIDIA/CUDA only: the SSE/CUDA variant
make gpu SDK_PATH=</path/to/cuda-samples>  

Note: the cuda-samples toolkit needs to be downloaded separately.
You can simply keep it in your home directory or wherever
convenient; no system-wide installation is necessary.
Also, you don't have to install the example packages, just
unpacking would suffice to access the Common/helper_cuda.h header file.

In the GPU2/Makefile, the CUDA_PATH is defaulted to /usr/local/cuda,
which is most common. Change it to the appropriate value for your
system if necessary. Alternatively, you can override the default
value by supplying CUDA_PATH from the command line.

# Apple Silicon only: the CPU/Metal variant
make metal

# Notes:

The above instructions assume that the Fortran and g++ compilers,
as well as CUDA, are installed in the system, outside of a conda
environment. In case of conda-only installations, activate
the appropriate virtual environment first, both for compiling and running the
executable.

By default, the successfully compiled executables will be located
in GPU2/run_versions with names nbody7b.<variant>. Change the value
of RUNDIR in GPU2/Makefile as necessary if you prefer a different
executable path. 

To avoid potential variant contamination, it is advisable to do a 'make clean'
before compiling each variant. Alternatively, you can supply a different BUILD=
from the command line during each make (as is done automatically during the AMUSE worker
compilation). See GPU2/Makefile.
