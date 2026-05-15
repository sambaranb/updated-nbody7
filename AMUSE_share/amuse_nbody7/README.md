# amuse-nbody7

AMUSE community-code interface for [NBODY7](https://github.com/sambaranb/updated-nbody7),
S. Aarseth's direct N-body integrator with ARCHAIN post-Newtonian chain
regularization. The interface wraps NBODY7 behind an MPI-IPC AMUSE
worker and ships four variants of the integrator side-by-side as
separate wheels.

References:

S. J. Aarseth, Gravitational N-Body Simulations, by Sverre J. Aarseth, pp. 430. ISBN 0521432723. Cambridge, UK: Cambridge University Press, November 2003. (2003) p. 430 [Handbook of numerical alorithms in NBODY6/7]

S. J. Aarseth, MNRAS 422, 841 (2012), arXiv:1202.4688 [astro-ph.SR] [Original NBODY7 paper]

S. Banerjee, K. Belczynski, C. L. Fryer, P. Berczik, J. R. Hurley, R. Spurzem, and L. Wang, A&A 639, A41 (2020), arXiv:1902.07718 [astro-ph.SR] [Various updates to NBODY7, included in this version]

S. Banerjee, MNRAS 500, 3002 (2021), arXiv:2004.07382 [astro-ph.HE] [Various updates to NBODY7, included in this version]


## Variants

| Wheel                  | Worker binary           | Build requirements                  |
|------------------------|-------------------------|-------------------------------------|
| `amuse-nbody7`         | `nbody7_worker`         | Fortran + C++ compiler              |
| `amuse-nbody7-sse`     | `nbody7_sse_worker`     | x86 with SSE2 (any modern Intel/AMD)|
| `amuse-nbody7-avx`     | `nbody7_avx_worker`     | x86 with AVX/AVX2                   |
| `amuse-nbody7-gpu`     | `nbody7_gpu_worker`     | CUDA toolkit + NVIDIA cuda-samples  |

The base `amuse-nbody7` wheel ships the Python interface and the cpu
worker. The three variant wheels (`-sse`, `-avx`, `-gpu`) are
binary-only — they install just their `nbody7_<v>_worker` binary into
the same `amuse_nbody7` namespace and declare `amuse-nbody7` as a
dependency. Once installed, a user picks a variant per simulation:

```python
from amuse.community.nbody7 import Nbody7
sim = Nbody7(mode="gpu")   # or "cpu" (default), "sse", "avx"
```

## Quick start

After extracting the source tarball into AMUSE's `src/` (or
`packages/`):

```bash
make demo-cpu                                  # build + test the cpu variant
make demo-sse                                  # cpu + sse
make demo-avx                                  # cpu + avx
make demo-gpu SDK_PATH=$HOME/src/cuda-samples  # cpu + gpu (see GPU note below)
make demo-all  [SDK_PATH=...]                  # all variants (gpu only if SDK_PATH set)
make demo-help                                 # short help
```

Each `demo-<variant>` target installs the matching wheel (plus the base
wheel for the SIMD/GPU variants) and runs the pytest sweep
(`tests/test_nbody7.py`, 11 tests) against that variant via
`pytest --mode=<variant>`. The base install is idempotent, so the
targets compose cleanly.

**Why `demo-sse / -avx / -gpu` build the cpu worker too.** The base
`amuse-nbody7` wheel ships the shared Python interface *together with*
the cpu worker (`nbody7_worker`) as a wheel artifact — variant wheels
declare `amuse-nbody7` as a `dependencies` entry, so the base install
must complete before any variant can be pip-installed, and that base
install builds `nbody7_worker`. So compiling cpu first under
`demo-gpu` (or `-sse`, `-avx`) is expected behaviour, not a sign of
variant contamination. The pytest sweep that follows runs strictly
against the requested variant's binary via `--mode=<variant>`.

## Manual build + test (à la carte)

If you want to install and test just one variant by hand — for example
to iterate on a single SIMD path or to run pytest with custom flags —
the underlying steps the demo targets wrap are:

```bash
# Base wheel (required, ships the Python interface + cpu worker):
make install-amuse-nbody7

# SIMD variants (each needs the base wheel installed first):
make install-amuse-nbody7-sse
make install-amuse-nbody7-avx

# GPU variant (also needs the base wheel; see SDK_PATH note below):
make install-amuse-nbody7-gpu SDK_PATH=$HOME/src/cuda-samples

# Run the pytest sweep against a specific variant. Choose your own
# workdir so NBODY7's diagnostic files don't land in the source tree.
mkdir -p /tmp/nb7_check && cd /tmp/nb7_check
python -m pytest --verbose --run-slow --mode=cpu /path/to/amuse_nbody7/tests/test_nbody7.py
python -m pytest --verbose --run-slow --mode=sse /path/to/amuse_nbody7/tests/test_nbody7.py
python -m pytest --verbose --run-slow --mode=avx /path/to/amuse_nbody7/tests/test_nbody7.py
python -m pytest --verbose --run-slow --mode=gpu /path/to/amuse_nbody7/tests/test_nbody7.py
```

The `--mode=<v>` flag (handled by `tests/conftest.py`) monkey-patches
`Nbody7Interface` / `Nbody7` to default to that variant, so the same
11-test sweep exercises whichever worker is selected.

## GPU variant: SDK_PATH is required

`lib/gpunb.velocity.cu` and `lib/gpupot.gpu.cu` `#include <helper_cuda.h>`
from NVIDIA's [cuda-samples](https://github.com/NVIDIA/cuda-samples)
repo. The conda CUDA toolkit does **not** ship cuda-samples, so the GPU
build needs you to pass an explicit path:

```bash
git clone https://github.com/NVIDIA/cuda-samples ~/src/cuda-samples
make demo-gpu SDK_PATH=~/src/cuda-samples
```

If `SDK_PATH` is not set, `make demo-gpu` (and `make demo-all` for the
GPU portion) prints an error and exits.

## Worker output and demo workdir

NBODY7 writes diagnostics (`CMPACT`, `NBSTAT`, `OUT3`, `OUT33`,
`fort.41`, `fort.992`, `STOP*`, `run.out`, `err.out`) on fixed unit
numbers in the worker's current directory. To keep the source tree
clean, the demo targets chdir into a per-variant workdir before
invoking pytest:

```
demo_workdir/
├── cpu/   # cpu run artifacts (+ test/ archive sub-dir)
├── sse/
├── avx/
└── gpu/
```

Override the location with `make demo-cpu DEMO_WORKDIR=/scratch/...`.
`make clean` removes `demo_workdir/` along with the build artifacts.

## Numerical-equivalence sign-off (May 2026)

The four variants have been runtime-verified for numerical equivalence
on four hosts (energy conservation `|ΔE/E|` in 1e-6 – 3e-5 over a
128-body Plummer to T=3, IC bit-identical across variants):

| Host                          | CPU | SSE | AVX | GPU |
|-------------------------------|-----|-----|-----|-----|
| NVIDIA A40 (sm_86)            | ✓   | ✓   | ✓   | ✓   |
| NVIDIA RTX 2080 (sm_75)       | ✓   | ✓   | ✓   | ✓   |
| Linux x86_64 (no GPU)         | ✓   | ✓   | ✓   | —   |
| Apple Silicon M4 Pro (macOS)  | ✓   | —   | —   | —   |

(SIMD variants are x86-only; Metal port for Apple Silicon is in
progress upstream but AMUSE has no Metal infrastructure yet.) The AVX
kernel produces `|ΔE/E|` ~5–10× larger than CPU/SSE/GPU but still well
below test tolerance.

## Known notes

* **GPU on Turing (sm_75) and older:** the GPU2 Makefile gates
  `-DWITH_MAX_THREADS` on `NVCC_ARCH >= 86`. Below that, the CUDA
  threads-per-block falls back to `NTHREAD=64`; without the gate,
  `NTHREAD=128` trips `cudaErrorIllegalAddress` on Turing.
* **`Opening GPUIRR lib. SSE ver.` under the gpu variant** is expected:
  upstream NBODY7 pairs `gpunb=GPU` with `gpuirr=SSE`, keeping
  irregular forces on the CPU SSE path even in GPU builds.
* **AMUSE's `AMUSE_CUDA` macro** detects nvcc and sets `CUDA_TK` /
  `NVCC`, but on some hosts leaves `CUDA_FLAGS` / `CUDA_LDFLAGS` empty.
  The gpu worker rule adds explicit `-L $(CUDA_TK)/lib
  -L $(CUDA_TK)/lib64 -lcudart` on top of `$(CUDA_LDFLAGS)` so the link
  succeeds regardless.
