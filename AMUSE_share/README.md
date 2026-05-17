# AMUSE_share/

Distribution-ready bundle for the amuse-nbody7 AMUSE community-code
interface. Devs unpack `amusenb7v1.tar.gz` into AMUSE's `src/` (or
`packages/`); the full package documentation, quick-start, and demo
target reference lives at:

* [amuse_nbody7/README.md](amuse_nbody7/README.md)

This top-level directory also holds:

* `amusenb7v1.tar.gz` — the shipped tarball (contents of `amuse_nbody7/`)
* `tests_extra/` — standalone smoke + equivalence test scripts kept
  outside the tarball; useful for reproducing the cross-variant
  numerical-equivalence check on new hardware. The python notebook
  is for demo to the audience.

Download amuse interface:

amuse-nbody7 is not yet an official package that is shipped
with the standard AMUSE bundle. Once you have successfully
installed AMUSE on your PC/Server/Box (install at least amuse-framework;
check sanity by installing amuse_hermite and, if you have
a CUDA platform, also by installing sapporo_light):

cd <amuse_checkout_dir>/src
(e.g. cd ~/amuse-2026.3.0/src)
Here all other shipped AMUSE interfaces stay. 

Download NBODY7's AMUSE interface manually from github public repository:
(user: sambaranb, package: updated-nbody7, branch: mpi-parallel)
curl -LO https://raw.githubusercontent.com/sambaranb/updated-nbody7/mpi-parallel/AMUSE_share/amusenb7v1.tar.gz
or
wget https://raw.githubusercontent.com/sambaranb/updated-nbody7/mpi-parallel/AMUSE_share/amusenb7v1.tar.gz

Unpack the tarball:
tar zxvf amusenb7v1.tar.gz

cd amuse_nbody7

Activate your AMUSE conda environment (say, Amuse-env) if it isn't activated already:
conda activate Amuse-env

Follow the instructions in README.md (linked above) for installing the amuse_nbody7 workers.
