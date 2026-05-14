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
  numerical-equivalence check on new hardware.
