* ============================================================================
*
*                    N B O D Y 7   for   A M U S E
*                    *****************************
*
*   Fortran-side bridge between the AMUSE worker (nbody7_worker) and the
*   refactored NBODY7 community code on the mpi-parallel branch.
*
*   Modelled on amuse_nbody6xx/interface.f90 (Maxwell X. TSAI/CAI, 2013-2014)
*   updated to use the /AMUSEBLK/ COMMON block introduced in Step 2
*   (Block/amuse.h) instead of modifying common6.h directly, and extended
*   to cover the full NBODY7 parameter set.
*
*   File format
*   -----------
*   Despite the .f90 extension, this file is FIXED-FORM Fortran (column 7
*   for code, * in column 1 for comments, & in column 6 of continuation
*   lines). The reason: common6.h and amuse.h are both fixed-form F77
*   includes that declare the NBODY7 COMMON blocks; including them from a
*   free-form file mis-parses the column-1 comment lines and silently
*   drops the COMMON declarations. The Makefile compiles this file with
*   -ffixed-form -ffixed-line-length-none.
*
*   The AMUSE Python framework reaches every public function in this
*   module over an MPI-IPC channel. Each function returns 0 on success
*   or a non-zero error code; results travel back through the
*   amusifier-generated nbody7_worker.f90 which copies the INTENT(out)
*   arguments into the MPI reply.
*
*   Initialization model
*   --------------------
*   AMUSE state machine:  UNINITIALISED -> INITIALISED -> EVOLVED ...
*     initialize_code()   -> install AMUSE-friendly defaults, set amusein=1
*     set_*()             -> override defaults via setters
*     new_particle() ...  -> stage particles into BODY/X/XDOT/NAME
*     commit_parameters() -> finalize parameter set (just N = ...)
*     commit_particles()  -> N has its final value
*     evolve_model(t_end) -> first call invokes CALL NBODY6 (which runs
*                            START / FPOLY0 / FPOLY2 / etc. with READs
*                            short-circuited via amusein=1); subsequent
*                            calls drive the integrator via INTAMUSE.
*
*   v0 scope
*   --------
*   CPU-only NBODY7 behind the MPI-IPC AMUSE worker. CUDA/Metal variants
*   share this same interface - only the worker build target differs.
*
* ============================================================================

      MODULE Nbody7Interface
*
*     No IMPLICIT NONE: common6.h declares IMPLICIT REAL*8 (A-H,O-Z)
*     and relies on the Fortran default I-N -> INTEGER rule for the
*     COMMON-block counters (NTOT, NZERO, NSUB, ...). IMPLICIT NONE
*     here would suppress that default and force every NBODY7 variable
*     to be redeclared explicitly.
*

*     -------------------------------------------------------------------
*     Module-scope state.
*
*     last_index      - next free particle slot when staging via
*                       new_particle.
*     nb6_initialised - 0 until the first evolve_model() invokes
*                       CALL NBODY6 (which is what turns staged
*                       particles into a running integrator).
*                       Subsequent evolve_model() calls just drive
*                       INTAMUSE.
*     -------------------------------------------------------------------
      INTEGER :: last_index = 1
      INTEGER :: nb6_initialised = 0

      CONTAINS

* ===========================================================================
* Lifecycle
* ===========================================================================

      FUNCTION initialize_code() RESULT(ierr)
      INCLUDE 'common6.h'
      INCLUDE 'amuse.h'
      INTEGER :: ierr

*     Switch the NBODY7 source tree into AMUSE mode: every gated READ
*     in input.f / data.f / scale.f is now short-circuited.
      amusein = 1

*     AMUSE-friendly defaults. These can be overridden via the set_*
*     functions before commit_parameters() / first evolve_model().
      KSTART_AMUSE = 1
      TCOMP_AMUSE  = 1.0d4
      NRAND_AMUSE  = 4353
      RS0_AMUSE    = 0.14d0
      TCRITp       = 1.0d6
      isernb       = 40
      iserreg      = 40

*     ARCHAIN parameters: speed of light in N-body units, BH flag,
*     and disruption-mode flag. Read from STDIN by ksint.f / chain.f
*     in the standalone path; injected here under AMUSE so that the
*     first KSINT or CHAIN invocation does not fault on a missing
*     fort.5 read. Override via set_archain_params() before commit.
      CLIGHT_AMUSE = 2.0d4
      NBH_AMUSE    = 0
      IDIS_AMUSE   = 0

*     Particle-setup parameters that are read from STDIN on the same
*     DATA line as ALPHAS/BODY1/BODYN. Re-injected into the matching
*     COMMON variables inside Block/data.f after ZERO has cleared them.
      NBIN0_AMUSE  = 0
      NHI0_AMUSE   = 0
      ZMET_AMUSE   = 0.02d0
      DTPLOT_AMUSE = 10.0d0
      EPOCH0_AMUSE = 0.0d0

*     SSE / BSE parameters that data.f normally reads from the file
*     'input_bse'. Defaults match the canonical input_bse shipped
*     with NBODY7 production runs. Active only when KZ(19) >= 3
*     (Tout-Hurley stellar evolution); harmless otherwise.
*     Line 1 of input_bse: neta bwind hewind wconst alpha1 lambda
      neta   = 0.5d0
      bwind  = 0.0d0
      hewind = 1.0d0
      wconst = 1.0d0
      ALPHA1 = 1.0d0
      LAMBDA = 0.5d0
*     Line 2: ceflag tflag ifflag wdflag bhflag nsflag mxns
      ceflag = 0
      tflag  = 0
      ifflag = 0
      wdflag = 1
      bhflag = 3
      nsflag = 3
      mxns   = 2.5d0
*     Line 3: psflag kmech ecflag edflag
      psflag = 1
      kmech  = 1
      ecflag = 1
      edflag = 0
*     Line 4: disp beta psii acc2 epsnov ftzacc fmrg eddfac gamm1
      disp   = 265.0d0
      beta   = 0.125d0
      psii   = 1.0d0
      acc2   = 1.5d0
      epsnov = 0.001d0
      ftzacc = 0.50d0
      fmrg   = 0.50d0
      eddfac = 1.0d0
      gamm1  = -1.0d0

*     Mirror standard NBODY7 input parameters that INPUT would normally
*     parse from STDIN. Users can override via setters.
      NFIX   = 1
      NCRIT  = 10
      NRUN   = 1
      NNBMAX = 300

      ETAI   = 0.05d0
      ETAR   = 0.05d0
      DTADJ  = 1.0d0
      DELTAT = 1.0d0
      TCRIT  = 1.0d10
      QE     = 2.0d-5
      RBAR   = 1.0d0
      ZMBAR  = 0.7d0

      DTMIN  = 1.0d-4
      RMIN   = 0.01d0
      ETAU   = 0.1d0
      ECLOSE = 1.0d0
      GMIN   = 1.0d-6
      GMAX   = 0.01d0

*     KZ defaults for the AMUSE-driven flow. Particles arrive via
*     new_particle (NOT fort.10) and stellar evolution is opt-in.
*     Users can rewrite any KZ slot via set_kz() before commit.
*       KZ(1)  = 1   COMMON-save dump on output (cheap, harmless)
*       KZ(3)  = 1   write energy diagnostics
*       KZ(5)  = 1   standard cluster initial conditions
*       KZ(11) = -1  ARCHAIN post-Newtonian chain (the v0 highlight)
*       KZ(19) = 0   no Tout-Hurley stellar evolution by default
*                    (opt in via set_kz(19, 3) once input_bse is wired)
*       KZ(22) = 0   particles supplied via the AMUSE new_particle
*                    interface; the DATA fort.10 / IMF / SETUP paths
*                    are bypassed by the amusein=1 gate in Block/data.f
      KZ(1:50) = 0
      KZ(1)  = 1
      KZ(3)  = 1
      KZ(5)  = 1
      KZ(11) = -1
      KZ(19) = 0
      KZ(22) = 0

*     Particle-staging cursor.
      last_index      = 1
      nb6_initialised = 0

      ierr = 0
      END FUNCTION initialize_code


      FUNCTION cleanup_code() RESULT(ierr)
      INTEGER :: ierr
      ierr = 0
      END FUNCTION cleanup_code


      FUNCTION commit_parameters() RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER :: ierr
      N = last_index - 1
      NTOT  = N
      NZERO = N
      ierr = 0
      END FUNCTION commit_parameters


      FUNCTION recommit_parameters() RESULT(ierr)
      INTEGER :: ierr
      ierr = 0
      END FUNCTION recommit_parameters


* ===========================================================================
* Particle management
* ===========================================================================

      FUNCTION new_particle(index_of_the_particle, mass,
     &                      x_in, y_in, z_in,
     &                      vx_in, vy_in, vz_in,
     &                      radius_in) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(OUT) :: index_of_the_particle
      DOUBLE PRECISION, INTENT(IN) :: mass
      DOUBLE PRECISION, INTENT(IN) :: x_in, y_in, z_in
      DOUBLE PRECISION, INTENT(IN) :: vx_in, vy_in, vz_in
      DOUBLE PRECISION, INTENT(IN) :: radius_in
      INTEGER :: ierr, i

      i = last_index
      BODY(i)    = mass
      X(1, i)    = x_in
      X(2, i)    = y_in
      X(3, i)    = z_in
      XDOT(1, i) = vx_in
      XDOT(2, i) = vy_in
      XDOT(3, i) = vz_in
      RADIUS(i)  = radius_in
      NAME(i)    = i

      index_of_the_particle = i
      last_index = i + 1
      ierr = 0
      END FUNCTION new_particle


      FUNCTION delete_particle(index_of_the_particle) RESULT(ierr)
      INTEGER, INTENT(IN) :: index_of_the_particle
      INTEGER :: ierr
*     Particle removal mid-run is non-trivial in NBODY7 (touches KS
*     pairs, chain, neighbour lists). Not supported in v0.
      ierr = -2
      END FUNCTION delete_particle


      FUNCTION commit_particles() RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER :: ierr
      N     = last_index - 1
      NTOT  = N
      NZERO = N
      ierr = 0
      END FUNCTION commit_particles


      FUNCTION recommit_particles() RESULT(ierr)
      INTEGER :: ierr
      ierr = 0
      END FUNCTION recommit_particles


      FUNCTION get_number_of_particles(num) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(OUT) :: num
      INTEGER :: ierr
      IF (N .LE. 0) THEN
          num = last_index - 1
      ELSE
          num = N
      END IF
      ierr = 0
      END FUNCTION get_number_of_particles


      FUNCTION get_index_of_first_particle(idx) RESULT(ierr)
      INTEGER, INTENT(OUT) :: idx
      INTEGER :: ierr
      idx = 1
      ierr = 0
      END FUNCTION get_index_of_first_particle


      FUNCTION get_index_of_next_particle(idx_in, idx_next)
     &                                    RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN)  :: idx_in
      INTEGER, INTENT(OUT) :: idx_next
      INTEGER :: ierr
      IF (idx_in .LT. N) THEN
          idx_next = idx_in + 1
          ierr = 0
      ELSE
          idx_next = 0
          ierr = 1
      END IF
      END FUNCTION get_index_of_next_particle


* ===========================================================================
* State accessors
* ===========================================================================

      FUNCTION get_state(idx, mass, x_o, y_o, z_o,
     &                   vx_o, vy_o, vz_o, r_o) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN)  :: idx
      DOUBLE PRECISION, INTENT(OUT) :: mass, x_o, y_o, z_o
      DOUBLE PRECISION, INTENT(OUT) :: vx_o, vy_o, vz_o, r_o
      INTEGER :: ierr
      mass = BODY(idx)
      x_o  = X(1, idx)
      y_o  = X(2, idx)
      z_o  = X(3, idx)
      vx_o = XDOT(1, idx)
      vy_o = XDOT(2, idx)
      vz_o = XDOT(3, idx)
      r_o  = RADIUS(idx)
      ierr = 0
      END FUNCTION get_state


      FUNCTION set_state(idx, mass, x_in, y_in, z_in,
     &                   vx_in, vy_in, vz_in, r_in) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN) :: idx
      DOUBLE PRECISION, INTENT(IN) :: mass, x_in, y_in, z_in
      DOUBLE PRECISION, INTENT(IN) :: vx_in, vy_in, vz_in, r_in
      INTEGER :: ierr
      BODY(idx)    = mass
      X(1, idx)    = x_in
      X(2, idx)    = y_in
      X(3, idx)    = z_in
      XDOT(1, idx) = vx_in
      XDOT(2, idx) = vy_in
      XDOT(3, idx) = vz_in
      RADIUS(idx)  = r_in
      ierr = 0
      END FUNCTION set_state


      FUNCTION get_mass(idx, mass) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN) :: idx
      DOUBLE PRECISION, INTENT(OUT) :: mass
      INTEGER :: ierr
      mass = BODY(idx)
      ierr = 0
      END FUNCTION get_mass


      FUNCTION set_mass(idx, mass) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN) :: idx
      DOUBLE PRECISION, INTENT(IN) :: mass
      INTEGER :: ierr
      BODY(idx) = mass
      ierr = 0
      END FUNCTION set_mass


      FUNCTION get_position(idx, x_o, y_o, z_o) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN) :: idx
      DOUBLE PRECISION, INTENT(OUT) :: x_o, y_o, z_o
      INTEGER :: ierr
      x_o = X(1, idx)
      y_o = X(2, idx)
      z_o = X(3, idx)
      ierr = 0
      END FUNCTION get_position


      FUNCTION set_position(idx, x_in, y_in, z_in) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN) :: idx
      DOUBLE PRECISION, INTENT(IN) :: x_in, y_in, z_in
      INTEGER :: ierr
      X(1, idx) = x_in
      X(2, idx) = y_in
      X(3, idx) = z_in
      ierr = 0
      END FUNCTION set_position


      FUNCTION get_velocity(idx, vx_o, vy_o, vz_o) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN) :: idx
      DOUBLE PRECISION, INTENT(OUT) :: vx_o, vy_o, vz_o
      INTEGER :: ierr
      vx_o = XDOT(1, idx)
      vy_o = XDOT(2, idx)
      vz_o = XDOT(3, idx)
      ierr = 0
      END FUNCTION get_velocity


      FUNCTION set_velocity(idx, vx_in, vy_in, vz_in) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN) :: idx
      DOUBLE PRECISION, INTENT(IN) :: vx_in, vy_in, vz_in
      INTEGER :: ierr
      XDOT(1, idx) = vx_in
      XDOT(2, idx) = vy_in
      XDOT(3, idx) = vz_in
      ierr = 0
      END FUNCTION set_velocity


      FUNCTION get_name(idx, val) RESULT(ierr)
*     NBODY7 keeps a per-particle integer NAME (in /PAIRS/ COMMON via
*     common6.h). It is preserved across KS / chain regularizations
*     and used to track the identity of a star through binary
*     formation, mergers, and ejections.
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN) :: idx
      INTEGER, INTENT(OUT) :: val
      INTEGER :: ierr
      val = NAME(idx)
      ierr = 0
      END FUNCTION get_name


      FUNCTION set_name(idx, val) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN) :: idx, val
      INTEGER :: ierr
      NAME(idx) = val
      ierr = 0
      END FUNCTION set_name


      FUNCTION get_kstar(idx, val) RESULT(ierr)
*     NBODY7 stellar-evolution type code (in /NBODY/ COMMON via
*     common6.h). The standard SSE/BSE encoding is used:
*       0  low-mass main sequence    7  helium main sequence
*       1  main sequence             8  helium giant
*       2  Hertzsprung gap           9  helium HG/AGB
*       3  red giant branch         10  helium WD
*       4  core helium burning      11  CO white dwarf
*       5  AGB                      12  ONe white dwarf
*       6  thermal-pulsing AGB      13  neutron star
*                                   14  black hole
*                                   15  massless remnant
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN) :: idx
      INTEGER, INTENT(OUT) :: val
      INTEGER :: ierr
      val = KSTAR(idx)
      ierr = 0
      END FUNCTION get_kstar


      FUNCTION set_kstar(idx, val) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN) :: idx, val
      INTEGER :: ierr
      KSTAR(idx) = val
      ierr = 0
      END FUNCTION set_kstar


      FUNCTION get_radius(idx, r_o) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN) :: idx
      DOUBLE PRECISION, INTENT(OUT) :: r_o
      INTEGER :: ierr
      r_o = RADIUS(idx)
      ierr = 0
      END FUNCTION get_radius


      FUNCTION set_radius(idx, r_in) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN) :: idx
      DOUBLE PRECISION, INTENT(IN) :: r_in
      INTEGER :: ierr
      RADIUS(idx) = r_in
      ierr = 0
      END FUNCTION set_radius


* ===========================================================================
* Diagnostics
* ===========================================================================

      FUNCTION get_acceleration(idx, ax_o, ay_o, az_o) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN) :: idx
      DOUBLE PRECISION, INTENT(OUT) :: ax_o, ay_o, az_o
      INTEGER :: ierr
      ax_o = FI(1, idx) + FR(1, idx)
      ay_o = FI(2, idx) + FR(2, idx)
      az_o = FI(3, idx) + FR(3, idx)
      ierr = 0
      END FUNCTION get_acceleration


      FUNCTION get_potential(idx, phi) RESULT(ierr)
*     The NBODY7 common6.h does not expose a per-particle potential
*     array (only the system-total POT). Compute it on demand via
*     direct N-body summation. This is O(N) per call - fine for the
*     diagnostic role this method plays in AMUSE.
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN) :: idx
      DOUBLE PRECISION, INTENT(OUT) :: phi
      INTEGER :: ierr, j
      DOUBLE PRECISION :: dx, dy, dz, r2
      phi = 0.0d0
      DO j = 1, N
          IF (j .NE. idx) THEN
              dx = X(1, j) - X(1, idx)
              dy = X(2, j) - X(2, idx)
              dz = X(3, j) - X(3, idx)
              r2 = dx*dx + dy*dy + dz*dz
              phi = phi - BODY(j) / SQRT(r2)
          END IF
      END DO
      ierr = 0
      END FUNCTION get_potential


      FUNCTION get_kinetic_energy(ke) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: ke
      INTEGER :: ierr
      ke = ZKIN
      ierr = 0
      END FUNCTION get_kinetic_energy


      FUNCTION get_potential_energy(pe) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: pe
      INTEGER :: ierr
*     NBODY7 stores POT as +|U|; AMUSE expects U <= 0.
      pe = -POT
      ierr = 0
      END FUNCTION get_potential_energy


      FUNCTION get_total_mass(mtot) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: mtot
      INTEGER :: ierr
      mtot = ZMASS
      ierr = 0
      END FUNCTION get_total_mass


      FUNCTION get_total_radius(rtot) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: rtot
      INTEGER :: ierr
      rtot = RBAR
      ierr = 0
      END FUNCTION get_total_radius


      FUNCTION get_center_of_mass_position(cmx, cmy, cmz) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: cmx, cmy, cmz
      INTEGER :: ierr, idx
      DOUBLE PRECISION :: mtot
      cmx = 0.0d0
      cmy = 0.0d0
      cmz = 0.0d0
      mtot = 0.0d0
      DO idx = 1, N
          mtot = mtot + BODY(idx)
          cmx  = cmx  + BODY(idx) * X(1, idx)
          cmy  = cmy  + BODY(idx) * X(2, idx)
          cmz  = cmz  + BODY(idx) * X(3, idx)
      END DO
      IF (mtot .GT. 0.0d0) THEN
          cmx = cmx / mtot
          cmy = cmy / mtot
          cmz = cmz / mtot
      END IF
      ierr = 0
      END FUNCTION get_center_of_mass_position


      FUNCTION get_center_of_mass_velocity(vx_o, vy_o, vz_o)
     &                                     RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: vx_o, vy_o, vz_o
      INTEGER :: ierr, idx
      DOUBLE PRECISION :: mtot
      vx_o = 0.0d0
      vy_o = 0.0d0
      vz_o = 0.0d0
      mtot = 0.0d0
      DO idx = 1, N
          mtot = mtot + BODY(idx)
          vx_o = vx_o + BODY(idx) * XDOT(1, idx)
          vy_o = vy_o + BODY(idx) * XDOT(2, idx)
          vz_o = vz_o + BODY(idx) * XDOT(3, idx)
      END DO
      IF (mtot .GT. 0.0d0) THEN
          vx_o = vx_o / mtot
          vy_o = vy_o / mtot
          vz_o = vz_o / mtot
      END IF
      ierr = 0
      END FUNCTION get_center_of_mass_velocity


* ===========================================================================
* Time / evolution
* ===========================================================================

      FUNCTION get_time(t_o) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: t_o
      INTEGER :: ierr
      t_o = TTOT
      ierr = 0
      END FUNCTION get_time


      FUNCTION get_time_step(dt_o) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: dt_o
      INTEGER :: ierr
      dt_o = ETAI
      ierr = 0
      END FUNCTION get_time_step


      FUNCTION get_begin_time(t_o) RESULT(ierr)
      DOUBLE PRECISION, INTENT(OUT) :: t_o
      INTEGER :: ierr
      t_o = 0.0d0
      ierr = 0
      END FUNCTION get_begin_time


      FUNCTION set_begin_time(t_in) RESULT(ierr)
      DOUBLE PRECISION, INTENT(IN) :: t_in
      INTEGER :: ierr
*     AMUSE allows shifting the simulation epoch; NBODY7 starts at
*     TTOT=0. v0: accept and ignore (placeholder for a future
*     TIME-offset wiring).
      ierr = 0
      END FUNCTION set_begin_time


      FUNCTION evolve_model(t_end) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: t_end
      INTEGER :: ierr

*     First call: do the heavy initialisation (CALL NBODY6 runs START,
*     which calls INPUT/DATA/SCALE/FPOLY0/FPOLY2 with READs short-
*     circuited because amusein=1). Subsequent calls just step the
*     integrator forward via INTAMUSE.
      IF (nb6_initialised .EQ. 0) THEN
          CALL NBODY6
          nb6_initialised = 1
      END IF

      DO WHILE (TTOT .LT. t_end)
          CALL INTAMUSE
      END DO

      ierr = 0
      END FUNCTION evolve_model


      FUNCTION synchronize_model() RESULT(ierr)
      INTEGER :: ierr
      ierr = 0
      END FUNCTION synchronize_model


* ===========================================================================
* Softening (NBODY7 has no Plummer softening; eps^2 = 0 always)
* ===========================================================================

      FUNCTION get_eps2(eps2_o) RESULT(ierr)
      DOUBLE PRECISION, INTENT(OUT) :: eps2_o
      INTEGER :: ierr
      eps2_o = 0.0d0
      ierr = 0
      END FUNCTION get_eps2


      FUNCTION set_eps2(eps2_in) RESULT(ierr)
      DOUBLE PRECISION, INTENT(IN) :: eps2_in
      INTEGER :: ierr
*     NBODY7 is a regularized integrator; explicit Plummer softening
*     is not supported. Accept and ignore.
      ierr = 0
      END FUNCTION set_eps2


* ===========================================================================
* NBODY7-specific parameters
* ===========================================================================

      FUNCTION get_kz(option, val) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN)  :: option
      INTEGER, INTENT(OUT) :: val
      INTEGER :: ierr
      IF (option .GE. 1 .AND. option .LE. 50) THEN
          val = KZ(option)
          ierr = 0
      ELSE
          val  = 0
          ierr = -1
      END IF
      END FUNCTION get_kz


      FUNCTION set_kz(option, val) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN) :: option, val
      INTEGER :: ierr
      IF (option .GE. 1 .AND. option .LE. 50) THEN
          KZ(option) = val
          ierr = 0
      ELSE
          ierr = -1
      END IF
      END FUNCTION set_kz


      FUNCTION get_nnbmax(val) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(OUT) :: val
      INTEGER :: ierr
      val = NNBMAX
      ierr = 0
      END FUNCTION get_nnbmax


      FUNCTION set_nnbmax(val) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN) :: val
      INTEGER :: ierr
      NNBMAX = val
      ierr = 0
      END FUNCTION set_nnbmax


      FUNCTION get_eta(eta_o) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: eta_o
      INTEGER :: ierr
      eta_o = ETAI
      ierr = 0
      END FUNCTION get_eta


      FUNCTION set_eta(eta_in) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: eta_in
      INTEGER :: ierr
      ETAI = eta_in
      ETAR = eta_in
      ierr = 0
      END FUNCTION set_eta


      FUNCTION get_etai(eta_o) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: eta_o
      INTEGER :: ierr
      eta_o = ETAI
      ierr = 0
      END FUNCTION get_etai


      FUNCTION set_etai(eta_in) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: eta_in
      INTEGER :: ierr
      ETAI = eta_in
      ierr = 0
      END FUNCTION set_etai


      FUNCTION get_etar(eta_o) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: eta_o
      INTEGER :: ierr
      eta_o = ETAR
      ierr = 0
      END FUNCTION get_etar


      FUNCTION set_etar(eta_in) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: eta_in
      INTEGER :: ierr
      ETAR = eta_in
      ierr = 0
      END FUNCTION set_etar


      FUNCTION get_etau(eta_o) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: eta_o
      INTEGER :: ierr
      eta_o = ETAU
      ierr = 0
      END FUNCTION get_etau


      FUNCTION set_etau(eta_in) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: eta_in
      INTEGER :: ierr
      ETAU = eta_in
      ierr = 0
      END FUNCTION set_etau


      FUNCTION get_rbar(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = RBAR
      ierr = 0
      END FUNCTION get_rbar


      FUNCTION set_rbar(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      RBAR = val
      ierr = 0
      END FUNCTION set_rbar


      FUNCTION get_zmbar(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = ZMBAR
      ierr = 0
      END FUNCTION get_zmbar


      FUNCTION set_zmbar(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      ZMBAR = val
      ierr = 0
      END FUNCTION set_zmbar


      FUNCTION get_qe(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = QE
      ierr = 0
      END FUNCTION get_qe


      FUNCTION set_qe(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      QE = val
      ierr = 0
      END FUNCTION set_qe


      FUNCTION get_dtadj(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = DTADJ
      ierr = 0
      END FUNCTION get_dtadj


      FUNCTION set_dtadj(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      DTADJ = val
      ierr = 0
      END FUNCTION set_dtadj


      FUNCTION get_deltat(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = DELTAT
      ierr = 0
      END FUNCTION get_deltat


      FUNCTION set_deltat(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      DELTAT = val
      ierr = 0
      END FUNCTION set_deltat


      FUNCTION get_dtmin(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = DTMIN
      ierr = 0
      END FUNCTION get_dtmin


      FUNCTION set_dtmin(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      DTMIN = val
      ierr = 0
      END FUNCTION set_dtmin


      FUNCTION get_rmin(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = RMIN
      ierr = 0
      END FUNCTION get_rmin


      FUNCTION set_rmin(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      RMIN = val
      ierr = 0
      END FUNCTION set_rmin


      FUNCTION get_eclose(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = ECLOSE
      ierr = 0
      END FUNCTION get_eclose


      FUNCTION set_eclose(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      ECLOSE = val
      ierr = 0
      END FUNCTION set_eclose


      FUNCTION get_gmin(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = GMIN
      ierr = 0
      END FUNCTION get_gmin


      FUNCTION set_gmin(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      GMIN = val
      ierr = 0
      END FUNCTION set_gmin


      FUNCTION get_gmax(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = GMAX
      ierr = 0
      END FUNCTION get_gmax


      FUNCTION set_gmax(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      GMAX = val
      ierr = 0
      END FUNCTION set_gmax


      FUNCTION get_tcrit(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = TCRIT
      ierr = 0
      END FUNCTION get_tcrit


      FUNCTION set_tcrit(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      TCRIT = val
      ierr = 0
      END FUNCTION set_tcrit


      FUNCTION get_tcomp(val) RESULT(ierr)
      INCLUDE 'amuse.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = TCOMP_AMUSE
      ierr = 0
      END FUNCTION get_tcomp


      FUNCTION set_tcomp(val) RESULT(ierr)
      INCLUDE 'amuse.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      TCOMP_AMUSE = val
      ierr = 0
      END FUNCTION set_tcomp


      FUNCTION get_kstart(val) RESULT(ierr)
      INCLUDE 'amuse.h'
      INTEGER, INTENT(OUT) :: val
      INTEGER :: ierr
      val = KSTART_AMUSE
      ierr = 0
      END FUNCTION get_kstart


      FUNCTION set_kstart(val) RESULT(ierr)
      INCLUDE 'amuse.h'
      INTEGER, INTENT(IN) :: val
      INTEGER :: ierr
      KSTART_AMUSE = val
      ierr = 0
      END FUNCTION set_kstart


      FUNCTION get_nrand(val) RESULT(ierr)
      INCLUDE 'amuse.h'
      INTEGER, INTENT(OUT) :: val
      INTEGER :: ierr
      val = NRAND_AMUSE
      ierr = 0
      END FUNCTION get_nrand


      FUNCTION set_nrand(val) RESULT(ierr)
      INCLUDE 'amuse.h'
      INTEGER, INTENT(IN) :: val
      INTEGER :: ierr
      NRAND_AMUSE = val
      ierr = 0
      END FUNCTION set_nrand


      FUNCTION set_archain_params(clight_in, nbh_in, idis_in)
     &                            RESULT(ierr)
*     Set the three ARCHAIN startup parameters (speed of light in
*     N-body units, BH flag, disruption-mode flag) that ksint.f and
*     chain.f read from STDIN in the standalone path. Stored in
*     /AMUSEBLK/ and copied into the relevant COMMON blocks the first
*     time KSINT or CHAIN is invoked.
      INCLUDE 'amuse.h'
      DOUBLE PRECISION, INTENT(IN) :: clight_in
      INTEGER, INTENT(IN) :: nbh_in, idis_in
      INTEGER :: ierr
      CLIGHT_AMUSE = clight_in
      NBH_AMUSE    = nbh_in
      IDIS_AMUSE   = idis_in
      ierr = 0
      END FUNCTION set_archain_params


      FUNCTION get_archain_params(clight_o, nbh_o, idis_o) RESULT(ierr)
*     Read the current ARCHAIN startup parameters back out of
*     /AMUSEBLK/. See set_archain_params for the semantics.
      INCLUDE 'amuse.h'
      DOUBLE PRECISION, INTENT(OUT) :: clight_o
      INTEGER, INTENT(OUT) :: nbh_o, idis_o
      INTEGER :: ierr
      clight_o = CLIGHT_AMUSE
      nbh_o    = NBH_AMUSE
      idis_o   = IDIS_AMUSE
      ierr = 0
      END FUNCTION get_archain_params


      FUNCTION get_tstar(val) RESULT(ierr)
*     Conversion factor from N-body time to physical Myr, populated
*     by Block/units.f as part of the START sequence. Returns 0
*     until the first evolve_model() invocation triggers CALL NBODY6,
*     which runs SCALE / UNITS to compute it from RBAR, ZMASS, ZMBAR.
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = TSTAR
      ierr = 0
      END FUNCTION get_tstar


      FUNCTION get_vstar(val) RESULT(ierr)
*     Conversion factor from N-body velocity to km/s. Populated by
*     Block/units.f during the START sequence (zero until then).
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = VSTAR
      ierr = 0
      END FUNCTION get_vstar


      FUNCTION get_rs0(val) RESULT(ierr)
      INCLUDE 'amuse.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = RS0_AMUSE
      ierr = 0
      END FUNCTION get_rs0


      FUNCTION set_rs0(val) RESULT(ierr)
      INCLUDE 'amuse.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      RS0_AMUSE = val
      ierr = 0
      END FUNCTION set_rs0


* ===========================================================================
* Particle-setup parameters (one DATA line of input_bse companion vars)
*
* These five live in /AMUSEBLK/ and are re-injected into the matching
* COMMON variables inside Block/data.f (because Block/zero.f clears
* NBIN0, NHI0 and EPOCH0 at every START). Setters write to the
* /AMUSEBLK/ slot; getters read it back, so what-you-set is what-you-
* get even before the first evolve_model triggers START.
* ===========================================================================

      FUNCTION get_nbin0(val) RESULT(ierr)
      INCLUDE 'amuse.h'
      INTEGER, INTENT(OUT) :: val
      INTEGER :: ierr
      val = NBIN0_AMUSE
      ierr = 0
      END FUNCTION get_nbin0


      FUNCTION set_nbin0(val) RESULT(ierr)
      INCLUDE 'amuse.h'
      INTEGER, INTENT(IN) :: val
      INTEGER :: ierr
      NBIN0_AMUSE = val
      ierr = 0
      END FUNCTION set_nbin0


      FUNCTION get_nhi0(val) RESULT(ierr)
      INCLUDE 'amuse.h'
      INTEGER, INTENT(OUT) :: val
      INTEGER :: ierr
      val = NHI0_AMUSE
      ierr = 0
      END FUNCTION get_nhi0


      FUNCTION set_nhi0(val) RESULT(ierr)
      INCLUDE 'amuse.h'
      INTEGER, INTENT(IN) :: val
      INTEGER :: ierr
      NHI0_AMUSE = val
      ierr = 0
      END FUNCTION set_nhi0


      FUNCTION get_zmet(val) RESULT(ierr)
      INCLUDE 'amuse.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = ZMET_AMUSE
      ierr = 0
      END FUNCTION get_zmet


      FUNCTION set_zmet(val) RESULT(ierr)
      INCLUDE 'amuse.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      ZMET_AMUSE = val
      ierr = 0
      END FUNCTION set_zmet


      FUNCTION get_dtplot(val) RESULT(ierr)
      INCLUDE 'amuse.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = DTPLOT_AMUSE
      ierr = 0
      END FUNCTION get_dtplot


      FUNCTION set_dtplot(val) RESULT(ierr)
      INCLUDE 'amuse.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      DTPLOT_AMUSE = val
      ierr = 0
      END FUNCTION set_dtplot


      FUNCTION get_epoch0(val) RESULT(ierr)
      INCLUDE 'amuse.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = EPOCH0_AMUSE
      ierr = 0
      END FUNCTION get_epoch0


      FUNCTION set_epoch0(val) RESULT(ierr)
      INCLUDE 'amuse.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      EPOCH0_AMUSE = val
      ierr = 0
      END FUNCTION set_epoch0


* ===========================================================================
* SSE / BSE stellar-evolution parameters (active when KZ(19) >= 3)
*
* All 26 of these are read by Block/data.f from the file 'input_bse'
* in the standalone path. Under AMUSE, the gated read is skipped and
* the COMMON-block values (in /VALUE1-7/, /FLAGS/, /FLAGS2/, /FLAGS3/
* via common6.h) are used as-is. Block/zero.f does not clear them, so
* values written here through commit_parameters survive into the
* integration. Defaults installed by initialize_code match the
* canonical input_bse shipped with NBODY7 production runs.
*
* Group A: 16 REAL*8 parameters
* ===========================================================================

      FUNCTION get_bse_neta(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = neta
      ierr = 0
      END FUNCTION get_bse_neta


      FUNCTION set_bse_neta(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      neta = val
      ierr = 0
      END FUNCTION set_bse_neta


      FUNCTION get_bse_bwind(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = bwind
      ierr = 0
      END FUNCTION get_bse_bwind


      FUNCTION set_bse_bwind(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      bwind = val
      ierr = 0
      END FUNCTION set_bse_bwind


      FUNCTION get_bse_hewind(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = hewind
      ierr = 0
      END FUNCTION get_bse_hewind


      FUNCTION set_bse_hewind(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      hewind = val
      ierr = 0
      END FUNCTION set_bse_hewind


      FUNCTION get_bse_wconst(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = wconst
      ierr = 0
      END FUNCTION get_bse_wconst


      FUNCTION set_bse_wconst(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      wconst = val
      ierr = 0
      END FUNCTION set_bse_wconst


      FUNCTION get_bse_alpha1(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = ALPHA1
      ierr = 0
      END FUNCTION get_bse_alpha1


      FUNCTION set_bse_alpha1(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      ALPHA1 = val
      ierr = 0
      END FUNCTION set_bse_alpha1


      FUNCTION get_bse_lambda(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = LAMBDA
      ierr = 0
      END FUNCTION get_bse_lambda


      FUNCTION set_bse_lambda(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      LAMBDA = val
      ierr = 0
      END FUNCTION set_bse_lambda


      FUNCTION get_bse_mxns(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = mxns
      ierr = 0
      END FUNCTION get_bse_mxns


      FUNCTION set_bse_mxns(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      mxns = val
      ierr = 0
      END FUNCTION set_bse_mxns


      FUNCTION get_bse_ftzacc(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = ftzacc
      ierr = 0
      END FUNCTION get_bse_ftzacc


      FUNCTION set_bse_ftzacc(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      ftzacc = val
      ierr = 0
      END FUNCTION set_bse_ftzacc


      FUNCTION get_bse_disp(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = disp
      ierr = 0
      END FUNCTION get_bse_disp


      FUNCTION set_bse_disp(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      disp = val
      ierr = 0
      END FUNCTION set_bse_disp


      FUNCTION get_bse_beta(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = beta
      ierr = 0
      END FUNCTION get_bse_beta


      FUNCTION set_bse_beta(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      beta = val
      ierr = 0
      END FUNCTION set_bse_beta


      FUNCTION get_bse_psii(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = psii
      ierr = 0
      END FUNCTION get_bse_psii


      FUNCTION set_bse_psii(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      psii = val
      ierr = 0
      END FUNCTION set_bse_psii


      FUNCTION get_bse_acc2(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = acc2
      ierr = 0
      END FUNCTION get_bse_acc2


      FUNCTION set_bse_acc2(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      acc2 = val
      ierr = 0
      END FUNCTION set_bse_acc2


      FUNCTION get_bse_epsnov(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = epsnov
      ierr = 0
      END FUNCTION get_bse_epsnov


      FUNCTION set_bse_epsnov(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      epsnov = val
      ierr = 0
      END FUNCTION set_bse_epsnov


      FUNCTION get_bse_eddfac(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = eddfac
      ierr = 0
      END FUNCTION get_bse_eddfac


      FUNCTION set_bse_eddfac(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      eddfac = val
      ierr = 0
      END FUNCTION set_bse_eddfac


      FUNCTION get_bse_gamm1(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = gamm1
      ierr = 0
      END FUNCTION get_bse_gamm1


      FUNCTION set_bse_gamm1(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      gamm1 = val
      ierr = 0
      END FUNCTION set_bse_gamm1


      FUNCTION get_bse_fmrg(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(OUT) :: val
      INTEGER :: ierr
      val = fmrg
      ierr = 0
      END FUNCTION get_bse_fmrg


      FUNCTION set_bse_fmrg(val) RESULT(ierr)
      INCLUDE 'common6.h'
      DOUBLE PRECISION, INTENT(IN) :: val
      INTEGER :: ierr
      fmrg = val
      ierr = 0
      END FUNCTION set_bse_fmrg


* ---------------------------------------------------------------------------
* Group B: 10 INTEGER parameters
* ---------------------------------------------------------------------------

      FUNCTION get_bse_ceflag(val) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(OUT) :: val
      INTEGER :: ierr
      val = ceflag
      ierr = 0
      END FUNCTION get_bse_ceflag


      FUNCTION set_bse_ceflag(val) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN) :: val
      INTEGER :: ierr
      ceflag = val
      ierr = 0
      END FUNCTION set_bse_ceflag


      FUNCTION get_bse_tflag(val) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(OUT) :: val
      INTEGER :: ierr
      val = tflag
      ierr = 0
      END FUNCTION get_bse_tflag


      FUNCTION set_bse_tflag(val) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN) :: val
      INTEGER :: ierr
      tflag = val
      ierr = 0
      END FUNCTION set_bse_tflag


      FUNCTION get_bse_ifflag(val) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(OUT) :: val
      INTEGER :: ierr
      val = ifflag
      ierr = 0
      END FUNCTION get_bse_ifflag


      FUNCTION set_bse_ifflag(val) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN) :: val
      INTEGER :: ierr
      ifflag = val
      ierr = 0
      END FUNCTION set_bse_ifflag


      FUNCTION get_bse_wdflag(val) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(OUT) :: val
      INTEGER :: ierr
      val = wdflag
      ierr = 0
      END FUNCTION get_bse_wdflag


      FUNCTION set_bse_wdflag(val) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN) :: val
      INTEGER :: ierr
      wdflag = val
      ierr = 0
      END FUNCTION set_bse_wdflag


      FUNCTION get_bse_bhflag(val) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(OUT) :: val
      INTEGER :: ierr
      val = bhflag
      ierr = 0
      END FUNCTION get_bse_bhflag


      FUNCTION set_bse_bhflag(val) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN) :: val
      INTEGER :: ierr
      bhflag = val
      ierr = 0
      END FUNCTION set_bse_bhflag


      FUNCTION get_bse_nsflag(val) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(OUT) :: val
      INTEGER :: ierr
      val = nsflag
      ierr = 0
      END FUNCTION get_bse_nsflag


      FUNCTION set_bse_nsflag(val) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN) :: val
      INTEGER :: ierr
      nsflag = val
      ierr = 0
      END FUNCTION set_bse_nsflag


      FUNCTION get_bse_psflag(val) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(OUT) :: val
      INTEGER :: ierr
      val = psflag
      ierr = 0
      END FUNCTION get_bse_psflag


      FUNCTION set_bse_psflag(val) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN) :: val
      INTEGER :: ierr
      psflag = val
      ierr = 0
      END FUNCTION set_bse_psflag


      FUNCTION get_bse_kmech(val) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(OUT) :: val
      INTEGER :: ierr
      val = kmech
      ierr = 0
      END FUNCTION get_bse_kmech


      FUNCTION set_bse_kmech(val) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN) :: val
      INTEGER :: ierr
      kmech = val
      ierr = 0
      END FUNCTION set_bse_kmech


      FUNCTION get_bse_ecflag(val) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(OUT) :: val
      INTEGER :: ierr
      val = ecflag
      ierr = 0
      END FUNCTION get_bse_ecflag


      FUNCTION set_bse_ecflag(val) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN) :: val
      INTEGER :: ierr
      ecflag = val
      ierr = 0
      END FUNCTION set_bse_ecflag


      FUNCTION get_bse_edflag(val) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(OUT) :: val
      INTEGER :: ierr
      val = edflag
      ierr = 0
      END FUNCTION get_bse_edflag


      FUNCTION set_bse_edflag(val) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN) :: val
      INTEGER :: ierr
      edflag = val
      ierr = 0
      END FUNCTION set_bse_edflag


* ===========================================================================
* External gravity (used by AMUSE bridge couplings)
* ===========================================================================

      FUNCTION get_gravity_at_point(eps1, x1, y1, z1,
     &                              fx, fy, fz, npts) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN) :: npts
      DOUBLE PRECISION, INTENT(IN)  :: eps1(npts)
      DOUBLE PRECISION, INTENT(IN)  :: x1(npts), y1(npts), z1(npts)
      DOUBLE PRECISION, INTENT(OUT) :: fx(npts), fy(npts), fz(npts)
      INTEGER :: ierr, i, j
      DOUBLE PRECISION :: dx, dy, dz, r2, r2i, ri, mr3i

      DO i = 1, npts
          fx(i) = 0.0d0
          fy(i) = 0.0d0
          fz(i) = 0.0d0
      END DO
      DO i = 1, npts
          DO j = 1, N
              dx = X(1, j) - x1(i)
              dy = X(2, j) - y1(i)
              dz = X(3, j) - z1(i)
              r2  = dx*dx + dy*dy + dz*dz + eps1(i)
              r2i = 1.0d0 / r2
              ri  = SQRT(r2i)
              mr3i = BODY(j) * ri * r2i
              fx(i) = fx(i) + mr3i * dx
              fy(i) = fy(i) + mr3i * dy
              fz(i) = fz(i) + mr3i * dz
          END DO
      END DO
      ierr = 0
      END FUNCTION get_gravity_at_point


      FUNCTION get_potential_at_point(eps1, x1, y1, z1,
     &                                phi, npts) RESULT(ierr)
      INCLUDE 'common6.h'
      INTEGER, INTENT(IN) :: npts
      DOUBLE PRECISION, INTENT(IN)  :: eps1(npts)
      DOUBLE PRECISION, INTENT(IN)  :: x1(npts), y1(npts), z1(npts)
      DOUBLE PRECISION, INTENT(OUT) :: phi(npts)
      INTEGER :: ierr, i, j
      DOUBLE PRECISION :: dx, dy, dz, r2, ri

      DO i = 1, npts
          phi(i) = 0.0d0
      END DO
      DO i = 1, npts
          DO j = 1, N
              dx = X(1, j) - x1(i)
              dy = X(2, j) - y1(i)
              dz = X(3, j) - z1(i)
              r2 = dx*dx + dy*dy + dz*dz + eps1(i)
              ri = 1.0d0 / SQRT(r2)
              phi(i) = phi(i) - BODY(j) * ri
          END DO
      END DO
      ierr = 0
      END FUNCTION get_potential_at_point

      END MODULE Nbody7Interface
