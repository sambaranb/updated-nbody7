AC_CONFIG_AUX_DIR([shared])
AC_CONFIG_MACRO_DIRS([shared/m4])

AC_INIT([amuse-nbody7], [1.0])

AMUSE_LOG_ENVVARS()

# Pick up Conda env or virtualenv
AMUSE_VENV()

# Set CPU_COUNT
AX_COUNT_CPUS()
AC_SUBST(CPU_COUNT)

# Find the compiler(s)
AC_PROG_FC()
AX_FC_ISO_C_BINDING()
AC_PROG_CXX()

# Find tools for creating static libraries
AC_CHECK_TARGET_TOOL([AR], [ar], [:])
AS_IF([test $AR = ":"], [AC_MSG_ERROR([ar command not found.])])

AC_PROG_RANLIB()
AS_IF([test $RANLIB = ":"], [AC_MSG_ERROR([ranlib command not found.])])


# Find tools to download and unpack with
AMUSE_DOWNLOAD()
AC_CHECK_TOOL(TAR, tar)
AC_CHECK_TOOL(GUNZIP, gunzip)

# Find AMUSE libraries
AMUSE_LIB_STOPCOND()
AMUSE_LIB_FORSOCKETS()

# Find external dependencies
AC_LANG_PUSH([Fortran])

##### MPI is always needed to build the worker #####
AX_MPI()

##### Remove anything not needed by the code #####

AC_OPENMP()
AC_SUBST(OPENMP_FCFLAGS)

# CUDA detection for the gpu-variant wheel (amuse-nbody7-gpu).
# Gracefully reports "CUDA not found" on hosts without nvcc + libcudart;
# the gpu wheel build will fail at compile time on such hosts but cpu/sse/avx
# wheels remain unaffected. CUDA_TK env var can override autodetection.
AMUSE_CUDA()


##### End remove anything not needed by the code #####

AC_LANG_POP([Fortran])

# Generate output
AC_CONFIG_FILES([config.mk])
AC_OUTPUT

