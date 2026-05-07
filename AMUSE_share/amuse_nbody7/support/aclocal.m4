# generated automatically by aclocal 1.16.5 -*- Autoconf -*-

# Copyright (C) 1996-2021 Free Software Foundation, Inc.

# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.

m4_ifndef([AC_CONFIG_MACRO_DIRS], [m4_defun([_AM_CONFIG_MACRO_DIRS], [])m4_defun([AC_CONFIG_MACRO_DIRS], [_AM_CONFIG_MACRO_DIRS($@)])])
# AMUSE_CUDA_VERIFY_HEADERS(PATH)
#
# Checks that the given path contains a CUDA installation with include/cuda.h.
#
# Sets amuse_cuda_verify_FLAGS to -I/path/to/cuda/include if successful.
AC_DEFUN([AMUSE_CUDA_VERIFY_HEADERS], [
    AS_IF([test "x$1" = x], [
        amuse_cuda_verify_msg="for cuda.h"
        amuse_cuda_flags=
    ], [
        amuse_cuda_verify_msg="for cuda.h in $1/include"
        amuse_cuda_flags="-I$1/include"
    ])

    AC_MSG_CHECKING([$amuse_cuda_verify_msg])
    amuse_cuda_verify_FLAGS=

    ax_save_CFLAGS="$CFLAGS"

    CFLAGS="$amuse_cuda_flags"
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([#include <cuda.h>])], [
        amuse_cuda_verify_FLAGS="$CFLAGS"
        amuse_cuda_verify_headers_found="yes"
        AC_MSG_RESULT([yes])
    ], [
        AC_MSG_RESULT([no])
    ])

    CFLAGS="$ax_save_CFLAGS"
])

# AMUSE_CUDA_VERIFY_LIBS(PATH)
#
# Checks that the given path contains a CUDA installation with lib/libcudart.so or
# lib64/libcudart.so.
#
# Sets amuse_cuda_verify_LDFLAGS to -L/path/to/cuda/lib{64} if successful.
AC_DEFUN([AMUSE_CUDA_VERIFY_LIBS], [
    AS_IF([test "x$1" = x], [
        amuse_cuda_verify_msg="for libcudart"
        amuse_cuda_ldflags=
    ], [
        amuse_cuda_verify_msg="for libcudart in $1/lib"
        amuse_cuda_verify_msg64="for libcudart in $1/lib64"
        amuse_cuda_ldflags="-L$1/lib"
        amuse_cuda_ldflags64="-L$1/lib64"
    ])

    AC_MSG_CHECKING([$amuse_cuda_verify_msg])
    amuse_cuda_verify_LDFLAGS=

    ax_save_LDFLAGS="$LDFLAGS"
    ax_save_LIBS="$LIBS"

    LDFLAGS="$amuse_cuda_ldflags"
    LIBS="-lcudart"
    AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <cuda.h>], [cudaFree(0);])], [
        amuse_cuda_verify_LDFLAGS="$LDFLAGS"
        amuse_cuda_verify_libs_found="yes"
        AC_MSG_RESULT([yes])
    ], [
        AC_MSG_RESULT([no])
    ])

    AS_IF([test "x$amuse_cuda_verify_libs_found" = x], [
        AC_MSG_CHECKING([$amuse_cuda_verify_msg64])
        LDFLAGS="$amuse_cuda_ldflags64"
        AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <cuda.h>], [cudaFree(0);])], [
            amuse_cuda_verify_LDFLAGS="$LDFLAGS"
            amuse_cuda_verify_libs_found="yes"
            AC_MSG_RESULT([yes])
        ], [
            AC_MSG_RESULT([no])
        ])
    ])

    LIBS="$ax_save_LIBS"
    LDFLAGS="$ax_save_LDFLAGS"
])

# AMUSE_CUDA_VERIFY_SET_VARS(CUDA_TK)
#
# Sets CUDA_TK, NVCC, CUDA_FLAGS and CUDA_LDFLAGS from the internal values, if all
# of them have been found.
AC_DEFUN([AMUSE_CUDA_VERIFY_SET_VARS], [
    AS_IF(
        [test \(x"$amuse_cuda_verify_NVCC" != x\) -a \(x"$amuse_cuda_verify_headers_found" != x\) -a \(x"$amuse_cuda_verify_libs_found" != x\)],
        [
            CUDA_TK="$1"
            NVCC="$amuse_cuda_verify_NVCC"
            CUDA_FLAGS="$amuse_cuda_verify_FLAGS"
            CUDA_LDFLAGS="$amuse_cuda_verify_LDFLAGS"
        ]
    )
])

# AMUSE_CUDA_VERIFY(PATH)
#
# Checks that the given path contains a CUDA installation with bin/nvcc,
# lib/libcudart.so, and include/cuda.h.
#
# Sets CUDA_TK=<path> if successful, and if successful also sets NVCC to the location
# of nvcc, CUDA_LDFLAGS to -L/path/to/cuda/libs, and CUDA_FLAGS to
# -I/path/to/cuda/include.
AC_DEFUN([AMUSE_CUDA_VERIFY], [
    amuse_cuda_verify_headers_found=
    amuse_cuda_verify_libs_found=

    AS_IF([test -d $1], [
        AC_PATH_PROG([amuse_cuda_verify_NVCC], [nvcc], [], [$1/bin])

        AS_IF([test x"$amuse_cuda_verify_NVCC" != x], [
            AMUSE_CUDA_VERIFY_HEADERS([$1])
        ])

        AS_IF([test x"$amuse_cuda_verify_headers_found" != x], [
            AMUSE_CUDA_VERIFY_LIBS([$1])
        ])

        AMUSE_CUDA_VERIFY_SET_VARS([$1])
    ], [
        AC_MSG_NOTICE([$1 does not exist or is not a directory, CUDA not found there])
    ])
])

# AMUSE_CUDA_VERIFY_DEFAULT()
#
# Checks for NVCC on the PATH, then finds the CUDA directory from there.
#
# Sets CUDA_TK=<path> if successful, and if successful also sets NVCC to the location
# of nvcc, CUDA_LDFLAGS to -L/path/to/cuda/libs, and CUDA_FLAGS to
# -I/path/to/cuda/include.
AC_DEFUN([AMUSE_CUDA_VERIFY_DEFAULT], [
    AC_PATH_PROG([amuse_cuda_verify_NVCC], [nvcc])

    AS_IF([test x"$amuse_cuda_verify_NVCC" != x], [
        # Got nvcc, verify that we have the rest too
        amuse_cuda_verify_ctk_rel=$(dirname -- "$amuse_cuda_verify_NVCC")/..
        # Canonicalise path portably so that it looks nicer
        amuse_cuda_verify_CUDA_TK=$(test -d "$amuse_cuda_verify_ctk_rel" && CDPATH= cd -P -- "$amuse_cuda_verify_ctk_rel" && pwd -P)

        AS_IF([test x"$amuse_cuda_verify_NVCC" != x], [
            AMUSE_CUDA_VERIFY_HEADERS()
        ])

        AS_IF([test x"$amuse_cuda_verify_headers_found" != x], [
            AMUSE_CUDA_VERIFY_LIBS()
        ])

        AMUSE_CUDA_VERIFY_SET_VARS([$amuse_cuda_verify_CUDA_TK])
    ], [
        AC_MSG_NOTICE([Could not find CUDA via PATH])
    ])
])

# AMUSE_CUDA
# ----------------------------------------------------------
# set up for CUDA
#
# will set:
#
# CUDA_TK
# NVCC
#
# CUDA_LDFLAGS  (with -L/path/to/cuda/lib)
# CUDA_FLAGS    (with -I/path/to/cuda/include)
#
# and call AC_SUBST on them.
#
AC_DEFUN([AMUSE_CUDA], [
    AC_LANG_PUSH([C])
    AS_IF([test x"$CUDA_TK" != x], [
        # User set CUDA_TK, verify and use it if it works
        AMUSE_CUDA_VERIFY($CUDA_TK)

        AS_IF([test x"$NVCC" = x], [
            AC_MSG_ERROR([CUDA_TK is set, but there is no nvcc in $CUDA_TK/bin. Please set CUDA_TK to the CUDA installation directory, or unset it to autodetect.])
        ])
    ], [
        # CUDA_TK not set, try to discover CUDA via PATH
        AMUSE_CUDA_VERIFY_DEFAULT()

        # Not in PATH, try default directory
        AS_IF([test x"$CUDA_TK" = x], [
            AMUSE_CUDA_VERIFY([/usr/local/cuda])
        ])

        # Try /usr/local/cuda-#.#, but only if there's exactly one match
        AS_IF([test x"$CUDA_TK" = x], [
            # List CUDA installations
            amuse_cuda_installs=$(ls -1 -d /usr/local/cuda-*.* 2>/dev/null)
            # If there's more than one, the above will have newlines, and change when we do this
            amuse_cuda_installs2=$(echo "x${amuse_cuda_installs}" | tr -d '[:space:]')
            AS_IF([test "x${amuse_cuda_installs}" != x], [
                AS_IF([test "x${amuse_cuda_installs}" = "${amuse_cuda_installs2}"], [
                    # Here, there's exactly one match
                    AMUSE_CUDA_VERIFY([$amuse_cuda_installs])
                ], [
                    AC_MSG_NOTICE([Multiple CUDA installations found in /usr/local, without a /usr/local/cuda symlink. Please set CUDA_TK to specify which one you want to use.])
                ])
            ], [
                AC_MSG_NOTICE([No CUDA installations found in /usr/local])
            ])
        ])

        # Try some other locations
        for dir in /opt/cuda /usr/local/cuda/cuda /opt/cuda/cuda ; do
            AS_IF([test x"$CUDA_TK" = x], [
                AMUSE_CUDA_VERIFY([$dir])
            ])
        done

        # These directories were checked by the old macro, but they're highly obsolete
        # so we're not trying them anymore:
        # /usr/lib/nvidia, /usr/include/nvidia, /usr/lib/nvidia-current, /usr/include/nvidia-current

        AS_IF([test x"$CUDA_TK" = x], [
            AC_MSG_NOTICE([CUDA not found. Set CUDA_TK if you have it in an odd location.])
        ])
    ])

    AS_IF([test x"$CUDA_TK" != x], [
        AC_MSG_NOTICE([CUDA found at $CUDA_TK])
        AC_SUBST(CUDA_TK)
        AC_SUBST(NVCC)
        AC_SUBST(CUDA_FLAGS)
        AC_SUBST(CUDA_LDFLAGS)
    ])

    AC_LANG_POP([C])
])


# Helper macros for detecting download tools
#
# Some of the community codes aren't included with AMUSE, but are downloaded at build
# time. We need a tool for that, and here is where we find one.
#
# AMUSE_DOWNLOAD()
#
# Searches for a download tool.
#
# This macro tries to find a downloader and sets DOWNLOAD to a command that will take a
# URL and download its contents to standard output. This makes it easier to write that
# output to a file with a known name, which you usually want in a Makefile.
#
# This used to support both wget and curl, but we had problems with wget with BitBucket
# and a MESA server, so that's been removed and it's curl-only now.
#
# To download a file, use $(DOWNLOAD) https://example.com >example.html
#
AC_DEFUN([AMUSE_DOWNLOAD], [
    AC_CHECK_TOOL(CURL, curl)

    AC_MSG_CHECKING([for a tool to download files with])
    if test "x$CURL" != "x"
    then
        DOWNLOAD="$CURL -L"
        AC_MSG_RESULT([yes])
    else
        AC_MSG_RESULT([no])
    fi

    AC_SUBST([DOWNLOAD])
])


# Helper macros for detecting AMUSE libraries.
#
# See https://stackoverflow.com/questions/10220946/pkg-check-modules-considered-harmful
# for some background on why we're not using PKG_CHECK_MODULES directly. Do note that
# PKG_CHECK_MODULES checks for pkg-config and gives a useful error these days, so that
# that page is partially outdated.

# AMUSE_LIB(prefix, module, library, function)
#
# Searches for an AMUSE library and sets ${prefix}_CFLAGS, ${prefix}_LIBS to the
# appropriate values if it is found.
#
# prefix: prefix for the variables to be set
# module: name of the pkg-config module to search for if needed
# library: name of the library to be searched for
# function: name of a function in the library to use for the link check
AC_DEFUN([AMUSE_LIB], [
    amuse_save_LIBS="$LIBS"
    amuse_save_LIB_CFLAGS="$[$1][_CFLAGS]"
    amuse_save_LIB_LIBS="$[$1][_LIBS]"
    amuse_save_PKG_CONFIG_PATH="$PKG_CONFIG_PATH"

    # If we have an active virtualenv, make sure pkg-config searches it
    if test "a${VIRTUAL_ENV}" != "a"
    then
        PKG_CONFIG_PATH="${VIRTUAL_ENV}/lib/pkgconfig:${PKG_CONFIG_PATH}"
    fi

    # All AMUSE libs export C symbols
    AC_LANG_PUSH([C])

    # Search for the library, first directly then fall back to pkg-config
    AC_SEARCH_LIBS([$4], [$3], [
        FOUND_$1="yes"
        $1_LIBS="$LIBS"
        $1_CFLAGS=""
    ], [
        PKG_CHECK_MODULES([$1], [$2], [
            FOUND_$1="yes"
        ], [
            FOUND_$1="no"
        ])
    ])

    AC_LANG_POP([C])

    PKG_CONFIG_PATH="$amuse_save_PKG_CONFIG_PATH"
    LIBS="$amuse_save_LIBS"

    # If we have an active CONDA environment, assume that the lib is coming from
    # there and add an additional flag so that .mod files can be found. Only really
    # needed for stopcond and forsockets, and hopefully conda-forge will give us a
    # better solution soon.
    if test "${FOUND_$1}" == "yes" -a "x$CONDA_PREFIX" != "x"
    then
        $1_CFLAGS="${$1_CFLAGS} -I${CONDA_PREFIX}/include"
    fi

    # If the user overrode the variables, go with what they set instead of
    # what we just detected.
    AS_IF([test "x$amuse_save_LIB_CFLAGS" != "x"], [
        $1_CFLAGS="$amuse_save_LIB_CFLAGS"
    ])
    AS_IF([test "x$amuse_save_LIB_LIBS" != "x"], [
        $1_LIBS="$amuse_save_LIB_LIBS"
    ])

    AC_SUBST([FOUND_][$1])
    AC_SUBST([$1][_CFLAGS])
    AC_SUBST([$1][_LIBS])
])


# AMUSE_LIB_STOPCOND()
#
# Searches for the AMUSE stopping conditions library and sets STOPCOND_CFLAGS
# and STOPCOND_LIBS to the appropriate values if it is found.
AC_DEFUN([AMUSE_LIB_STOPCOND], [
    AMUSE_LIB([STOPCOND], [stopcond], [stopcond], [is_condition_enabled])
])


# AMUSE_LIB_STOPCONDMPI()
#
# Searches for the AMUSE stopping conditions library and sets STOPCONDMPI_CFLAGS
# and STOPCONDMPI_LIBS to the appropriate values if it is found.
#
# We use a different function here, to avoid getting a cached value from
# AMUSE_LIB_STOPCOND if both are used. Seems like it indexes by function name.
AC_DEFUN([AMUSE_LIB_STOPCONDMPI], [
    AMUSE_LIB([STOPCONDMPI], [stopcondmpi], [stopcondmpi], [get_set_conditions_])
])


# AMUSE_LIB_AMUSE_MPI()
#
# Searches for the AMUSE MPI helper library and sets AMUSE_MPI_CFLAGS and
# AMUSE_MPI_LIBS to the appropriate values if it is found.
AC_DEFUN([AMUSE_LIB_AMUSE_MPI], [
    AMUSE_LIB([AMUSE_MPI], [amuse_mpi], [amuse_mpi], [get_comm_world])
])


# AMUSE_LIB_FORSOCKETS()
#
# Searches for the AMUSE forsockets library and sets FORSOCKETS_CFLAGS and
# FORSOCKETS_LIBS to the appropriate values if it is found.
AC_DEFUN([AMUSE_LIB_FORSOCKETS], [
    AMUSE_LIB([FORSOCKETS], [forsockets], [forsockets], [forsockets_close])
])


# AMUSE_LIB_SIMPLE_HASH()
#
# Searches for the AMUSE simple hash library and sets SIMPLE_HASH_CFLAGS and
# SIMPLE_HASH_LIBS to the appropriate values if it is found.
AC_DEFUN([AMUSE_LIB_SIMPLE_HASH], [
    AMUSE_LIB([SIMPLE_HASH], [simple_hash], [simple_hash], [init_hash])
])


# AMUSE_LIB_G6LIB()
#
# Searches for the g6 library and sets G6LIB_CFLAGS and G6LIB_LIBS to
# the appropriate values if it is found.
AC_DEFUN([AMUSE_LIB_G6LIB], [
    AMUSE_LIB([G6LIB], [g6lib], [g6], [g6_npipes])
])


# AMUSE_LIB_SAPPORO_LIGHT()
#
# Searches for the Sapporo light library and sets SAPPORO_LIGHT_CFLAGS and
# SAPPORO_LIGHT_LIBS to the appropriate values if it is found.
AC_DEFUN([AMUSE_LIB_SAPPORO_LIGHT], [
    AMUSE_LIB([SAPPORO_LIGHT], [sapporo_light], [sapporo], [get_device_count])
])


# Macro for logging important environment variables

# AMUSE_LOG_ENVVARS()
#
AC_DEFUN([AMUSE_LOG_ENVVARS], [
    amuse_le_vars="SHELL CC CXX FC F77 MPICC MPICXX MPIFC CFLAGS CXXFLAGS CPPFLAGS LDFLAGS LIBS OMPI_CC OMPI_CXX OMPI_FC OMPI_CPPFLAGS OMPI_LDFLAGS OMPI_LIBS"

    AS_ECHO(["=== Shell environment ==="])

    for amuse_le_var in ${amuse_le_vars}
    do
        eval "amuse_le_val=\${$amuse_le_var}"
        AS_ECHO(["${amuse_le_var}=${amuse_le_val}"])
    done
    AS_ECHO(["=== End shell environment ==="])
])


# Helper macro for installing into Conda envs and virtualenvs.

# AMUSE_VENV()
#
# Detect the active Conda env or virtualenv and extend LDFLAGS and
# PKG_CONFIG_PATH to point to it. This then makes it possible to detect the AMUSE
# libraries as installed in the environment and link to them.
#
AC_DEFUN([AMUSE_VENV], [
    AS_IF([test "x$VIRTUAL_ENV" != x], [
        CFLAGS="$CFLAGS -I${VIRTUAL_ENV}/include"
        CXXFLAGS="$CXXFLAGS -I${VIRTUAL_ENV}/include"
        FFLAGS="$FFLAGS -I${VIRTUAL_ENV}/include"
        FCFLAGS="$FCFLAGS -I${VIRTUAL_ENV}/include"
        LDFLAGS="$LDFLAGS -L${VIRTUAL_ENV}/lib -Wl,-rpath,${VIRTUAL_ENV}/lib"
        PKG_CONFIG_PATH="$VIRTUAL_ENV/lib/pkgconfig:$PKG_CONFIG_PATH"
    ])

    AS_IF([test "x$CONDA_PREFIX" != x], [
        # Conda does not set FCFLAGS, so we copy from FFLAGS here
        FCFLAGS="$FFLAGS"
        # Conda pkg-config includes this already, but in case we have one from
        # the system...
        PKG_CONFIG_PATH="$PKG_CONFIG_PATH:${CONDA_PREFIX}/lib/pkgconfig"
    ])
    # Needs to be exported or the PKG_CHECK_MODULES macro won't see it
    export PKG_CONFIG_PATH
    AC_SUBST([FFLAGS])
])


# ===========================================================================
#      https://www.gnu.org/software/autoconf-archive/ax_count_cpus.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_COUNT_CPUS([ACTION-IF-DETECTED],[ACTION-IF-NOT-DETECTED])
#
# DESCRIPTION
#
#   Attempt to count the number of logical processor cores (including
#   virtual and HT cores) currently available to use on the machine and
#   place detected value in CPU_COUNT variable.
#
#   On successful detection, ACTION-IF-DETECTED is executed if present. If
#   the detection fails, then ACTION-IF-NOT-DETECTED is triggered. The
#   default ACTION-IF-NOT-DETECTED is to set CPU_COUNT to 1.
#
# LICENSE
#
#   Copyright (c) 2014,2016 Karlson2k (Evgeny Grin) <k2k@narod.ru>
#   Copyright (c) 2012 Brian Aker <brian@tangent.org>
#   Copyright (c) 2008 Michael Paul Bailey <jinxidoru@byu.net>
#   Copyright (c) 2008 Christophe Tournayre <turn3r@users.sourceforge.net>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 22

  AC_DEFUN([AX_COUNT_CPUS],[dnl
      AC_REQUIRE([AC_CANONICAL_HOST])dnl
      AC_REQUIRE([AC_PROG_EGREP])dnl
      AC_MSG_CHECKING([the number of available CPUs])
      CPU_COUNT="0"

      # Try generic methods

      # 'getconf' is POSIX utility, but '_NPROCESSORS_ONLN' and
      # 'NPROCESSORS_ONLN' are platform-specific
      command -v getconf >/dev/null 2>&1 && \
        CPU_COUNT=`getconf _NPROCESSORS_ONLN 2>/dev/null || getconf NPROCESSORS_ONLN 2>/dev/null` || CPU_COUNT="0"
      AS_IF([[test "$CPU_COUNT" -gt "0" 2>/dev/null || ! command -v nproc >/dev/null 2>&1]],[[: # empty]],[dnl
        # 'nproc' is part of GNU Coreutils and is widely available
        CPU_COUNT=`OMP_NUM_THREADS='' nproc 2>/dev/null` || CPU_COUNT=`nproc 2>/dev/null` || CPU_COUNT="0"
      ])dnl

      AS_IF([[test "$CPU_COUNT" -gt "0" 2>/dev/null]],[[: # empty]],[dnl
        # Try platform-specific preferred methods
        AS_CASE([[$host_os]],dnl
          [[*linux*]],[[CPU_COUNT=`lscpu -p 2>/dev/null | $EGREP -e '^@<:@0-9@:>@+,' -c` || CPU_COUNT="0"]],dnl
          [[*darwin*]],[[CPU_COUNT=`sysctl -n hw.logicalcpu 2>/dev/null` || CPU_COUNT="0"]],dnl
          [[freebsd*]],[[command -v sysctl >/dev/null 2>&1 && CPU_COUNT=`sysctl -n kern.smp.cpus 2>/dev/null` || CPU_COUNT="0"]],dnl
          [[netbsd*]], [[command -v sysctl >/dev/null 2>&1 && CPU_COUNT=`sysctl -n hw.ncpuonline 2>/dev/null` || CPU_COUNT="0"]],dnl
          [[solaris*]],[[command -v psrinfo >/dev/null 2>&1 && CPU_COUNT=`psrinfo 2>/dev/null | $EGREP -e '^@<:@0-9@:>@.*on-line' -c 2>/dev/null` || CPU_COUNT="0"]],dnl
          [[mingw*]],[[CPU_COUNT=`ls -qpU1 /proc/registry/HKEY_LOCAL_MACHINE/HARDWARE/DESCRIPTION/System/CentralProcessor/ 2>/dev/null | $EGREP -e '^@<:@0-9@:>@+/' -c` || CPU_COUNT="0"]],dnl
          [[msys*]],[[CPU_COUNT=`ls -qpU1 /proc/registry/HKEY_LOCAL_MACHINE/HARDWARE/DESCRIPTION/System/CentralProcessor/ 2>/dev/null | $EGREP -e '^@<:@0-9@:>@+/' -c` || CPU_COUNT="0"]],dnl
          [[cygwin*]],[[CPU_COUNT=`ls -qpU1 /proc/registry/HKEY_LOCAL_MACHINE/HARDWARE/DESCRIPTION/System/CentralProcessor/ 2>/dev/null | $EGREP -e '^@<:@0-9@:>@+/' -c` || CPU_COUNT="0"]]dnl
        )dnl
      ])dnl

      AS_IF([[test "$CPU_COUNT" -gt "0" 2>/dev/null || ! command -v sysctl >/dev/null 2>&1]],[[: # empty]],[dnl
        # Try less preferred generic method
        # 'hw.ncpu' exist on many platforms, but not on GNU/Linux
        CPU_COUNT=`sysctl -n hw.ncpu 2>/dev/null` || CPU_COUNT="0"
      ])dnl

      AS_IF([[test "$CPU_COUNT" -gt "0" 2>/dev/null]],[[: # empty]],[dnl
      # Try platform-specific fallback methods
      # They can be less accurate and slower then preferred methods
        AS_CASE([[$host_os]],dnl
          [[*linux*]],[[CPU_COUNT=`$EGREP -e '^processor' -c /proc/cpuinfo 2>/dev/null` || CPU_COUNT="0"]],dnl
          [[*darwin*]],[[CPU_COUNT=`system_profiler SPHardwareDataType 2>/dev/null | $EGREP -i -e 'number of cores:'|cut -d : -f 2 -s|tr -d ' '` || CPU_COUNT="0"]],dnl
          [[freebsd*]],[[CPU_COUNT=`dmesg 2>/dev/null| $EGREP -e '^cpu@<:@0-9@:>@+: '|sort -u|$EGREP -e '^' -c` || CPU_COUNT="0"]],dnl
          [[netbsd*]], [[CPU_COUNT=`command -v cpuctl >/dev/null 2>&1 && cpuctl list 2>/dev/null| $EGREP -e '^@<:@0-9@:>@+ .* online ' -c` || \
                           CPU_COUNT=`dmesg 2>/dev/null| $EGREP -e '^cpu@<:@0-9@:>@+ at'|sort -u|$EGREP -e '^' -c` || CPU_COUNT="0"]],dnl
          [[solaris*]],[[command -v kstat >/dev/null 2>&1 && CPU_COUNT=`kstat -m cpu_info -s state -p 2>/dev/null | $EGREP -c -e 'on-line'` || \
                           CPU_COUNT=`kstat -m cpu_info 2>/dev/null | $EGREP -c -e 'module: cpu_info'` || CPU_COUNT="0"]],dnl
          [[mingw*]],[AS_IF([[CPU_COUNT=`reg query 'HKLM\\Hardware\\Description\\System\\CentralProcessor' 2>/dev/null | $EGREP -e '\\\\@<:@0-9@:>@+$' -c`]],dnl
                        [[: # empty]],[[test "$NUMBER_OF_PROCESSORS" -gt "0" 2>/dev/null && CPU_COUNT="$NUMBER_OF_PROCESSORS"]])],dnl
          [[msys*]],[[test "$NUMBER_OF_PROCESSORS" -gt "0" 2>/dev/null && CPU_COUNT="$NUMBER_OF_PROCESSORS"]],dnl
          [[cygwin*]],[[test "$NUMBER_OF_PROCESSORS" -gt "0" 2>/dev/null && CPU_COUNT="$NUMBER_OF_PROCESSORS"]]dnl
        )dnl
      ])dnl

      AS_IF([[test "x$CPU_COUNT" != "x0" && test "$CPU_COUNT" -gt 0 2>/dev/null]],[dnl
          AC_MSG_RESULT([[$CPU_COUNT]])
          m4_ifvaln([$1],[$1],)dnl
        ],[dnl
          m4_ifval([$2],[dnl
            AS_UNSET([[CPU_COUNT]])
            AC_MSG_RESULT([[unable to detect]])
            $2
          ], [dnl
            CPU_COUNT="1"
            AC_MSG_RESULT([[unable to detect (assuming 1)]])
          ])dnl
        ])dnl
      ])dnl

# ===========================================================================
#          http://www.gnu.org/software/autoconf-archive/ax_mpi.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_MPI([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro tries to find out how to compile programs that use MPI
#   (Message Passing Interface), a standard API for parallel process
#   communication (see http://www-unix.mcs.anl.gov/mpi/)
#
#   On success, it sets the MPICC, MPICXX, MPIF77, or MPIFC output variable
#   to the name of the MPI compiler, depending upon the current language.
#   (This may just be $CC/$CXX/$F77/$FC, but is more often something like
#   mpicc/mpiCC/mpif77/mpif90.) It also sets MPILIBS to any libraries that
#   are needed for linking MPI (e.g. -lmpi or -lfmpi, if a special
#   MPICC/MPICXX/MPIF77/MPIFC was not found).
#
#   If you want to compile everything with MPI, you should use something
#   like this for C:
#
#     if test -z "$CC" && test -n "$MPICC"; then
#       CC="$MPICC"
#     fi
#     AC_PROG_CC
#     AX_MPI
#     CC="$MPICC"
#     LIBS="$MPILIBS $LIBS"
#
#   and similar for C++ (change all instances of CC to CXX), Fortran 77
#   (with F77 instead of CC) or Fortran (with FC instead of CC).
#
#   NOTE: The above assumes that you will use $CC (or whatever) for linking
#   as well as for compiling. (This is the default for automake and most
#   Makefiles.)
#
#   The user can force a particular library/compiler by setting the
#   MPICC/MPICXX/MPIF77/MPIFC and/or MPILIBS environment variables.
#
#   ACTION-IF-FOUND is a list of shell commands to run if an MPI library is
#   found, and ACTION-IF-NOT-FOUND is a list of commands to run if it is not
#   found. If ACTION-IF-FOUND is not specified, the default action will
#   define HAVE_MPI.
#
# LICENSE
#
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#   Copyright (c) 2008 Julian C. Cummings <cummings@cacr.caltech.edu>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 7

AU_ALIAS([ACX_MPI], [AX_MPI])
AC_DEFUN([AX_MPI], [
AC_PREREQ([2.71]) dnl for AC_LANG_CASE

AC_LANG_CASE([C], [
	AC_REQUIRE([AC_PROG_CC])
	AC_ARG_VAR(MPICC,[MPI C compiler command])
	AC_CHECK_PROGS(MPICC, mpicc hcc mpxlc_r mpxlc mpcc cmpicc, $CC)
	ax_mpi_save_CC="$CC"
	CC="$MPICC"
	AC_SUBST(MPICC)
 	AC_MSG_CHECKING([checking MPI C flags])
 	ax_mpi_c_flags="`$MPICC -showme:compile 2>/dev/null| cut -d\  -f2-`"
 	ax_mpi_c_libs="`$MPICC -showme:link 2>/dev/null| cut -d\  -f2-`"
        AS_IF([test "x$ax_mpi_c_flags" = "x"],[
          ax_mpi_c_flags="`$MPICC -show -c 2>/dev/null| cut -d\  -f2-|sed s/-c\ //`"
          ax_mpi_c_libs="`$MPICC -show 2>/dev/null| cut -d\  -f2-`"
	   AS_IF([test "x$ax_mpi_c_flags" = "x"],[AC_MSG_RESULT([could not determine c flags from show functions])],[AC_MSG_RESULT([flags found])])
          
        ], [
          AC_MSG_RESULT([flags found])])
	MPI_CFLAGS="$ax_mpi_c_flags"
	MPI_CLIBS="$ax_mpi_c_libs"
	AC_SUBST(MPI_CFLAGS)
	AC_SUBST(MPI_CLIBS)

],
[C++], [
	AC_REQUIRE([AC_PROG_CXX])
	AC_ARG_VAR(MPICXX,[MPI C++ compiler command])
	AC_CHECK_PROGS(MPICXX, mpicxx mpiCC mpic++ hcp mpxlC_r mpxlC mpCC cmpic++, $CXX)
	ax_mpi_save_CXX="$CXX"
	CXX="$MPICXX"
	AC_SUBST(MPICXX)

 	AC_MSG_CHECKING([checking MPI C++ flags])
 	ax_mpi_cc_flags="`$MPICXX -showme:compile 2>/dev/null| cut -d\  -f2-`"
 	ax_mpi_cc_libs="`$MPICXX -showme:link 2>/dev/null| cut -d\  -f2-`"
        AS_IF([test "x$ax_mpi_cc_flags" = "x"],[
          ax_mpi_cc_flags="`$MPICXX -show -c 2>/dev/null| cut -d\  -f2-|sed s/-c\ //`"
          ax_mpi_cc_libs="`$MPICXX -show 2>/dev/null| cut -d\  -f2-`"
	   AS_IF([test "x$ax_mpi_cc_flags" = "x"],[AC_MSG_RESULT([could not determine C++ flags from show functions])],[AC_MSG_RESULT([flags found])])
          
        ], [
          AC_MSG_RESULT([flags found])])
	MPI_CXXFLAGS="$ax_mpi_cc_flags"
	MPI_CXXLIBS="$ax_mpi_cc_libs"
	AC_SUBST(MPI_CXXFLAGS)
	AC_SUBST(MPI_CXXLIBS)

],
[Fortran 77], [
	AC_REQUIRE([AC_PROG_F77])
	AC_ARG_VAR(MPIF77,[MPI Fortran 77 compiler command])
	AC_CHECK_PROGS(MPIF77, mpif77 hf77 mpxlf_r mpxlf mpf77 cmpifc, $F77)
	ax_mpi_save_F77="$F77"
	F77="$MPIF77"
	AC_SUBST(MPIF77)
],
[Fortran], [
	AC_REQUIRE([AC_PROG_FC])
	AC_ARG_VAR(MPIFC,[MPI Fortran compiler command])
	AC_CHECK_PROGS(MPIFC, mpif90 mpxlf95_r mpxlf90_r mpxlf95 mpxlf90 mpf90 cmpif90c, $FC)
	ax_mpi_save_FC="$FC"
	FC="$MPIFC"
	AC_SUBST(MPIFC)
	AC_MSG_CHECKING([checking MPI Fortran flags])
 	ax_mpi_fc_flags="`$MPIFC -showme:compile 2>/dev/null| cut -d\  -f2-`"
 	ax_mpi_fc_libs="`$MPIFC -showme:link 2>/dev/null| cut -d\  -f2-`"
        AS_IF([test "x$ax_mpi_fc_flags" = "x"],[
          ax_mpi_fc_flags="`$MPIFC -show -c 2>/dev/null| cut -d\  -f2-|sed s/-c\ //`"
          ax_mpi_fc_libs="`$MPIFC -show 2>/dev/null| cut -d\  -f2-`"
	   AS_IF([test "x$ax_mpi_fc_flags" = "x"],[AC_MSG_RESULT([could not determine c flags from show functions])],[AC_MSG_RESULT([flags found])])
          
        ], [
          AC_MSG_RESULT([flags found])])
	MPI_FCFLAGS="$ax_mpi_fc_flags"
	MPI_FCLIBS="$ax_mpi_fc_libs"
	AC_SUBST(MPI_FCFLAGS)
	AC_SUBST(MPI_FCLIBS)
])

if test x = x"$MPILIBS"; then
	AC_LANG_CASE([C], [AC_CHECK_FUNC(MPI_Init, [MPILIBS=" "])],
		[C++], [AC_CHECK_FUNC(MPI_Init, [MPILIBS=" "])],
		[Fortran 77], [AC_MSG_CHECKING([for MPI_Init])
			AC_LINK_IFELSE([AC_LANG_PROGRAM([],[      call MPI_Init])],[MPILIBS=" "
				AC_MSG_RESULT(yes)], [AC_MSG_RESULT(no)])],
		[Fortran], [AC_MSG_CHECKING([for MPI_Init])
			AC_LINK_IFELSE([AC_LANG_PROGRAM([],[      call MPI_Init])],[MPILIBS=" "
				AC_MSG_RESULT(yes)], [AC_MSG_RESULT(no)])])
fi
AC_LANG_CASE([Fortran 77], [
	if test x = x"$MPILIBS"; then
		AC_CHECK_LIB(fmpi, MPI_Init, [MPILIBS="-lfmpi"])
	fi
	if test x = x"$MPILIBS"; then
		AC_CHECK_LIB(fmpich, MPI_Init, [MPILIBS="-lfmpich"])
	fi
],
[Fortran], [
	if test x = x"$MPILIBS"; then
		AC_CHECK_LIB(fmpi, MPI_Init, [MPILIBS="-lfmpi"])
	fi
	if test x = x"$MPILIBS"; then
		AC_CHECK_LIB(mpichf90, MPI_Init, [MPILIBS="-lmpichf90"])
	fi

])
if test x = x"$MPILIBS"; then
	AC_CHECK_LIB(mpi, MPI_Init, [MPILIBS="-lmpi"])
fi
if test x = x"$MPILIBS"; then
	AC_CHECK_LIB(mpich, MPI_Init, [MPILIBS="-lmpich"])
fi

dnl We have to use AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[]], [[]])],[],[]) and not AC_CHECK_HEADER because the
dnl latter uses $CPP, not $CC (which may be mpicc).
AC_LANG_CASE([C], [if test x != x"$MPILIBS"; then
	AC_MSG_CHECKING([for mpi.h])
	AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]], [[]])],[AC_MSG_RESULT(yes)],[MPILIBS=""
		AC_MSG_RESULT(no)])
fi],
[C++], [if test x != x"$MPILIBS"; then
	AC_MSG_CHECKING([for mpi.h])
	AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]], [[]])],[AC_MSG_RESULT(yes)],[MPILIBS=""
		AC_MSG_RESULT(no)])
fi],
[Fortran 77], [if test x != x"$MPILIBS"; then
	AC_MSG_CHECKING([for mpif.h])
	AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[      include 'mpif.h'])],[AC_MSG_RESULT(yes)], [MPILIBS=""
		AC_MSG_RESULT(no)])
fi],
[Fortran], [if test x != x"$MPILIBS"; then
	AC_MSG_CHECKING([for mpif.h])
	AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[      include 'mpif.h'])],[AC_MSG_RESULT(yes)], [MPILIBS=""
		AC_MSG_RESULT(no)])
fi])

AC_LANG_CASE([C], [CC="$ax_mpi_save_CC"],
	[C++], [CXX="$ax_mpi_save_CXX"],
	[Fortran 77], [F77="$ax_mpi_save_F77"],
	[Fortran], [FC="$ax_mpi_save_FC"])

AC_SUBST(MPILIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x = x"$MPILIBS"; then
        $2
        :
else
        ifelse([$1],,[AC_DEFINE(HAVE_MPI,1,[Define if you have the MPI library.])],[$1])
        :
fi
])dnl AX_MPI

AC_DEFUN([AX_FC_WORKS],[
AC_LANG_PUSH(Fortran)
AC_MSG_CHECKING([whether the Fortran 90 compiler ($FC $FCFLAGS $LDFLAGS) works])
AC_LINK_IFELSE([
    AC_LANG_SOURCE([
        program conftest
        integer, dimension(10) :: n
        end
    ])
],[
    ax_cv_prog_fc_works="yes"
    AC_MSG_RESULT([$ax_cv_prog_fc_works])
],[
    ax_cv_prog_fc_works="no"
    AC_MSG_WARN([installation or configuration problem: Fortran 90 compiler cannot create executables.])
])
# The intel compiler sometimes generates these work.pc and .pcl files
rm -f work.pc work.pcl
AC_LANG_POP(Fortran)
])


AC_DEFUN([AX_FC_ISO_C_BINDING],[
AC_LANG_PUSH(Fortran)
AC_MSG_CHECKING([if the fortran compiler supports iso c binding])
AC_LINK_IFELSE([
    AC_LANG_SOURCE([
        program conftest
        use ISO_C_BINDING
        integer, dimension(10) :: n
        end
    ])
],[
    FC_ISO_C_BINDINGS="yes"
    AC_MSG_RESULT([$ax_cv_fc_iso_c_bindings])
],[
    FC_ISO_C_BINDINGS="no"
    AC_MSG_RESULT([no, fortran codes and sockets or embedding will not work.])
])
# The intel compiler sometimes generates these work.pc and .pcl files
rm -f work.pc work.pcl
AC_LANG_POP(Fortran)
])

# pkg.m4 - Macros to locate and utilise pkg-config.   -*- Autoconf -*-
# serial 11 (pkg-config-0.29.1)

dnl Copyright © 2004 Scott James Remnant <scott@netsplit.com>.
dnl Copyright © 2012-2015 Dan Nicholson <dbn.lists@gmail.com>
dnl
dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2 of the License, or
dnl (at your option) any later version.
dnl
dnl This program is distributed in the hope that it will be useful, but
dnl WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
dnl General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
dnl 02111-1307, USA.
dnl
dnl As a special exception to the GNU General Public License, if you
dnl distribute this file as part of a program that contains a
dnl configuration script generated by Autoconf, you may include it under
dnl the same distribution terms that you use for the rest of that
dnl program.

dnl PKG_PREREQ(MIN-VERSION)
dnl -----------------------
dnl Since: 0.29
dnl
dnl Verify that the version of the pkg-config macros are at least
dnl MIN-VERSION. Unlike PKG_PROG_PKG_CONFIG, which checks the user's
dnl installed version of pkg-config, this checks the developer's version
dnl of pkg.m4 when generating configure.
dnl
dnl To ensure that this macro is defined, also add:
dnl m4_ifndef([PKG_PREREQ],
dnl     [m4_fatal([must install pkg-config 0.29 or later before running autoconf/autogen])])
dnl
dnl See the "Since" comment for each macro you use to see what version
dnl of the macros you require.
m4_defun([PKG_PREREQ],
[m4_define([PKG_MACROS_VERSION], [0.29.1])
m4_if(m4_version_compare(PKG_MACROS_VERSION, [$1]), -1,
    [m4_fatal([pkg.m4 version $1 or higher is required but ]PKG_MACROS_VERSION[ found])])
])dnl PKG_PREREQ

dnl PKG_PROG_PKG_CONFIG([MIN-VERSION])
dnl ----------------------------------
dnl Since: 0.16
dnl
dnl Search for the pkg-config tool and set the PKG_CONFIG variable to
dnl first found in the path. Checks that the version of pkg-config found
dnl is at least MIN-VERSION. If MIN-VERSION is not specified, 0.9.0 is
dnl used since that's the first version where most current features of
dnl pkg-config existed.
AC_DEFUN([PKG_PROG_PKG_CONFIG],
[m4_pattern_forbid([^_?PKG_[A-Z_]+$])
m4_pattern_allow([^PKG_CONFIG(_(PATH|LIBDIR|SYSROOT_DIR|ALLOW_SYSTEM_(CFLAGS|LIBS)))?$])
m4_pattern_allow([^PKG_CONFIG_(DISABLE_UNINSTALLED|TOP_BUILD_DIR|DEBUG_SPEW)$])
AC_ARG_VAR([PKG_CONFIG], [path to pkg-config utility])
AC_ARG_VAR([PKG_CONFIG_PATH], [directories to add to pkg-config's search path])
AC_ARG_VAR([PKG_CONFIG_LIBDIR], [path overriding pkg-config's built-in search path])

if test "x$ac_cv_env_PKG_CONFIG_set" != "xset"; then
	AC_PATH_TOOL([PKG_CONFIG], [pkg-config])
fi
if test -n "$PKG_CONFIG"; then
	_pkg_min_version=m4_default([$1], [0.9.0])
	AC_MSG_CHECKING([pkg-config is at least version $_pkg_min_version])
	if $PKG_CONFIG --atleast-pkgconfig-version $_pkg_min_version; then
		AC_MSG_RESULT([yes])
	else
		AC_MSG_RESULT([no])
		PKG_CONFIG=""
	fi
fi[]dnl
])dnl PKG_PROG_PKG_CONFIG

dnl PKG_CHECK_EXISTS(MODULES, [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
dnl -------------------------------------------------------------------
dnl Since: 0.18
dnl
dnl Check to see whether a particular set of modules exists. Similar to
dnl PKG_CHECK_MODULES(), but does not set variables or print errors.
dnl
dnl Please remember that m4 expands AC_REQUIRE([PKG_PROG_PKG_CONFIG])
dnl only at the first occurence in configure.ac, so if the first place
dnl it's called might be skipped (such as if it is within an "if", you
dnl have to call PKG_CHECK_EXISTS manually
AC_DEFUN([PKG_CHECK_EXISTS],
[AC_REQUIRE([PKG_PROG_PKG_CONFIG])dnl
if test -n "$PKG_CONFIG" && \
    AC_RUN_LOG([$PKG_CONFIG --exists --print-errors "$1"]); then
  m4_default([$2], [:])
m4_ifvaln([$3], [else
  $3])dnl
fi])

dnl _PKG_CONFIG([VARIABLE], [COMMAND], [MODULES])
dnl ---------------------------------------------
dnl Internal wrapper calling pkg-config via PKG_CONFIG and setting
dnl pkg_failed based on the result.
m4_define([_PKG_CONFIG],
[if test -n "$$1"; then
    pkg_cv_[]$1="$$1"
 elif test -n "$PKG_CONFIG"; then
    PKG_CHECK_EXISTS([$3],
                     [pkg_cv_[]$1=`$PKG_CONFIG --[]$2 "$3" 2>/dev/null`
		      test "x$?" != "x0" && pkg_failed=yes ],
		     [pkg_failed=yes])
 else
    pkg_failed=untried
fi[]dnl
])dnl _PKG_CONFIG

dnl _PKG_SHORT_ERRORS_SUPPORTED
dnl ---------------------------
dnl Internal check to see if pkg-config supports short errors.
AC_DEFUN([_PKG_SHORT_ERRORS_SUPPORTED],
[AC_REQUIRE([PKG_PROG_PKG_CONFIG])
if $PKG_CONFIG --atleast-pkgconfig-version 0.20; then
        _pkg_short_errors_supported=yes
else
        _pkg_short_errors_supported=no
fi[]dnl
])dnl _PKG_SHORT_ERRORS_SUPPORTED


dnl PKG_CHECK_MODULES(VARIABLE-PREFIX, MODULES, [ACTION-IF-FOUND],
dnl   [ACTION-IF-NOT-FOUND])
dnl --------------------------------------------------------------
dnl Since: 0.4.0
dnl
dnl Note that if there is a possibility the first call to
dnl PKG_CHECK_MODULES might not happen, you should be sure to include an
dnl explicit call to PKG_PROG_PKG_CONFIG in your configure.ac
AC_DEFUN([PKG_CHECK_MODULES],
[AC_REQUIRE([PKG_PROG_PKG_CONFIG])dnl
AC_ARG_VAR([$1][_CFLAGS], [C compiler flags for $1, overriding pkg-config])dnl
AC_ARG_VAR([$1][_LIBS], [linker flags for $1, overriding pkg-config])dnl

pkg_failed=no
AC_MSG_CHECKING([for $1])

_PKG_CONFIG([$1][_CFLAGS], [cflags], [$2])
_PKG_CONFIG([$1][_LIBS], [libs], [$2])

m4_define([_PKG_TEXT], [Alternatively, you may set the environment variables $1[]_CFLAGS
and $1[]_LIBS to avoid the need to call pkg-config.
See the pkg-config man page for more details.])

if test $pkg_failed = yes; then
   	AC_MSG_RESULT([no])
        _PKG_SHORT_ERRORS_SUPPORTED
        if test $_pkg_short_errors_supported = yes; then
	        $1[]_PKG_ERRORS=`$PKG_CONFIG --short-errors --print-errors --cflags --libs "$2" 2>&1`
        else 
	        $1[]_PKG_ERRORS=`$PKG_CONFIG --print-errors --cflags --libs "$2" 2>&1`
        fi
	# Put the nasty error message in config.log where it belongs
	echo "$$1[]_PKG_ERRORS" >&AS_MESSAGE_LOG_FD

	m4_default([$4], [AC_MSG_ERROR(
[Package requirements ($2) were not met:

$$1_PKG_ERRORS

Consider adjusting the PKG_CONFIG_PATH environment variable if you
installed software in a non-standard prefix.

_PKG_TEXT])[]dnl
        ])
elif test $pkg_failed = untried; then
     	AC_MSG_RESULT([no])
	m4_default([$4], [AC_MSG_FAILURE(
[The pkg-config script could not be found or is too old.  Make sure it
is in your PATH or set the PKG_CONFIG environment variable to the full
path to pkg-config.

_PKG_TEXT

To get pkg-config, see <http://pkg-config.freedesktop.org/>.])[]dnl
        ])
else
	$1[]_CFLAGS=$pkg_cv_[]$1[]_CFLAGS
	$1[]_LIBS=$pkg_cv_[]$1[]_LIBS
        AC_MSG_RESULT([yes])
	$3
fi[]dnl
])dnl PKG_CHECK_MODULES


dnl PKG_CHECK_MODULES_STATIC(VARIABLE-PREFIX, MODULES, [ACTION-IF-FOUND],
dnl   [ACTION-IF-NOT-FOUND])
dnl ---------------------------------------------------------------------
dnl Since: 0.29
dnl
dnl Checks for existence of MODULES and gathers its build flags with
dnl static libraries enabled. Sets VARIABLE-PREFIX_CFLAGS from --cflags
dnl and VARIABLE-PREFIX_LIBS from --libs.
dnl
dnl Note that if there is a possibility the first call to
dnl PKG_CHECK_MODULES_STATIC might not happen, you should be sure to
dnl include an explicit call to PKG_PROG_PKG_CONFIG in your
dnl configure.ac.
AC_DEFUN([PKG_CHECK_MODULES_STATIC],
[AC_REQUIRE([PKG_PROG_PKG_CONFIG])dnl
_save_PKG_CONFIG=$PKG_CONFIG
PKG_CONFIG="$PKG_CONFIG --static"
PKG_CHECK_MODULES($@)
PKG_CONFIG=$_save_PKG_CONFIG[]dnl
])dnl PKG_CHECK_MODULES_STATIC


dnl PKG_INSTALLDIR([DIRECTORY])
dnl -------------------------
dnl Since: 0.27
dnl
dnl Substitutes the variable pkgconfigdir as the location where a module
dnl should install pkg-config .pc files. By default the directory is
dnl $libdir/pkgconfig, but the default can be changed by passing
dnl DIRECTORY. The user can override through the --with-pkgconfigdir
dnl parameter.
AC_DEFUN([PKG_INSTALLDIR],
[m4_pushdef([pkg_default], [m4_default([$1], ['${libdir}/pkgconfig'])])
m4_pushdef([pkg_description],
    [pkg-config installation directory @<:@]pkg_default[@:>@])
AC_ARG_WITH([pkgconfigdir],
    [AS_HELP_STRING([--with-pkgconfigdir], pkg_description)],,
    [with_pkgconfigdir=]pkg_default)
AC_SUBST([pkgconfigdir], [$with_pkgconfigdir])
m4_popdef([pkg_default])
m4_popdef([pkg_description])
])dnl PKG_INSTALLDIR


dnl PKG_NOARCH_INSTALLDIR([DIRECTORY])
dnl --------------------------------
dnl Since: 0.27
dnl
dnl Substitutes the variable noarch_pkgconfigdir as the location where a
dnl module should install arch-independent pkg-config .pc files. By
dnl default the directory is $datadir/pkgconfig, but the default can be
dnl changed by passing DIRECTORY. The user can override through the
dnl --with-noarch-pkgconfigdir parameter.
AC_DEFUN([PKG_NOARCH_INSTALLDIR],
[m4_pushdef([pkg_default], [m4_default([$1], ['${datadir}/pkgconfig'])])
m4_pushdef([pkg_description],
    [pkg-config arch-independent installation directory @<:@]pkg_default[@:>@])
AC_ARG_WITH([noarch-pkgconfigdir],
    [AS_HELP_STRING([--with-noarch-pkgconfigdir], pkg_description)],,
    [with_noarch_pkgconfigdir=]pkg_default)
AC_SUBST([noarch_pkgconfigdir], [$with_noarch_pkgconfigdir])
m4_popdef([pkg_default])
m4_popdef([pkg_description])
])dnl PKG_NOARCH_INSTALLDIR


dnl PKG_CHECK_VAR(VARIABLE, MODULE, CONFIG-VARIABLE,
dnl [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
dnl -------------------------------------------
dnl Since: 0.28
dnl
dnl Retrieves the value of the pkg-config variable for the given module.
AC_DEFUN([PKG_CHECK_VAR],
[AC_REQUIRE([PKG_PROG_PKG_CONFIG])dnl
AC_ARG_VAR([$1], [value of $3 for $2, overriding pkg-config])dnl

_PKG_CONFIG([$1], [variable="][$3]["], [$2])
AS_VAR_COPY([$1], [pkg_cv_][$1])

AS_VAR_IF([$1], [""], [$5], [$4])dnl
])dnl PKG_CHECK_VAR

dnl PKG_WITH_MODULES(VARIABLE-PREFIX, MODULES,
dnl   [ACTION-IF-FOUND],[ACTION-IF-NOT-FOUND],
dnl   [DESCRIPTION], [DEFAULT])
dnl ------------------------------------------
dnl
dnl Prepare a "--with-" configure option using the lowercase
dnl [VARIABLE-PREFIX] name, merging the behaviour of AC_ARG_WITH and
dnl PKG_CHECK_MODULES in a single macro.
AC_DEFUN([PKG_WITH_MODULES],
[
m4_pushdef([with_arg], m4_tolower([$1]))

m4_pushdef([description],
           [m4_default([$5], [build with ]with_arg[ support])])

m4_pushdef([def_arg], [m4_default([$6], [auto])])
m4_pushdef([def_action_if_found], [AS_TR_SH([with_]with_arg)=yes])
m4_pushdef([def_action_if_not_found], [AS_TR_SH([with_]with_arg)=no])

m4_case(def_arg,
            [yes],[m4_pushdef([with_without], [--without-]with_arg)],
            [m4_pushdef([with_without],[--with-]with_arg)])

AC_ARG_WITH(with_arg,
     AS_HELP_STRING(with_without, description[ @<:@default=]def_arg[@:>@]),,
    [AS_TR_SH([with_]with_arg)=def_arg])

AS_CASE([$AS_TR_SH([with_]with_arg)],
            [yes],[PKG_CHECK_MODULES([$1],[$2],$3,$4)],
            [auto],[PKG_CHECK_MODULES([$1],[$2],
                                        [m4_n([def_action_if_found]) $3],
                                        [m4_n([def_action_if_not_found]) $4])])

m4_popdef([with_arg])
m4_popdef([description])
m4_popdef([def_arg])

])dnl PKG_WITH_MODULES

dnl PKG_HAVE_WITH_MODULES(VARIABLE-PREFIX, MODULES,
dnl   [DESCRIPTION], [DEFAULT])
dnl -----------------------------------------------
dnl
dnl Convenience macro to trigger AM_CONDITIONAL after PKG_WITH_MODULES
dnl check._[VARIABLE-PREFIX] is exported as make variable.
AC_DEFUN([PKG_HAVE_WITH_MODULES],
[
PKG_WITH_MODULES([$1],[$2],,,[$3],[$4])

AM_CONDITIONAL([HAVE_][$1],
               [test "$AS_TR_SH([with_]m4_tolower([$1]))" = "yes"])
])dnl PKG_HAVE_WITH_MODULES

dnl PKG_HAVE_DEFINE_WITH_MODULES(VARIABLE-PREFIX, MODULES,
dnl   [DESCRIPTION], [DEFAULT])
dnl ------------------------------------------------------
dnl
dnl Convenience macro to run AM_CONDITIONAL and AC_DEFINE after
dnl PKG_WITH_MODULES check. HAVE_[VARIABLE-PREFIX] is exported as make
dnl and preprocessor variable.
AC_DEFUN([PKG_HAVE_DEFINE_WITH_MODULES],
[
PKG_HAVE_WITH_MODULES([$1],[$2],[$3],[$4])

AS_IF([test "$AS_TR_SH([with_]m4_tolower([$1]))" = "yes"],
        [AC_DEFINE([HAVE_][$1], 1, [Enable ]m4_tolower([$1])[ support])])
])dnl PKG_HAVE_DEFINE_WITH_MODULES

