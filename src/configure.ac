#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

m4_include([./ensc_annotation_tools_version.m4])
m4_include([./ax_lib_mysql.m4])

AC_PREREQ([2.61])
AC_INIT([EnsC-annotation-tools], [VERSION_NUMBER], [dev@ensembl.org])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_SRCDIR([OBDA/BioIndex.h])
AM_INIT_AUTOMAKE([1.9 foreign])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([limits.h malloc.h stdlib.h string.h strings.h sys/param.h unistd.h])

AC_CONFIG_FILES([Makefile OBDA/Makefile Programs/Makefile])

# Check if this is a Git repository and, if so, include info about the current
# status in the version info. Otherwise use PACKAGE_VERSION (which will get
# set from VERSION_NUMBER in version.m4).
if test -d "$srcdir/../.git"; then
  AC_DEFINE_UNQUOTED([ENSC_ANNOTATION_TOOLS_VERSION], [git describe --always --abbrev=1], [ensC-annotation-tools package version])
else
  AC_DEFINE([ENSC_ANNOTATION_TOOLS_VERSION], [PACKAGE_VERSION], [ensC-annotation-tools package version])
fi

AC_CONFIG_MACRO_DIR([m4])

LT_INIT([disable-shared])

LT_LANG([C])

# keep libtool up to date automatically.
AC_SUBST([LIBTOOL_DEPS])


# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_INT8_T
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([bzero getcwd memmove memset mkdir realpath stpcpy strcasecmp strchr strpbrk strrchr strstr strtol])

# Fail if EnsC-core is not present
AC_CHECK_HEADER([ensc/EnsC.h], [AC_MSG_NOTICE([EnsC-core is available])], [AC_MSG_ERROR([EnsC is not available you can find it in github at https://github.com/Ensembl/ensc-core])])

# Standard macro to test for mysql
AX_LIB_MYSQL()

AM_CONDITIONAL([HAVE_MYSQL_LIB], [test -n "$MYSQL_VERSION"])
AM_COND_IF([HAVE_MYSQL_LIB], [AC_SUBST([HAVE_MYSQL_LIB], [1])], [AC_MSG_WARN([MySQL is not present])])
AM_COND_IF([HAVE_MYSQL_LIB], [[ CFLAGS="$CFLAGS $MYSQL_CFLAGS" ]])

# Also need to check for headers
AC_CHECK_HEADER([mysql.h], [HAVE_MYSQL_HEADER=1])
AM_CONDITIONAL([HAVE_MYSQL_HEADER], [ test -n "$HAVE_MYSQL_HEADER" && test -n "$HAVE_MYSQL_LIB" ])
AM_COND_IF([HAVE_MYSQL_HEADER], [AC_DEFINE([MYSQL_STRICT], [1], [define to 1 if mysql is in strict mode])])

AC_CHECK_HEADERS([htslib/sam.h htslib/hts.h], [HAVE_SAMTOOLS_HEADER=1])
AM_CONDITIONAL([HAVE_SAMTOOLS], [ test -n "$HAVE_SAMTOOLS_HEADER" ])
AM_COND_IF([HAVE_SAMTOOLS], [AC_MSG_NOTICE([htslib is available])])
AM_COND_IF([HAVE_SAMTOOLS], [AC_DEFINE([HAVE_SAMTOOLS], [1], [define to 1 if htslib is available])])

AC_CHECK_HEADERS([libconfig.h], [HAVE_LIBCONFIG_HEADER=1])
AM_CONDITIONAL([HAVE_LIBCONFIG], [ test -n "$HAVE_LIBCONFIG_HEADER" ])
AM_COND_IF([HAVE_LIBCONFIG], [AC_MSG_NOTICE([libconfig is available])])
AM_COND_IF([HAVE_LIBCONFIG], [AC_DEFINE([HAVE_LIBCONFIG], [1], [define to 1 if libconfig is available])])

AC_CHECK_HEADERS([gperftools/tcmalloc.h], [HAVE_LIBTCMALLOC_HEADER=1])
AM_CONDITIONAL([HAVE_LIBTCMALLOC], [ test -n "$HAVE_LIBTCMALLOC_HEADER" ])
AM_COND_IF([HAVE_LIBTCMALLOC], [AC_MSG_NOTICE([libtcmalloc is available])])
AM_COND_IF([HAVE_LIBTCMALLOC], [AC_DEFINE([HAVE_LIBTCMALLOC], [1], [define to 1 if libtcmalloc is available])])



AC_SUBST(PKG_CONFIG)
AC_SUBST(CFLAGS)
AC_SUBST(CPPFLAGS)
AC_SUBST(LDFLAGS)


AC_OUTPUT


echo \
"-------------------------------------------------

 ${PACKAGE_NAME} Version ${PACKAGE_VERSION}

 Prefix: '${prefix}'
 C++ Compiler: '${CXX} ${CXXFLAGS} ${CPPFLAGS}'
   C Compiler: '${CC} ${CFLAGS} ${CPPFLAGS}'

 Now type 'make @<:@<target>@:>@'
   where the optional @<:@<target>@:>@ is:
     all       - build all libraries (default)
     install   - build and install everything
"

if test -z "$HAVE_MYSQL_HEADER"
then
  echo " Warning: MySql is not installed, only indicate fastasplit_random translate and countsnps will be compiled"
  echo
fi

if test -z "$HAVE_SAMTOOLS_HEADER"
then
  echo " samtools not available: bam* programs and RefineSolexaGenes will not be compiled"; echo
fi

if test -z "$HAVE_LIBCONFIG_HEADER"
then
  echo " libconfig not available: RefineSolexaGenes will not be compiled"; echo
fi

if test -z "$HAVE_LIBTCMALLOC_HEADER"
then
  echo " libtcmalloc not available: RefineSolexaGenes will be compiled using malloc"; echo
fi


echo "---------------------------------------------------------------------
"
