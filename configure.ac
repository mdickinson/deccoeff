AC_PREREQ(2.61)
AC_INIT(deccoeff, 0.2)

AC_CONFIG_HEADERS(deccoeff_config.h)

dnl checks for header files

AC_HEADER_STDBOOL
AC_CHECK_HEADERS(stdint.h)

dnl ideally, we'd like to use uint32_t for our Deccoeff limbs, which
dnl means that we also need uint64_t to hold intermediate results when
dnl doing multiplications and divisions.  But these types aren't part
dnl of C89, so we need to check for them (and provide fallbacks if
dnl they don't exist).

AC_TYPE_UINT32_T
AC_TYPE_UINT64_T

AC_OUTPUT

