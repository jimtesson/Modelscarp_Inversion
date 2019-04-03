AC_PREREQ([2.67])
AC_INIT([Modelscarp Inversion], [1.0.0], [https://github.com/jimtesson/Modelscarp_Inversion])
AM_INIT_AUTOMAKE([1.10 no-define subdir-objects])
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_F77
AC_PROG_FC
AC_PROG_INSTALL

PKG_CHECK_MODULES([RJMCMC], [rjmcmc])

AC_SUBST([RJMCMC_CFLAGS])
AC_SUBST([RJMCMC_LDFLAGS])

# Checks for header files.


# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([
	Makefile
	])
AC_OUTPUT
