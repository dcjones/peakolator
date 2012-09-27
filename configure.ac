AC_INIT([peakolator],
        m4_esyscmd([./git-version-gen | tr -d '\n']),
        [dcjones@cs.washington.edu])
AM_INIT_AUTOMAKE([foreign -Wall])
AM_SILENT_RULES([yes])
AC_CONFIG_MACRO_DIR([m4])

AM_PROG_AR
LT_INIT
AC_PROG_CC

opt_CFLAGS="-std=gnu99 -Wall -Wextra -pedantic -g -O3 -DNDEBUG"
dbg_CFLAGS="-std=gnu99 -Wall -Wextra -pedantic -g -O0"

AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],
                              [enable debugging info (default is no)])],
              [], [enable_debug=no])

AS_IF([test "x$enable_debug" = xyes],
      [CFLAGS="$dbg_CFLAGS"],
      [CFLAGS="$opt_CFLAGS"])

ACX_PTHREAD()
LIBS="$PTHREAD_LIBS $LIBS"
CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
CC="$PTHREAD_CC"

AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile examples/Makefile])
AC_OUTPUT
