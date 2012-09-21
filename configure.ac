AC_INIT([peakolator],
        m4_esyscmd([./git-version-gen | tr -d '\n']),
        [dcjones@cs.washington.edu])
AM_INIT_AUTOMAKE([foreign -Wall])
AM_SILENT_RULES([yes])
AC_CONFIG_MACRO_DIR([m4])

AM_PROG_AR
LT_INIT
AC_PROG_CC

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
