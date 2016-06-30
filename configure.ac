AC_INIT([mbpp], [1.0], [pikelner@theor.jinr.ru])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_PROG_CC

dnl AC_SEARCH_LIBS([Vegas], [cuba])

AC_SEARCH_LIBS([erf], [m])
AC_SEARCH_LIBS([Vegas], [cuba],[],
                        [
                          AC_MSG_ERROR([Could  not find CUBA library])
                        ])
AC_SEARCH_LIBS([Cuhre], [cuba],[],
                        [
                          AC_MSG_ERROR([Could  not find CUBA library])
                        ])

AC_SUBST([AM_CFLAGS], [-O3 -fpic -Wall -Wextra])

dnl AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
