AC_INIT([mbc], [1.1], [pikelner@theor.jinr.ru])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_PROG_CC

AC_SEARCH_LIBS([erf], [m])

dnl
dnl Integrator libraries
dnl
AC_SEARCH_LIBS([llVegas], [cubal],[],[
                          AC_MSG_ERROR([Could  not find CUBA library])
])
AC_SEARCH_LIBS([llCuhre], [cubal],[],[
                          AC_MSG_ERROR([Could  not find CUBA library])
])

dnl
dnl Gamma function with complex argument
dnl

dnl AC_SEARCH_LIBS([acb_gamma], [arb],[mpfr],[
dnl                             AC_MSG_ERROR([Could  not find Arb library])
dnl ])


AC_CONFIG_FILES([Makefile mbc.m])
AC_OUTPUT
