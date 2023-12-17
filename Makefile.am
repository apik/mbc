
COMPILEFLAGS= -O3

.PRECIOUS:

.c.x:
	$(CC) $(COMPILEFLAGS) $(CFLAGS) -I$(abs_top_srcdir) $? $(LDFLAGS) -lcuba -lm -o $@
