
COMPILEFLAGS= -O3

.PRECIOUS:

.c.x:
	$(CC) $(COMPILEFLAGS) $(CFLAGS) -I$(abs_top_srcdir) $? $(LDFLAGS) -Wl,-Bstatic -lflint -Wl,-Bdynamic -lmpfr -lgmp -lcubal -lm -o $@

#$(CC) $(COMPILEFLAGS) $(CFLAGS) -I$(abs_top_srcdir) $? $(LDFLAGS) -Wl,-Bstatic -lflint -lcubal -Wl,-Bdynamic -lm -o $@
