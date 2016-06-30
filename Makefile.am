
COMPILEFLAGS= -O3

.c.x:
	$(CC) $(COMPILEFLAGS) $(CFLAGS) -I$(abs_top_srcdir) $? $(LDFLAGS) -lcuba -lm -o $@
