OBJS=main.o

genpad: $(OBJS)
	gcc -s -o genpad $(OBJS) -lm -lpopt -lfftw3 -lsndfile

$(OBJS): %.o: %.c
	gcc -O2 -c -o $@ $<
