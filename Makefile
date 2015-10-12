CFLAGS = -Wall -g -I/usr/X11R6/include

.SUFFIX: .c .o

SRCS = pdf417decode.c \
	pdf417_dham.c \
	pdf417rs.c

OBJS = $(SRCS:.c=.o)

.c.o:
	gcc $(CFLAGS) -c $<

all: pdf417decode

pdf417decode: $(OBJS)
	gcc -g -o $@ $(OBJS) -L/usr/X11R6/lib -lpbm

clean:
	-rm -f *.o *~ pdf417decode

check:
	@for i in test/*.pbm*; do \
		echo Decoding image $$i ; \
		./pdf417decode -rs -e $$i ; \
		echo ; \
	done
