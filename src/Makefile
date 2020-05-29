# Makefile for ACOThOP

OPTIM_FLAGS=-O3
WARN_FLAGS=-w -ansi -pedantic
CFLAGS=$(WARN_FLAGS) $(OPTIM_FLAGS)
LDLIBS=-lm

acothop: acothop.o thop.o utilities.o ants.o inout.o unix_timer.o ls.o parse.o

all: clean acothop

clean:
	@$(RM) *.o acothop

acothop: acothop.o thop.o utilities.o ants.o inout.o unix_timer.o ls.o parse.o

acothop.o: acothop.c

thop.o: thop.c thop.h

ants.o: ants.c ants.h

inout.o: inout.c inout.h

utilities.o: utilities.c utilities.h

ls.o: ls.c ls.h

parse.o: parse.c parse.h

unix_timer.o: unix_timer.c timer.h

