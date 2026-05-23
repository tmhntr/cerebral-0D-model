# Makefile for cerebral-0D-model
# Requires SUNDIALS (brew install sundials)

CC       = cc
CFLAGS   = -g -O2
SUNDIALS = $(shell brew --prefix sundials)
OPENMPI  = $(shell brew --prefix open-mpi)
INCLUDE  = $(SUNDIALS)/include
MY_APP   = cbf
LIB      = -L$(SUNDIALS)/lib -L$(OPENMPI)/lib

LIBS     = -lsundials_cvodes -lsundials_nvecserial -lsundials_sunlinsoldense \
           -lsundials_sunmatrixdense -lsundials_sunnonlinsolnewton \
           -lsundials_core -lmpi -lm

cbf: ursino.c ursino.h f_ursino.c cerebral.c baroreflex.c pumpingfxn.c p_ursino.c
	${CC} ${CFLAGS} -I${INCLUDE} -c ursino.c -o ursino.o
	${CC} ${CFLAGS} ursino.o ${LIB} ${LIBS} -o ${MY_APP}

run:
	./${MY_APP}

clean:
	rm -f ${MY_APP} *.o

veryclean:
	rm -f *.dat ${MY_APP} *.o *~ *.txt *.out
