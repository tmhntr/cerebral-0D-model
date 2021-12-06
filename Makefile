# makefile.
# customized for TH/pm3
# change all DRIVER instances to whatever it needs to be.
whoisthis=$USER
#
#
CC       	= mpicc
CFLAGS   	= -g -O2
INCLUDE  	= /home/${whoisthis}/software/sundials/instdir/include
MY_APP	 	= cbf
LIB	 			= -L/home/${whoisthis}/software/sundials/instdir/lib

cbf:	ursino.c
	${CC} ${CFLAGS} -I${INCLUDE} -c ursino.c -o ursino.o
	${CC} ${CFLAGS} ursino.o -I${INCLUDE} -lm ${LIB} -lsundials_cvodes -lsundials_nvecserial -o ${MY_APP}

run:
	./${MY_APP}

clean:
	rm  ${MY_APP}

veryclean:
	rm -r *.dat ${MY_APP} ${MY_APPI} *.o *~ *.txt *.out
