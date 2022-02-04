# makefile.
# customized for TH/pm3
# change all DRIVER instances to whatever it needs to be.
whoisthis=${USER}
#
#
CC       	= mpicc
CFLAGS   	= -g -O2
MY_APP	 	= cbf

ifeq (${USER}, pm3user)
	INC_DIRS  	+= /home/pm3user/software/sundials/instdir/include
endif
# Add a prefix to INC_DIRS. So moduleA would become -ImoduleA. GCC understands this -I flag
INC_FLAGS := $(addprefix -I,$(INC_DIRS))


ifeq (${USER}, pm3user)
	LDFLAGS	+= -L/home/pm3user/software/sundials/instdir/lib
endif
LDFLAGS += -lm
LDFLAGS += -lsundials_cvodes
LDFLAGS += -lsundials_nvecserial

cbf:	ursino.c
	$(CC) $(CFLAGS) $(INC_FLAGS) -c ursino.c -o ursino.o
	$(CC) $(CFLAGS) ursino.o $(LDFLAGS) -o $(MY_APP)

run:
	./${MY_APP}

clean:
	rm  ${MY_APP}

veryclean:
	rm -r *.dat ${MY_APP} ${MY_APPI} *.o *~ *.txt *.out
