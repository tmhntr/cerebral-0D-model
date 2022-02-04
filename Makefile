# makefile.
# customized for TH/pm3
# change all DRIVER instances to whatever it needs to be.
whoisthis=${USER}
#
#
CC       	= mpicc
CFLAGS   	= -g -O2 -Wall
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

cbf:	CBF_parameters.o CBF_RHS.o CBF_driver.o
	$(CC) $(CFLAGS) $(INC_FLAGS) $^ $(LDFLAGS) -o $@

clean:
	rm  ${MY_APP}

veryclean:
	rm -r *.dat ${MY_APP} ${MY_APPI} *.o *~ *.txt *.out
