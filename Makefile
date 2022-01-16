# makefile.
# customized for TH/pm3


# Thanks to Job Vranish (https://spin.atomicobject.com/2016/08/26/makefile-c-projects/)
TARGET_EXEC := cbf

BUILD_DIR := build
SRC_DIRS := src

# Find all the C and C++ files we want to compile
# Note the single quotes around the * expressions. Make will incorrectly expand these otherwise.
SRCS := $(shell find $(SRC_DIRS) -name '*.cpp' -or -name '*.c' -or -name '*.s')

# String substitution for every C/C++ file.
# As an example, hello.cpp turns into ./build/hello.cpp.o
OBJS := $(SRCS:$(SRC_DIRS)/%=$(BUILD_DIR)/%.o)

# String substitution (suffix version without %).
# As an example, ./build/hello.cpp.o turns into ./build/hello.cpp.d
DEPS := $(OBJS:.o=.d)

# Every folder in ./src will need to be passed to GCC so that it can find header files
INC_DIRS := $(shell find $(SRC_DIRS) -type d)
ifeq (${USER}, pm3user)
	INC_DIRS  	+= /home/pm3user/software/sundials/instdir/include
else
	INC_DIRS 	+= ${PWD}/extern/sundials/instdir/include
endif
# Add a prefix to INC_DIRS. So moduleA would become -ImoduleA. GCC understands this -I flag
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

# The -MMD and -MP flags together generate Makefiles for us!
# These files will have .d instead of .o as the output.
CPPFLAGS := $(INC_FLAGS) -MMD -MP

CC       	= mpicc
CFLAGS   	= -g -O2
MY_APP	 	= cbf


ifeq (${USER}, pm3user)
	LDFLAGS	 		+= -L/home/pm3user/software/sundials/instdir/lib
else
	LDFLAGS 		+= -L${PWD}/extern/sundials/instdir/lib
endif

LDFLAGS += -lm
LDFLAGS += -lsundials_cvodes
LDFLAGS += -lsundials_nvecserial

# The final build step.
$(BUILD_DIR)/$(TARGET_EXEC): $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDFLAGS)

# Build step for C source
$(BUILD_DIR)/%.c.o: $(SRC_DIRS)/%.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -r $(BUILD_DIR)
