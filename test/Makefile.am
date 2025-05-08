# Autoconf/automake file

bin_PROGRAMS = test_uwlinca

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS} ${PROJ_INCLUDES} -I${UCVM_SRC_PATH}/src/ucvm -I../src
AM_LDFLAGS = ${LDFLAGS} ${PROJ_LDFLAGS}  -L../src -luwlinca -lm

objects = test_api.o
TARGETS = $(bin_PROGRAMS)

.PHONY = run_unit

all: $(bin_PROGRAMS)

############################################
# Executables
############################################

install: all
	mkdir -p ${prefix}/test
	cp test_uwlinca ${prefix}/test

test_uwlinca: $(objects)
	$(CC) -o $@ $^ $(AM_LDFLAGS) 

run_unit : test_uwlinca
	./test_uwlinca

clean: 
	rm -rf test_uwlinca test_api.o
