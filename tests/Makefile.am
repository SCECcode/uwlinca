# Autoconf/automake file

#bin_PROGRAMS = test_uwlinca
bin_PROGRAMS = 

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS} -I${UCVM_SRC_PATH}/src/ucvm
AM_LDFLAGS = ${LDFLAGS}

objects = test_api.o
TARGETS = $(bin_PROGRAMS)

all: $(bin_PROGRAMS)

install:
	mkdir -p ${prefix}/tests
#	cp test_uwlinca ${prefix}/tests

test_uwlinca$(EXEEXT): $(objects)
	$(CC) -o $@ $^ $(AM_CFLAGS) -L../src -luwlinca $(AM_LDFLAGS)

$(objects): %.o: %.c
	$(CC) -o $@ -c $^ $(AM_CFLAGS) -I../src/ $(AM_CFLAGS)
