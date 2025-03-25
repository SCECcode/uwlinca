# Autoconf/automake file

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS} -I${UCVM_SRC_PATH}/src/ucvm
AM_LDFLAGS = ${LDFLAGS}

TARGETS = libuwlinca.a libuwlinca.so

all: $(TARGETS)

install:
	mkdir -p ${prefix}
	mkdir -p ${prefix}/lib
	mkdir -p ${prefix}/include
	cp libuwlinca.so ${prefix}/lib
	cp libuwlinca.a ${prefix}/lib
	cp uwlinca.h ${prefix}/include

libuwlinca.a: uwlinca_static.o 
	$(AR) rcs $@ $^

libuwlinca.so: uwlinca.o
	$(CC) -shared $(AM_CFLAGS) -o libuwlinca.so $^ $(AM_LDFLAGS)

uwlinca.o: uwlinca.c
	$(CC) -fPIC -DDYNAMIC_LIBRARY -o $@ -c $^ $(AM_CFLAGS)
	
uwlinca_static.o: uwlinca.c
	$(CC) -o $@ -c $^ $(AM_CFLAGS)
	
clean:
	rm -rf $(TARGETS)
	rm -rf *.o 

