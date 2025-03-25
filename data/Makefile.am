# GNU Automake config


# Autoconf/automake file

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS}
AM_FCFLAGS = ${FCFLAGS}
AM_LDFLAGS = ${LDFLAGS}

install:
if WITH_UWLINCA_LARGEDATA_DIR
	$(MKDIR_P) ${prefix}/model
 	$(LN_S) -f $(UWLINCA_LARGEDATA_DIR) ${prefix}/model
	echo "USING WITH_UWLINCA_LARGEDATA_DIR in uwlinca"
else
	mkdir -p ${prefix}/data
	cp ./config ${prefix}/data
	cp -R ./uwlinca ${prefix}/data
	chmod og+r ${prefix}/data/config
	chmod -R og+r ${prefix}/data/uwlinca
	chmod og+x ${prefix}/data/uwlinca
endif

clean:
	echo "clean"

