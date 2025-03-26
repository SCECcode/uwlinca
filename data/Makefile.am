# GNU Automake config


# Autoconf/automake file

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS}
AM_FCFLAGS = ${FCFLAGS}
AM_LDFLAGS = ${LDFLAGS}

all:
if WITH_UWLINCA_LARGEDATA_DIR
	echo "USING WITH_UWLINCA_LARGEDATA_DIR in uwlinca"
else
	if [ ! -d "uwlinca" ]; then \
		./make_data_files.py; \
	fi
endif


install:
if WITH_UWLINCA_LARGEDATA_DIR
	$(MKDIR_P) ${prefix}/model
	$(LN_S) -f $(UWLINCA_LARGEDATA_DIR) ${prefix}/model
	echo "USING WITH_UWLINCA_LARGEDATA_DIR in uwlinca"
	cp config ${prefix}/data
else
	if [ ! -d "${prefix}/data/uwlinca" ]; then \
		$(MKDIR_P) ${prefix}/data; \
		mv uwlinca ${prefix}/data; \
		$(LN_S) -f ${prefix}/data/uwlinca .; \
	fi
        cp config ${prefix}/data
endif

clean:
	echo "clean"

