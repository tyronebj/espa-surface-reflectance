#-----------------------------------------------------------------------------
# Makefile
#
# Simple makefile for building and installing L8 SR.
#-----------------------------------------------------------------------------
.PHONY: all install clean

all:
	echo "make all in src..."; \
        (cd src; $(MAKE))
	echo "make all in src_pre_collection..."; \
        (cd src_pre_collection; $(MAKE))

install: all
	echo "make install in scripts..."; \
        (cd scripts; $(MAKE) install)
	echo "make install in src..."; \
        (cd src; $(MAKE) install)
	echo "make install in src_pre_collection..."; \
        (cd src_pre_collection; $(MAKE) install)

clean:
	echo "make clean in src..."; \
        (cd src; $(MAKE) clean)
	echo "make clean in src_pre_collection..."; \
        (cd src_pre_collection; $(MAKE) clean)

