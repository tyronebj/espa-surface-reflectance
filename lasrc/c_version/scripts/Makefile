#-----------------------------------------------------------------------------
# Makefile
#
# Simple makefile for building and installing Landsat SR Code (LaSRC).
#-----------------------------------------------------------------------------
.PHONY: all install clean

# Inherit from upper-level make.config
TOP = ../../..
include $(TOP)/make.config

# Define scripts
SCRIPTS = do_lasrc.py

#-----------------------------------------------------------------------------
all:

install:
	install -d $(link_path)
	install -d $(lasrc_bin_install_path)
	@for script in $(SCRIPTS); do \
		cmd="ln -sf $(lasrc_link_source_path)/$$script $(link_path)/$$script"; \
		echo "Running ... $$cmd"; $$cmd || exit 1; \
		cmd="ln -sf $(lasrc_link_source_path)/$$script $(link_path)/$$script"; \
		echo "Running ... $$cmd"; $$cmd || exit 1; \
	done;

clean:

