SHELL    =/bin/bash

# run checks for validity of CERN Root installation
have_root=$(shell $(ROOTSYS)/bin/root-config --version | head -1 | awk '{print $$1}')
$(if $(have_root),,$(error CERN Root installation has not been completed. You need to \
        source $$ROOTSYS/bin/thisroot.sh (for bash) or source $$ROOTSYS/bin/thisroot.csh (for C-shell).))
have_rootcint=$(shell rootcint -h 2>&1 | grep -i linkdef | head -1 | awk '{print $$1}')
$(if $(have_rootcint),,$(error command 'rootcint' is either missing or not working. Re-install CERN Root.))

# MAIN DIRECTORY
$(if $(wildcard $(SPECFIT)),,\
        $(info SPECFIT ($(SPECFIT)) was not found, using currect directory))

$(if $(wildcard $(SPECFIT)),,\
        $(eval SPECFIT=$(shell readlink -f $(shell pwd))))

# makefile module that builds the specfit package
include $(SPECFIT)/specfit.mk
