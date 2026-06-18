MODULE_TOPDIR = ../..

include $(MODULE_TOPDIR)/include/Make/Other.make

DSTDIR = $(HOME)/.grass7/addons/etc/r.slope.stability.rcode

# R
R_SRC := $(wildcard *.R)
R_DST := $(patsubst %,$(DSTDIR)/%,$(R_SRC))

default: $(R_DST)

$(DSTDIR)/%.R: %.R | $(DSTDIR)
	$(INSTALL_DATA) $< $@

$(DSTDIR):
	$(MKDIR) $@
