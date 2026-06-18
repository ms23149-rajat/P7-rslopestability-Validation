
MODULE_TOPDIR = ../..

SUBDIRS1 = \
	R_code \
	r.slope.stability \
	r.slope.stability.main \
	r.slope.stability.mult

SUBDIRS2 = \
	r.slope.stability \
	r.slope.stability.main \
	r.slope.stability.mult

SUBDIRS = $(SUBDIRS1)

include $(MODULE_TOPDIR)/include/Make/Dir.make

default: parsubdirs

SUBDIRS = $(SUBDIRS2)

install: installsubdirs
