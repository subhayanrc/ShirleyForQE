# Makefile for PW
sinclude ../make.inc

default: all

all: shirley

shirley:
	( cd src ; $(MAKE) all || exit 1 )

clean : 
	( cd src ; $(MAKE) clean )

distclean: clean 
