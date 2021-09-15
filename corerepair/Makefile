# Makefile for corerepair utilities

include ../make.inc


# location of needed modules and included files (if any)
MODFLAGS = $(BASEMOD_FLAGS) \
           $(MOD_FLAG)$(TOPDIR)/upftools

PWOBJS = ../PW/src/libpw.a
QEOBJS = ../Modules/libqemod.a ../KS_Solvers/libks_solvers.a \
         ../FFTXlib/libqefft.a ../LAXlib/libqela.a ../UtilXlib/libutil.a \
         ../dft-d3/libdftd3qe.a

MODULES = $(QEOBJS) $(MODS) ../upftools/write_upf.o

# local files
MODS = \
splines.o \
pointers.o \
fileio.o \
atomic_waves.o \
wigner3j.o 

ADDLOCALPROJ_OBJ= \
addlocalproj_upf.o \
intradialprod.o 

VALENCECOREPAIR_OBJ= \
valence_corerepair.o \
renormalize_ae_from_ps.o \
momentum_matrix.o

COREVALENCE_POSITION_OBJ= \
corevalence_position.o \
position_matrix.o

VALENCE_OVERLAP_OBJ= \
valence_overlap.o 

#TLDEPS= pwlibs
TLDEPS = mods libs libfox upf

all : tldeps addlocalproj_upf.x valence_corerepair.x corevalence_position.x valence_overlap.x

addlocalproj_upf.x : $(ADDLOCALPROJ_OBJ) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		$(ADDLOCALPROJ_OBJ) $(MODULES) $(LIBOBJS) $(QELIBS)

valence_corerepair.x : $(VALENCECOREPAIR_OBJ) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		$(VALENCECOREPAIR_OBJ) $(MODULES) $(LIBOBJS) $(QELIBS)

corevalence_position.x : $(COREVALENCE_POSITION_OBJ) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		$(COREVALENCE_POSITION_OBJ) $(MODULES) $(LIBOBJS) $(QELIBS)

valence_overlap.x : $(VALENCE_OVERLAP_OBJ) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		$(VALENCE_OVERLAP_OBJ) $(MODULES) $(LIBOBJS) $(QELIBS)

tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f  *.x *.o *~ *.F90 *.mod *.d *.i work.pc

include make.depend
