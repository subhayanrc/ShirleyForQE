# Makefile for postprocessing (PP)

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS = $(BASEMOD_FLAGS) \
           $(MOD_FLAG)../../PW/src \
           $(MOD_FLAG)../../FFTXlib

QEMODS=../../PW/src/libpw.a ../../Modules/libqemod.a \
       ../../KS_Solvers/libks_solvers.a \
       ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a \
       ../../UtilXlib/libutil.a ../../dft-d3/libdftd3qe.a

#MODFLAGS= $(MOD_FLAG)../iotk/src $(MOD_FLAG)../Modules $(MOD_FLAG)../PW $(MOD_FLAG).

# Important vaiables 
# DEBUG : for use in modifying the form of executable names EXE$(DEBUG).x
# COPY  : changes how the exe is transferred to the ../bin directory

DEBUG=
#DEBUG=_dbg
COPY=ln -fs
#COPY=cp --remove-destination

SHIRLEY_BASIS = \
shirley_basis.o \
shirley_basis_input.o \
start_shirley.o \
stop_shirley.o \
append_missing_slash.o \
print_clock_shirley.o \
build_optimal_basis.o \
scalapack_module.o \
wfc_shirley.o \
basis_shirley.o \
pot_rho_G.o \
deallocate_fft.o \
freeunit.o \
gk_l2gmap.o \
kpoint_expand.o \
rgrid.o 

SHIRLEY_HAM = \
shirley_ham.o \
shirley_ham_input.o \
start_shirley.o \
stop_shirley.o \
read_file_shirley.o \
append_missing_slash.o \
hamq.o \
print_clock_shirley.o \
mp_scatt.o \
hamq_shirley.o \
hamq_pool.o \
splines.o \
bspline90_22.o \
mpio.o \
freeunit.o \
init_us_2_shirley.o \
scalapack_module.o \
orthoatwfc_shirley.o \
atomic_wfc_shirley.o \
update_pseudo.o \
reallocate_nlpot.o 

SHIRLEY_QDIAGP = shirley_qdiagp.o $(NEWDIAGOBJ)

NEWDIAGOBJ= \
bspline90_22.o \
corerepair.o \
diag.o \
freeunit.o \
hamq_pool.o \
hamq_shirley.o \
kpt.o \
mpio.o \
mp_scatt.o \
plot.o \
scalapack_module.o \
shirley_constants.o \
shirley_input.o \
start_shirley.o

BANDSTRUCT = bandstruct.o start_shirley.o freeunit.o mpio.o

SHIRLEY_XAS = shirley_xas.o $(NEWDIAGOBJ)
XAS_PARA = xas_para.o start_shirley.o freeunit.o mpio.o
EFERMI = efermi.o fermi.o \
$(NEWDIAGOBJ)
#
SHIRLEY_OVERLAP = \
append_missing_slash.o \
print_clock_shirley.o \
shirley_overlap.o \
read_file_shirley.o \
reallocate_nlpot.o \
shirley_overlap_input.o \
start_shirley.o \
stop_shirley.o 
#
#DIAGOBJ= \
#bspline90_22.o \
#freeunit.o \
#hamq_shirley.o \
#start_shirley.o
#
#SHIRLEY_PLOT = shirley_plot.o $(NEWDIAGOBJ)
#SHIRLEY_EFERMI = shirley_efermi.o  fermi.o ../PW/efermit.o ../PW/sumkt.o $(NEWDIAGOBJ)
#SHIRLEY_EPSILON = shirley_epsilon.o scalapack_module.o fermi.o $(NEWDIAGOBJ)
#SHIRLEY_EPSTCI = shirley_epstci.o scalapack_module.o fermi.o stop_shirley.o print_clock_shirley.o $(NEWDIAGOBJ)
#SHIRLEY_EPSILON_TRIPLES = shirley_epsilon_triples.o fermi.o stop_shirley.o print_clock_shirley.o splines.o $(NEWDIAGOBJ)
#SHIRLEY_EPSILON00 = shirley_epsilon00.o fermi.o stop_shirley.o print_clock_shirley.o splines.o $(NEWDIAGOBJ)
#SHIRLEY_EPSILONIJ = shirley_epsilonij.o fermi.o stop_shirley.o print_clock_shirley.o splines.o $(NEWDIAGOBJ)
#SHIRLEY_EPSILONI0 = shirley_epsiloni0.o fermi.o stop_shirley.o print_clock_shirley.o splines.o $(NEWDIAGOBJ)
#SHIRLEY_EPSILONQ = shirley_epsilonq.o paulipairs.o fermi.o stop_shirley.o print_clock_shirley.o splines.o $(NEWDIAGOBJ)
#TRIPLES_TEST = triples_test.o fermi.o stop_shirley.o print_clock_shirley.o splines.o $(NEWDIAGOBJ)
#SHIRLEY_SELFENERGY_TRIPLES = shirley_selfenergy_triples.o fermi.o stop_shirley.o print_clock_shirley.o splines.o shirley_output.o $(NEWDIAGOBJ)
#SHIRLEY_EPSILON_CONSTANT = shirley_epsilon_constant.o stop_shirley.o print_clock_shirley.o splines.o $(NEWDIAGOBJ)
#
#
##SHIRLEY_QDIAG = shirley_qdiag.o $(DIAGOBJ)
#
#SHIRLEY_XASRESP = shirley_xasresp.o gmres.o zPackgmres.o $(NEWDIAGOBJ)
#SHIRLEY_EIG = shirley_eig.o $(NEWDIAGOBJ)
#SHIRLEY_PROJ = shirley_proj.o splines.o $(NEWDIAGOBJ)
#SHIRLEY_DHDK = shirley_dhdk.o mpio.o \
#               ../PW/tweights.o ../PW/efermit.o ../PW/sumkt.o ../flib/sort.o \
#               $(NEWDIAGOBJ)
#SHIRLEY_MCSAMP = shirley_mcsamp.o fermi.o $(NEWDIAGOBJ)
#SHIRLEY_MOM = shirley_mom.o kpt.o corerepair.o \
#              $(DIAGOBJ)
#GETTEST = gettest.o start_shirley.o
#DOS_PARA = dos_para.o start_shirley.o freeunit.o
#PDOS_PARA = pdos_para.o start_shirley.o proj_shirley.o freeunit.o
#JDOS_PARA = jdos_para.o start_shirley.o freeunit.o
#XARS_PARA = xars_para.o start_shirley.o freeunit.o mpio.o
#COND_PARA = cond_para.o fermi.o start_shirley.o freeunit.o mpio.o \
#            ../PW/efermit.o ../PW/sumkt.o ../flib/dost.o ../flib/sort.o
#PLASMA_FREQ = plasma_freq.o fermi.o start_shirley.o freeunit.o mpio.o
#
#SHIRLEY_OPT = shirley_opt.o $(NEWDIAGOBJ)
#
#SHIRLEY_QMOMP = shirley_qmomp.o print_clock_shirley.o diag_shirley.o $(DIAGOBJ)
#
#SHIRLEY_LIFETIME = shirley_lifetime.o \
#                   diag_shirley.o elph_shirley.o print_clock_shirley.o \
#                   $(DIAGOBJ)
#
#SHIRLEY_XLINEWIDTH = shirley_xlinewidth.o \
#                     diag_shirley.o elph_shirley.o print_clock_shirley.o \
#                     $(DIAGOBJ)
#
#SHIRLEY_JDOS = shirley_jdos.o \
#               diag_shirley.o elph_shirley.o print_clock_shirley.o \
#               $(DIAGOBJ)
#
#PWMAT_TRANSPOSE = \
#freeunit.o \
#pwmat_transpose.o
#
#DIAG_EPSILON = \
#start_shirley.o \
#stop_shirley.o \
#
#CHANGE_RHOGRID = \
#append_missing_slash.o \
#change_rhogrid.o \
#deallocate_fft.o \
#freeunit.o \
#kpoint_expand.o \
#pot_rho_G.o \
#print_clock_shirley.o \
#read_file_shirley.o \
#rgrid.o \
#start_shirley.o \
#stop_shirley.o \
#wfc_shirley.o 
#
#SHIRLEY_STIELTJES = shirley_stieltjes.o $(NEWDIAGOBJ)
#
#WRITETEST = \
#bspline90_22.o \
#hamq_pool.o \
#hamq_shirley.o \
#mpio.o \
#print_clock_shirley.o \
#scalapack_module.o \
#start_shirley.o \
#stop_shirley.o \
#writetest.o

SHIRLEY_BASIS_EXE=shirley_basis$(DEBUG).x
SHIRLEY_HAM_EXE=shirley_ham$(DEBUG).x
SHIRLEY_QDIAGP_EXE=shirley_qdiagp$(DEBUG).x
BANDSTRUCT_EXE=bandstruct$(DEBUG).x
SHIRLEY_XAS_EXE=shirley_xas$(DEBUG).x
XAS_PARA_EXE=xas_para$(DEBUG).x
EFERMI_EXE=efermi$(DEBUG).x
SHIRLEY_OVERLAP_EXE=shirley_overlap$(DEBUG).x

#SHIRLEY_OVERLAP_EXE=shirley_overlap$(DEBUG).x
##SHIRLEY_QDIAG_EXE=shirley_qdiag$(DEBUG).x
#SHIRLEY_PLOT_EXE=shirley_plot$(DEBUG).x
#SHIRLEY_XASRESP_EXE=shirley_xasresp$(DEBUG).x
#SHIRLEY_EIG_EXE=shirley_eig$(DEBUG).x
#SHIRLEY_PROJ_EXE=shirley_proj$(DEBUG).x
#SHIRLEY_DHDK_EXE=shirley_dhdk$(DEBUG).x
#SHIRLEY_MCSAMP_EXE=shirley_mcsamp$(DEBUG).x
#SHIRLEY_MOM_EXE=shirley_mom$(DEBUG).x
#DOS_PARA_EXE=dos_para$(DEBUG).x
#PDOS_PARA_EXE=pdos_para$(DEBUG).x
#JDOS_PARA_EXE=jdos_para$(DEBUG).x
#XARS_PARA_EXE=xars_para$(DEBUG).x
#COND_PARA_EXE=cond_para$(DEBUG).x
#PLASMA_FREQ_EXE=plasma_freq$(DEBUG).x
#GETTEST_EXE=gettest$(DEBUG).x
#SHIRLEY_OPT_EXE=shirley_opt$(DEBUG).x
#SHIRLEY_QMOMP_EXE=shirley_qmomp$(DEBUG).x
#SHIRLEY_EFERMI_EXE=shirley_efermi$(DEBUG).x
#SHIRLEY_LIFETIME_EXE=shirley_lifetime$(DEBUG).x
#SHIRLEY_XLINEWIDTH_EXE=shirley_xlinewidth$(DEBUG).x
#SHIRLEY_JDOS_EXE=shirley_jdos$(DEBUG).x
#PWMAT_TRANSPOSE_EXE=pwmat_transpose$(DEBUG).x
#SHIRLEY_EPSILON_EXE=shirley_epsilon$(DEBUG).x
#SHIRLEY_EPSTCI_EXE=shirley_epstci$(DEBUG).x
#SHIRLEY_EPSILON_TRIPLES_EXE=shirley_epsilon_triples$(DEBUG).x
#SHIRLEY_EPSILON00_EXE=shirley_epsilon00$(DEBUG).x
#SHIRLEY_EPSILONIJ_EXE=shirley_epsilonij$(DEBUG).x
#SHIRLEY_EPSILONI0_EXE=shirley_epsiloni0$(DEBUG).x
#SHIRLEY_EPSILONQ_EXE=shirley_epsilonq$(DEBUG).x
#TRIPLES_TEST_EXE=triples_test$(DEBUG).x
#SHIRLEY_SELFENERGY_TRIPLES_EXE=shirley_selfenergy_triples$(DEBUG).x
#SHIRLEY_EPSILON_CONSTANT_EXE=shirley_epsilon_constant$(DEBUG).x
#DIAG_EPSILON_EXE=diag_epsilon$(DEBUG).x
#CHANGE_RHOGRID_EXE=change_rhogrid$(DEBUG).x
#SHIRLEY_STIELTJES_EXE=shirley_stieltjes$(DEBUG).x
#WRITETEST_EXE=writetest$(DEBUG).x

TLDEPS= bindir mods libs pw

all : tldeps \
$(SHIRLEY_BASIS_EXE) \
$(SHIRLEY_HAM_EXE) \
$(SHIRLEY_QDIAGP_EXE) \
$(BANDSTRUCT_EXE) \
$(EFERMI_EXE) \
$(SHIRLEY_XAS_EXE) \
$(XAS_PARA_EXE) \
$(SHIRLEY_OVERLAP_EXE) \
#$(XARS_PARA_EXE) \
#$(SHIRLEY_EIG_EXE) \
#$(SHIRLEY_PROJ_EXE) \
#$(SHIRLEY_PLOT_EXE) \
#$(PDOS_PARA_EXE) \
#$(SHIRLEY_EFERMI_EXE) \
#$(SHIRLEY_DHDK_EXE) \
#$(COND_PARA_EXE) \
#$(CHANGE_RHOGRID_EXE) \
#$(PLASMA_FREQ_EXE) \
#$(SHIRLEY_EPSILON_TRIPLES_EXE) \
#$(SHIRLEY_EPSILON00_EXE) \
#$(SHIRLEY_EPSILONIJ_EXE) \
#$(SHIRLEY_EPSILONI0_EXE) \
#$(SHIRLEY_EPSILONQ_EXE) \
#$(SHIRLEY_SELFENERGY_TRIPLES_EXE) \
#$(SHIRLEY_XASRESP_EXE)

#$(SHIRLEY_QMOMP_EXE) \
#$(SHIRLEY_EFERMI_EXE) \
#$(SHIRLEY_LIFETIME_EXE) \
#$(SHIRLEY_XLINEWIDTH_EXE) \
#$(SHIRLEY_JDOS_EXE) \
#$(PWMAT_TRANSPOSE_EXE) 
#
##$(SHIRLEY_EPSILON_EXE) \
##$(DIAG_EPSILON_EXE)


#$(SHIRLEY_HAM_EXE) : $(SHIRLEY_HAM) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(SHIRLEY_HAM) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(SHIRLEY_OVERLAP_EXE) : $(SHIRLEY_OVERLAP) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(SHIRLEY_OVERLAP) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
##$(SHIRLEY_QDIAG_EXE) : $(SHIRLEY_QDIAG) $(PW_QDIAG)
##	$(MPIF90) $(LDFLAGS) -o $@ \
##	       $(SHIRLEY_QDIAG) $(PW_QDIAG) $(LIBS)
##	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(SHIRLEY_QDIAGP_EXE) : $(SHIRLEY_QDIAGP) $(PW_QDIAGP) $(MODULES) $(LIBOBJS)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(SHIRLEY_QDIAGP) $(PW_QDIAGP) $(MODULES) $(LIBOBJS) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(SHIRLEY_PLOT_EXE) : $(SHIRLEY_PLOT) $(PW_QDIAGP) $(PPOBJS)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(SHIRLEY_PLOT) $(PW_QDIAGP) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(SHIRLEY_XAS_EXE) : $(SHIRLEY_XAS) $(PW_QDIAGP) $(MODULES) $(LIBOBJS)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(SHIRLEY_XAS) $(PW_QDIAGP) $(MODULES) $(LIBOBJS) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(SHIRLEY_XASRESP_EXE) : $(SHIRLEY_XASRESP) $(PW_QDIAGP)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(SHIRLEY_XASRESP) $(PW_QDIAGP) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(SHIRLEY_EIG_EXE) : $(SHIRLEY_EIG) $(PW_QDIAGP)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(SHIRLEY_EIG) $(PW_QDIAGP) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(SHIRLEY_PROJ_EXE) : $(SHIRLEY_PROJ) $(PW_QDIAGP)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(SHIRLEY_PROJ) $(PW_QDIAGP) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(SHIRLEY_MCSAMP_EXE) : $(SHIRLEY_MCSAMP) $(PW_QDIAGP)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(SHIRLEY_MCSAMP) $(PW_QDIAGP) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(SHIRLEY_MOM_EXE) : $(SHIRLEY_MOM) $(PW_QDIAGP)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(SHIRLEY_MOM) $(PW_QDIAGP) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(DOS_PARA_EXE) : $(DOS_PARA) $(PW_QDIAGP)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(DOS_PARA) $(PW_QDIAGP) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(PDOS_PARA_EXE) : $(PDOS_PARA) $(PW_QDIAGP)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(PDOS_PARA) $(PW_QDIAGP) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(JDOS_PARA_EXE) : $(JDOS_PARA) $(PW_QDIAGP)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(JDOS_PARA) $(PW_QDIAGP) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(XAS_PARA_EXE) : $(XAS_PARA) $(PW_QDIAGP) $(MODULES) $(LIBOBJS)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(XAS_PARA) $(PW_QDIAGP) $(MODULES) $(LIBOBJS) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(XARS_PARA_EXE) : $(XARS_PARA) $(PW_QDIAGP) $(MODULES) $(LIBOBJS)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(XARS_PARA) $(PW_QDIAGP) $(MODULES) $(LIBOBJS) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(BANDSTRUCT_EXE) : $(BANDSTRUCT) $(PW_QDIAGP) $(MODULES) $(LIBOBJS)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(BANDSTRUCT) $(PW_QDIAGP) $(MODULES) $(LIBOBJS) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(COND_PARA_EXE) : $(COND_PARA) $(PW_QDIAGP)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(COND_PARA) $(PW_QDIAGP) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(EFERMI_EXE) : $(EFERMI) $(PW_QDIAGP) $(MODULES) $(LIBOBJS)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(EFERMI) $(PW_QDIAGP) $(MODULES) $(LIBOBJS) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(PLASMA_FREQ_EXE) : $(PLASMA_FREQ) $(PW_QDIAGP)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(PLASMA_FREQ) $(PW_QDIAGP) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(GETTEST_EXE) : $(GETTEST) $(PW_QDIAGP)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(GETTEST) $(PW_QDIAGP) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(SHIRLEY_OPT_EXE) : $(SHIRLEY_OPT) $(PW_QDIAGP)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(SHIRLEY_OPT) $(PW_QDIAGP) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(SHIRLEY_QMOMP_EXE) : $(SHIRLEY_QMOMP) $(PW_QDIAGP)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(SHIRLEY_QMOMP) $(PW_QDIAGP) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(SHIRLEY_EFERMI_EXE) : $(SHIRLEY_EFERMI) $(PW_QDIAGP)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(SHIRLEY_EFERMI) $(PW_QDIAGP) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(SHIRLEY_DHDK_EXE) : $(SHIRLEY_DHDK) $(PW_QDIAGP)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(SHIRLEY_DHDK) $(PW_QDIAGP) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(SHIRLEY_LIFETIME_EXE) : $(SHIRLEY_LIFETIME) $(PWOBJS) $(MODULES) $(LIBOBJS)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(SHIRLEY_LIFETIME) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(SHIRLEY_XLINEWIDTH_EXE) : $(SHIRLEY_XLINEWIDTH) $(PWOBJS) $(MODULES) $(LIBOBJS)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(SHIRLEY_XLINEWIDTH) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(SHIRLEY_JDOS_EXE) : $(SHIRLEY_JDOS) $(PWOBJS) $(MODULES) $(LIBOBJS)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(SHIRLEY_JDOS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(PWMAT_TRANSPOSE_EXE) : $(PWMAT_TRANSPOSE) $(MODULES) $(LIBOBJS)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(PWMAT_TRANSPOSE) $(MODULES) $(LIBOBJS) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(SHIRLEY_EPSILON_EXE) : $(SHIRLEY_EPSILON) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(SHIRLEY_EPSILON) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(SHIRLEY_EPSTCI_EXE) : $(SHIRLEY_EPSTCI) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(SHIRLEY_EPSTCI) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(SHIRLEY_EPSILON_TRIPLES_EXE) : $(SHIRLEY_EPSILON_TRIPLES) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(SHIRLEY_EPSILON_TRIPLES) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(SHIRLEY_EPSILON00_EXE) : $(SHIRLEY_EPSILON00) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(SHIRLEY_EPSILON00) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(SHIRLEY_EPSILONIJ_EXE) : $(SHIRLEY_EPSILONIJ) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(SHIRLEY_EPSILONIJ) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(SHIRLEY_EPSILONI0_EXE) : $(SHIRLEY_EPSILONI0) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(SHIRLEY_EPSILONI0) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(SHIRLEY_EPSILONQ_EXE) : $(SHIRLEY_EPSILONQ) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(SHIRLEY_EPSILONQ) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(TRIPLES_TEST_EXE) : $(TRIPLES_TEST) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(TRIPLES_TEST) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(SHIRLEY_SELFENERGY_TRIPLES_EXE) : $(SHIRLEY_SELFENERGY_TRIPLES) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(SHIRLEY_SELFENERGY_TRIPLES) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(SHIRLEY_EPSILON_CONSTANT_EXE) : $(SHIRLEY_EPSILON_CONSTANT) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(SHIRLEY_EPSILON_CONSTANT) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(DIAG_EPSILON_EXE) : $(DIAG_EPSILON) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(SHIRLEY_EPSILON) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(CHANGE_RHOGRID_EXE) : $(CHANGE_RHOGRID) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(CHANGE_RHOGRID) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(SHIRLEY_STIELTJES_EXE) : $(SHIRLEY_STIELTJES) $(PW_QDIAGP)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(SHIRLEY_STIELTJES) $(PW_QDIAGP) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )
#
#$(WRITETEST_EXE) : $(WRITETEST) $(PWOBJS) $(MODULES) $(LIBOBJS)
#	$(MPIF90) $(LDFLAGS) -o $@ \
#	       $(WRITETEST) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
#	- ( cd ../bin ; $(COPY) ../SHIRLEY/$@ . )

$(SHIRLEY_BASIS_EXE) : $(SHIRLEY_BASIS) $(QEMODS) $(LIBOBJS) 
	$(LD) $(LDFLAGS) -o $@ \
	       $(SHIRLEY_BASIS) $(QEMODS) $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; $(COPY) ../SHIRLEY/src/$@ . )

$(SHIRLEY_HAM_EXE) : $(SHIRLEY_HAM) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
	       $(SHIRLEY_HAM) $(QEMODS) $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; $(COPY) ../SHIRLEY/src/$@ . )

$(SHIRLEY_QDIAGP_EXE) : $(SHIRLEY_QDIAGP) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
	       $(SHIRLEY_QDIAGP) $(QEMODS) $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; $(COPY) ../SHIRLEY/src/$@ . )

$(BANDSTRUCT_EXE) : $(BANDSTRUCT) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
	       $(BANDSTRUCT) $(QEMODS) $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; $(COPY) ../SHIRLEY/src/$@ . )

$(SHIRLEY_XAS_EXE) : $(SHIRLEY_XAS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
	       $(SHIRLEY_XAS) $(QEMODS) $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; $(COPY) ../SHIRLEY/src/$@ . )

$(XAS_PARA_EXE) : $(XAS_PARA) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
	       $(XAS_PARA) $(QEMODS) $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; $(COPY) ../SHIRLEY/src/$@ . )

$(EFERMI_EXE) : $(EFERMI) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
	       $(EFERMI) $(QEMODS) $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; $(COPY) ../SHIRLEY/src/$@ . )

$(SHIRLEY_OVERLAP_EXE) : $(SHIRLEY_OVERLAP) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
	       $(SHIRLEY_OVERLAP) $(QEMODS) $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; $(COPY) ../SHIRLEY/src/$@ . )

tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i work.pc

include make.depend
# DO NOT DELETE
