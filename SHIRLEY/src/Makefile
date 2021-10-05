# Makefile for postprocessing (PP)

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS = $(BASEMOD_FLAGS) \
           $(MOD_FLAG)../../PW/src \
           $(MOD_FLAG)../../FFTXlib

QEMODS=../../Modules/libqemod.a \
       ../../PW/src/libpw.a \
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

#bz.o \
check_energies.o \
read_file_shirley.o \

SHIRLEY_HAM = \
shirley_ham.o \
shirley_ham_input.o \
start_shirley.o \
stop_shirley.o \
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

#B0.o \
#elph_shirley.o \
#get_dnl.o \
#get_proj_ham.o \
#gk_l2gmap.o \
#interaction.o \
#interaction_test.o \
#kpt.o \
#qinteraction.o \
#plot.o \
#plot_basis.o \
#pole_strength.o \
#proj_shirley.o \
#pwmat.o \
#read_file_shirley.o \
#rgrid.o \
#triples.o \
#trnvect_wrapper.o \
#
#SHIRLEY_OVERLAP = \
#append_missing_slash.o \
#print_clock_shirley.o \
#shirley_overlap.o \
#shirley_overlap_input.o \
#start_shirley.o \
#stop_shirley.o 
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
#SHIRLEY_XAS = shirley_xas.o $(NEWDIAGOBJ)
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
#XAS_PARA = xas_para.o start_shirley.o freeunit.o mpio.o
#XARS_PARA = xars_para.o start_shirley.o freeunit.o mpio.o
#COND_PARA = cond_para.o fermi.o start_shirley.o freeunit.o mpio.o \
#            ../PW/efermit.o ../PW/sumkt.o ../flib/dost.o ../flib/sort.o
#EFERMI = efermi.o fermi.o \
#            ../PW/efermit.o ../PW/sumkt.o ../flib/dost.o ../flib/sort.o \
#            $(NEWDIAGOBJ)
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
#
#MODULES = \
#../Modules/atom.o \
#../Modules/basic_algebra_routines.o \
#../Modules/cell_base.o \
#../Modules/check_stop.o \
#../Modules/clocks.o \
#../Modules/constants.o \
#../Modules/constraints_module.o \
#../Modules/control_flags.o \
#../Modules/coulomb_vcut.o \
#../Modules/descriptors.o \
#../Modules/dspev_drv.o \
#../Modules/electrons_base.o \
#../Modules/environment.o \
#../Modules/error_handler.o \
#../Modules/fft_base.o \
#../Modules/fft_interfaces.o \
#../Modules/fft_parallel.o \
#../Modules/fft_scalar.o \
#../Modules/fft_types.o \
#../Modules/funct.o \
#../Modules/griddim.o \
#../Modules/image_io_routines.o \
#../Modules/input_parameters.o \
#../Modules/io_files.o \
#../Modules/io_global.o \
#../Modules/ions_base.o \
#../Modules/kind.o \
#../Modules/mm_dispersion.o \
#../Modules/mp.o \
#../Modules/mp_base.o \
#../Modules/mp_global.o \
#../Modules/mp_wave.o \
#../Modules/parallel_include.o \
#../Modules/parameters.o \
#../Modules/parser.o \
#../Modules/paw_variables.o \
#../Modules/printout_base.o \
#../Modules/pseudo_types.o \
#../Modules/ptoolkit.o \
#../Modules/radial_grids.o \
#../Modules/random_numbers.o \
#../Modules/read_ncpp.o \
#../Modules/read_upf_v1.o \
#../Modules/read_upf_v2.o \
#../Modules/read_uspp.o \
#../Modules/recvec.o \
#../Modules/recvec_subs.o \
#../Modules/splinelib.o \
#../Modules/stick_base.o \
#../Modules/stick_set.o \
#../Modules/timestep.o \
#../Modules/uspp.o \
#../Modules/upf.o \
#../Modules/upf_to_internal.o \
#../Modules/version.o \
#../Modules/wannier.o \
#../Modules/wannier_new.o \
#../Modules/wavefunctions.o \
#../Modules/wrappers.o \
#../Modules/write_upf_v2.o \
#../Modules/ws_base.o \
#../Modules/xml_io_base.o \
#../Modules/zhpev_drv.o 
#
#PWOBJS = \
#../PW/a2fmod.o \
#../PW/add_bfield.o \
#../PW/add_efield.o \
#../PW/add_paw_to_deeq.o \
#../PW/add_vuspsi.o \
#../PW/addusdens.o \
#../PW/allocate_bp_efield.o \
#../PW/allocate_fft.o \
#../PW/allocate_locpot.o \
#../PW/allocate_nlpot.o \
#../PW/allocate_wfc.o \
#../PW/atomic_rho.o \
#../PW/atomic_wfc.o \
#../PW/average_pp.o \
#../PW/becmod.o \
#../PW/bp_c_phase.o \
#../PW/bp_calc_btq.o \
#../PW/bp_qvan3.o \
#../PW/bp_strings.o \
#../PW/buffers.o \
#../PW/ccgdiagg.o \
#../PW/cdiagh.o \
#../PW/cdiaghg.o \
#../PW/cegterg.o \
#../PW/clean_pw.o \
#../PW/close_files.o \
#../PW/compute_dip.o \
#../PW/compute_rho.o \
#../PW/compute_ux.o \
#../PW/coset.o \
#../PW/d_matrix.o \
#../PW/data_structure.o \
#../PW/deriv_drhoc.o \
#../PW/divide.o \
#../PW/divide_class.o \
#../PW/divide_class_so.o \
#../PW/divide_et_impera.o \
#../PW/dqvan2.o \
#../PW/drhoc.o \
#../PW/dvloc_of_g.o \
#../PW/dynamics_module.o \
#../PW/efermig.o \
#../PW/efermit.o \
#../PW/eqvect.o \
#../PW/ewald.o \
#../PW/ewald_dipole.o \
#../PW/exx.o \
#../PW/find_group.o \
#../PW/g_psi.o \
#../PW/g_psi_mod.o \
#../PW/g2_kin.o \
#../PW/gen_at_dj.o \
#../PW/gen_at_dy.o \
#../PW/gen_us_dj.o \
#../PW/gen_us_dy.o \
#../PW/get_locals.o \
#../PW/gk_sort.o \
#../PW/gradcorr.o \
#../PW/gweights.o \
#../PW/h_1psi.o \
#../PW/h_epsi_her_apply.o \
#../PW/h_epsi_her_set.o \
#../PW/h_psi.o \
#../PW/h_psi_meta.o \
#../PW/hinit0.o \
#../PW/hinit1.o \
#../PW/init_at_1.o \
#../PW/init_ns.o \
#../PW/init_run.o \
#../PW/init_us_1.o \
#../PW/init_us_2.o \
#../PW/init_vloc.o \
#../PW/interpolate.o \
#../PW/io_rho_xml.o \
#../PW/irrek.o \
#../PW/iweights.o \
#../PW/kernel_table.o \
#../PW/kpoint_grid.o \
#../PW/lchk_tauxk.o \
#../PW/make_pointlists.o \
#../PW/martyna_tuckerman.o \
#../PW/memory_report.o \
#../PW/multable.o \
#../PW/n_plane_waves.o \
#../PW/new_ns.o \
#../PW/new_occ.o \
#../PW/newd.o \
#../PW/noncol.o \
#../PW/ns_adj.o \
#../PW/offset_atom_wfc.o \
#../PW/openfil.o \
#../PW/ortho_wfc.o \
#../PW/orthoatwfc.o \
#../PW/output_tau.o \
#../PW/para.o \
#../PW/paw_init.o \
#../PW/paw_onecenter.o \
#../PW/paw_symmetry.o \
#../PW/potinit.o \
#../PW/print_clock_pw.o \
#../PW/punch.o \
#../PW/pwcom.o \
#../PW/pw_restart.o \
#../PW/qvan2.o \
#../PW/rcgdiagg.o \
#../PW/rdiaghg.o \
#../PW/read_pseudo.o \
#../PW/realus.o \
#../PW/report_mag.o \
#../PW/rho2zeta.o \
#../PW/rotate_wfc.o \
#../PW/rotate_wfc_gamma.o \
#../PW/rotate_wfc_k.o \
#../PW/ruotaijk.o \
#../PW/s_1psi.o \
#../PW/s_psi.o \
#../PW/save_in_cbands.o \
#../PW/save_in_electrons.o \
#../PW/save_in_ions.o \
#../PW/scale_h.o \
#../PW/scf_mod.o \
#../PW/seqopn.o \
#../PW/set_hubbard_l.o \
#../PW/set_kplusq.o \
#../PW/set_kup_and_kdw.o \
#../PW/set_rhoc.o \
#../PW/set_vrs.o \
#../PW/setlocal.o \
#../PW/setqf.o \
#../PW/setup.o \
#../PW/spinor.o \
#../PW/sph_ind.o \
#../PW/struct_fact.o \
#../PW/sum_band.o \
#../PW/sumkg.o \
#../PW/sumkt.o \
#../PW/summary.o \
#../PW/symm_base.o \
#../PW/symme.o \
#../PW/tabd.o \
#../PW/transform_becsum_nc.o \
#../PW/transform_becsum_so.o \
#../PW/trnvecc.o \
#../PW/tweights.o \
#../PW/update_pot.o \
#../PW/usnldiag.o \
#../PW/v_of_rho.o \
#../PW/vhpsi.o \
#../PW/vloc_of_g.o \
#../PW/vloc_psi.o \
#../PW/wannier_check.o \
#../PW/wannier_clean.o \
#../PW/wannier_init.o \
#../PW/weights.o \
#../PW/wfcinit.o \
#../PW/write_ns.o \
#../PW/wsweight.o \
#../PW/xc_vdW_DF.o
#
#PW_QDIAG = \
#../flib/date_and_tim.o \
#../flib/flush_unit.o \
#../PW/cdiagh.o 
#
#PPOBJS = ../PP/xsf.o ../PP/cube.o
#
##../PW/cdiagh.o \
#
#PW_QDIAGP = \
#../flib/date_and_tim.o \
#../flib/flush_unit.o \
#../flib/inpfile.o 

SHIRLEY_BASIS_EXE=shirley_basis$(DEBUG).x
SHIRLEY_HAM_EXE=shirley_ham$(DEBUG).x
SHIRLEY_QDIAGP_EXE=shirley_qdiagp$(DEBUG).x
BANDSTRUCT_EXE=bandstruct$(DEBUG).x
#SHIRLEY_OVERLAP_EXE=shirley_overlap$(DEBUG).x
##SHIRLEY_QDIAG_EXE=shirley_qdiag$(DEBUG).x
#SHIRLEY_PLOT_EXE=shirley_plot$(DEBUG).x
#SHIRLEY_XAS_EXE=shirley_xas$(DEBUG).x
#SHIRLEY_XASRESP_EXE=shirley_xasresp$(DEBUG).x
#SHIRLEY_EIG_EXE=shirley_eig$(DEBUG).x
#SHIRLEY_PROJ_EXE=shirley_proj$(DEBUG).x
#SHIRLEY_DHDK_EXE=shirley_dhdk$(DEBUG).x
#SHIRLEY_MCSAMP_EXE=shirley_mcsamp$(DEBUG).x
#SHIRLEY_MOM_EXE=shirley_mom$(DEBUG).x
#DOS_PARA_EXE=dos_para$(DEBUG).x
#PDOS_PARA_EXE=pdos_para$(DEBUG).x
#JDOS_PARA_EXE=jdos_para$(DEBUG).x
#XAS_PARA_EXE=xas_para$(DEBUG).x
#XARS_PARA_EXE=xars_para$(DEBUG).x
#COND_PARA_EXE=cond_para$(DEBUG).x
#EFERMI_EXE=efermi$(DEBUG).x
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
#$(SHIRLEY_OVERLAP_EXE) \
#$(EFERMI_EXE) \
#$(SHIRLEY_XAS_EXE) \
#$(XAS_PARA_EXE) \
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

$(SHIRLEY_BASIS_EXE) : $(SHIRLEY_BASIS) $(LIBOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	       $(SHIRLEY_BASIS) $(LIBOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../../bin ; $(COPY) ../SHIRLEY/src/$@ . )

$(SHIRLEY_HAM_EXE) : $(SHIRLEY_HAM) $(LIBOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	       $(SHIRLEY_HAM) $(LIBOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../../bin ; $(COPY) ../SHIRLEY/src/$@ . )

$(SHIRLEY_QDIAGP_EXE) : $(SHIRLEY_QDIAGP) $(LIBOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	       $(SHIRLEY_QDIAGP) $(LIBOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../../bin ; $(COPY) ../SHIRLEY/src/$@ . )

$(BANDSTRUCT_EXE) : $(BANDSTRUCT) $(LIBOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	       $(BANDSTRUCT) $(LIBOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../../bin ; $(COPY) ../SHIRLEY/src/$@ . )

tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i work.pc

include make.depend
# DO NOT DELETE
