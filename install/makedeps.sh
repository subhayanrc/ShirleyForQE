#!/bin/sh
# compute dependencies for the PWscf directory tree

# make sure there is no locale setting creating unneeded differences.
LC_ALL=C
export LC_ALL

# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
TOPDIR=`pwd`

if test $# = 0
then
    dirs=" LAXlib FFTXlib UtilXlib Modules clib LR_Modules upftools \
           KS_Solvers/Davidson KS_Solvers/Davidson_RCI KS_Solvers/CG KS_Solvers/PPCG \
           PW/src CPV/src PW/tools upftools PP/src PWCOND/src \
           PHonon/Gamma PHonon/PH PHonon/FD HP/src atomic/src \
           EPW/src XSpectra/src ACFDT/src NEB/src TDDFPT/src \
           GWW/pw4gww GWW/gww GWW/head GWW/bse GWW/simple \
	   GWW/simple_bse GWW/simple_ip \
           SHIRLEY/src corerepair" 
          
elif
    test $1 = "-addson" 
then
    echo "The script for adding new dependencies is running"
    echo "Usage: $0 -addson DIR DEPENDENCY_DIRS"
    echo "$0 assumes that the new dependencies are in $TOPDIR/../"
#    ninput=$#
#    echo "number of input arguments: $ninput"
    dirs=$2
    shift
    shift
    add_deps=$*
    echo "dependencies in $add_deps will be searched for $dirs"
else
    dirs=$*
fi


for dir in $dirs; do

    # the following command removes a trailing slash
    DIR=`echo ${dir%/}`

    # the following would also work
    #DIR=`echo $dir | sed "s,/$,,"`

    # set inter-directory dependencies - only directories containing
    # modules that are used, or files that are included, by routines
    # in directory DIR should be listed in DEPENDS
    LEVEL1=..
    LEVEL2=../..
    # default
    DEPENDS="$LEVEL1/include" 
    # for convenience, used later
    DEPEND1="$LEVEL1/include $LEVEL1/iotk/src $LEVEL1/FFTXlib $LEVEL1/LAXlib $LEVEL1/UtilXlib"
    DEPEND2="$LEVEL2/include $LEVEL2/iotk/src $LEVEL2/FFTXlib $LEVEL2/LAXlib $LEVEL2/UtilXlib \
             $LEVEL2/Modules"
    DEPEND3="$LEVEL2/include $LEVEL2/FFTXlib $LEVEL2/LAXlib $LEVEL2/UtilXlib"
    case $DIR in 
        Modules )
             DEPENDS="$DEPEND1 $LEVEL1/UtilXlib" ;;
        LAXlib )
             DEPENDS="$LEVEL1/UtilXlib " ;;
        upftools )
             DEPENDS="$DEPEND1 $LEVEL1/Modules" ;;
        LR_Modules )
             DEPENDS="$DEPEND1 $LEVEL1/Modules $LEVEL1/PW/src" ;;
	ACFDT/src ) 
             DEPENDS="$DEPEND2 $LEVEL2/PW/src $LEVEL2/PHonon/PH $LEVEL2/LR_Modules" ;;
	atomic/src | GWW/gww )
	     DEPENDS="$DEPEND2" ;;
	PW/src | CPV/src )
	     DEPENDS="$DEPEND2 ../../KS_Solvers/Davidson ../../KS_Solvers/CG ../../KS_Solvers/PPCG ../../dft-d3" ;;
	KS_Solvers/Davidson | KS_Solvers/Davidson_RCI | KS_Solvers/CG | KS_Solvers/PPCG )
	     DEPENDS="$DEPEND3" ;;
	PW/tools | PP/src | PWCOND/src | GWW/pw4gww | NEB/src )
	     DEPENDS="$DEPEND2 $LEVEL2/PW/src" ;;
	PHonon/FD | PHonon/PH | PHonon/Gamma | HP/src | TDDFPT/src | XSpectra/src  | GIPAW/src )
	     DEPENDS="$DEPEND2 $LEVEL2/PW/src $LEVEL2/LR_Modules" ;;
        EPW/src )
             DEPENDS="$DEPEND2 $LEVEL2/PW/src $LEVEL2/LR_Modules $LEVEL2/PHonon/PH $LEVEL2/Modules" ;; 
	GWW/head )
	     DEPENDS="$DEPEND2 $LEVEL2/PW/src $LEVEL2/PHonon/PH $LEVEL2/LR_Modules" ;;	
	GWW/bse )
	 DEPENDS="$DEPEND2 $LEVEL2/PW/src $LEVEL2/PHonon/PH $LEVEL2/LR_Modules $LEVEL2/GWW/pw4gww $LEVEL2/GWW/gww" ;;	
	GWW/simple )
	 DEPENDS="$DEPEND2 $LEVEL2/PW/src $LEVEL2/GWW/pw4gww $LEVEL2/GWW/gww" ;;
	GWW/simple_bse )
	 DEPENDS="$DEPEND2 $LEVEL2/GWW/gww" ;;
	GWW/simple_ip)
	DEPENDS="$DEPEND2" ;;
	SHIRLEY/src)
	DEPENDS="$LEVEL2/Modules $LEVEL2/UtilXlib $LEVEL2/FFTXlib $LEVEL2/PW/src" ;;
	corerepair)
	DEPENDS="$LEVEL1/Modules $LEVEL1/UtilXlib $LEVEL1/FFTXlib $LEVEL1/upftools" ;;
    *)
# if addson needs a make.depend file
	DEPENDS="$DEPENDS $add_deps"

    esac

    # generate dependencies file (only for directories that are present)
    if test -d $TOPDIR/../$DIR
    then
	cd $TOPDIR/../$DIR
       
	$TOPDIR/moduledep.sh $DEPENDS > make.depend
	$TOPDIR/includedep.sh $DEPENDS >> make.depend

        # handle special cases: modules for C-fortran binding, hdf5, MPI
        sed '/@iso_c_binding@/d' make.depend > make.depend.tmp
        sed '/@hdf5@/d;/@mpi@/d' make.depend.tmp > make.depend
        sed '/@fox_dom@/d;/@fox_wxml@/d'  make.depend > make.depend.tmp 
        sed '/@m_common_io@/d'   make.depend.tmp > make.depend

        if test "$DIR" = "FFTXlib"
        then
            sed '/@mkl_dfti/d' make.depend > make.depend.tmp
            sed '/@fftw3.f/d;s/@fftw.c@/fftw.c/' make.depend.tmp > make.depend
        fi

        if test "$DIR" = "LAXlib"
        then
            sed '/@elpa1@/d' make.depend > make.depend.tmp
            cp make.depend.tmp make.depend
        fi

        if test "$DIR" = "UtilXlib"
        then
            sed '/@elpa1@/d' make.depend > make.depend.tmp
            sed '/@ifcore@/d' make.depend.tmp > make.depend
        fi

        if test "$DIR" = "KS_Solvers/Davidson"
        then

            sed '/@elpa1@/d' make.depend > make.depend.tmp
            sed '/@ifcore@/d' make.depend.tmp > make.depend
        fi

        if test "$DIR" = "KS_Solvers/Davidson_RCI"
        then
            sed '/@elpa1@/d' make.depend > make.depend.tmp
            sed '/@ifcore@/d' make.depend.tmp > make.depend
        fi

        if test "$DIR" = "KS_Solvers/CG"
        then

            sed '/@elpa1@/d' make.depend > make.depend.tmp
            sed '/@ifcore@/d' make.depend.tmp > make.depend
        fi

        if test "$DIR" = "KS_Solvers/PPCG"
        then

            sed '/@elpa1@/d' make.depend > make.depend.tmp
            sed '/@ifcore@/d' make.depend.tmp > make.depend
        fi

        if test "$DIR" = "Modules"
        then

            sed '/@elpa1@/d' make.depend > make.depend.tmp
            sed '/@ifcore@/d' make.depend.tmp > make.depend
        fi

        if test "$DIR" = "PW/src" || test "$DIR" = "TDDFPT/src"
        then
            sed '/@environ_/d;/@solvent_tddfpt@/d' make.depend > make.depend.tmp
            cp make.depend.tmp make.depend
        fi

        if test "$DIR" = "CPV/src"
        then
            sed '/@f90_unix_proc@/d' make.depend > make.depend.tmp
            cp make.depend.tmp make.depend
        fi

        if test "$DIR" = "EPW/src"
        then
            sed '/@f90_unix_io@/d' make.depend > make.depend.tmp
            cp make.depend.tmp make.depend
            sed '/@f90_unix_env@/d' make.depend > make.depend.tmp
            cp make.depend.tmp make.depend
            sed '/@w90_io@/d' make.depend > make.depend.tmp
            cp make.depend.tmp make.depend
            sed '/@ifport@/d' make.depend > make.depend.tmp
            cp make.depend.tmp make.depend
        fi

        rm -f make.depend.tmp

        # check for missing dependencies 
        if grep @ make.depend
        then
	   notfound=1
	   echo WARNING: dependencies not found in directory $DIR
       else
           echo directory $DIR : ok
       fi
    else
       echo directory $DIR : not present in $TOPDIR 
    fi
done
if test "$notfound" = ""
then
    echo all dependencies updated successfully
fi
