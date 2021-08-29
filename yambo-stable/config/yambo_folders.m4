#
chmod u+x sbin/*
cp config/Makefile .
#
if ! test -d bin      ; then mkdir bin      ; fi
if ! test -d include  ; then mkdir include  ; fi
if ! test -d lib      ; then mkdir lib      ; fi
#
if   test -d include/system ; then rm -r include/system ; fi
if ! test -d include/system ; then mkdir include/system ; fi
#
# Copy system headers
#
IOTK_INCDIRS=`echo $IOTK_INCS       | sed "s/$IFLAG/ /g"` ; 
for includedir in $IOTK_INCDIRS;   do
  count=`ls -1 $includedir/*iotk*.h   2>/dev/null | wc -l`
  if test $count != 0; then ln -fs $includedir/*iotk*.h   include/system/  ; fi
  count=`ls -1 $includedir/*iotk*.mod 2>/dev/null | wc -l`
  if test $count != 0; then ln -fs $includedir/*iotk*.mod include/system/  ; fi
done
#
ETSF_INCDIRS=`echo "$ETSF_INCS"     | sed "s/$IFLAG/ /g"`
for includedir in $ETSF_INCDIRS;   do
  count=`ls -1 $includedir/*etsf*.h   2>/dev/null | wc -l`
  if test $count != 0; then ln -fs $includedir/*etsf*.h   include/system/  ; fi
  count=`ls -1 $includedir/*etsf*.mod 2>/dev/null | wc -l`
  if test $count != 0; then ln -fs $includedir/*etsf*.mod include/system/  ; fi
done
#
NETCDF_INCDIRS=`echo "$NETCDF_INCS" | sed "s/$IFLAG/ /g"`
NETCDFF_INCDIRS=`echo "$NETCDFF_INCS" | sed "s/$IFLAG/ /g"`
for includedir in $NETCDF_INCDIRS $NETCDFF_INCDIRS; do
  count=`ls -1 $includedir/*netcdf*.h   2>/dev/null | wc -l`
  if test $count != 0; then ln -fs $includedir/*netcdf*.h   include/system/  ; fi
  count=`ls -1 $includedir/*netcdf*.mod   2>/dev/null | wc -l`
  if test $count != 0; then ln -fs $includedir/*netcdf*.mod include/system/  ; fi
  if test -e $includedir/typesizes.h; then ln -fs $includedir/typesizes.h  include/system/  ; fi
done
#
HDF5_INCDIRS=`echo "$HDF5_INCS"     | sed "s/$IFLAG/ /g"`
for includedir in $HDF5_INCDIRS;   do
  count=`ls -1 $includedir/*H5*.h   2>/dev/null | wc -l`
  if test $count != 0; then ln -fs $includedir/*H5*.h       include/system/  ; fi
  count=`ls -1 $includedir/*H5*.mod 2>/dev/null | wc -l`
  if test $count != 0; then ln -fs $includedir/*H5*.mod     include/system/  ; fi
done
#
FFT_INCDIRS=`echo "$FFT_INCS"       | sed "s/$IFLAG/ /g"`
for includedir in $FFT_INCDIRS;    do
  count=`ls -1 $includedir/*fft*.h   2>/dev/null | wc -l`
  if test $count != 0; then ln -fs $includedir/*fft*.h      include/system/  ; fi
  count=`ls -1 $includedir/*fft*.mod 2>/dev/null | wc -l`
  if test $count != 0; then ln -fs $includedir/*fft*.mod    include/system/  ; fi
done
#
PETSC_INCDIRS=`echo "$PETSC_INCS"   | sed "s/$IFLAG/ /g"`
for includedir in $PETSC_INCDIRS;  do
  count=`ls -1 $includedir/*petsc*.h   2>/dev/null | wc -l`
  if test $count != 0; then ln -fs $includedir/*petsc*.h    include/system/  ; fi
  count=`ls -1 $includedir/*petsc*.mod 2>/dev/null | wc -l`
  if test $count != 0; then ln -fs $includedir/*petsc*.mod  include/system/  ; fi
done
#
SLEPC_INCDIRS=`echo "$SLEPC_INCS"   | sed "s/$IFLAG/ /g"`
for includedir in $SLEPC_INCDIRS;  do
  count=`ls -1 $includedir/*slepc*.h   2>/dev/null | wc -l`
  if test $count != 0; then ln -fs $includedir/*slepc*.h    include/system/  ; fi
  count=`ls -1 $includedir/*slepc*.mod 2>/dev/null | wc -l`
  if test $count != 0; then ln -fs $includedir/*slepc*.mod  include/system/  ; fi
done
#
LIBXC_INCDIRS=`echo "$LIBXC_INCS"   | sed "s/$IFLAG/ /g"`
for includedir in $LIBXC_INCDIRS;  do
  count=`ls -1 $includedir/*xc*.h   2>/dev/null | wc -l`
  if test $count != 0; then ln -fs $includedir/*xc*.h       include/system/  ; fi
  count=`ls -1 $includedir/*xc*.mod 2>/dev/null | wc -l`
  if test $count != 0; then ln -fs $includedir/*xc*.mod     include/system/  ; fi
done
#
MPI_INCDIRS=`echo "$MPI_INCS"       | sed "s/$IFLAG/ /g"`
for includedir in $MPI_INCDIRS;    do
  count=`ls -1 $includedir/*mpi*.h   2>/dev/null | wc -l`
  if test $count != 0; then ln -fs $includedir/*mpi*.h     include/system/  ; fi
  count=`ls -1 $includedir/*mpi*.mod 2>/dev/null | wc -l`
  if test $count != 0; then ln -fs $includedir/*mpi*.mod   include/system/  ; fi
done
#
if ! test -d "$extlibs_path/${FCKIND}/${FC}";         then mkdir -p "$extlibs_path/${FCKIND}/${FC}";         fi
if ! test -d "$extlibs_path/${FCKIND}/${FC}/lib";     then mkdir    "$extlibs_path/${FCKIND}/${FC}/lib";     fi
if ! test -d "$extlibs_path/${FCKIND}/${FC}/include"; then mkdir    "$extlibs_path/${FCKIND}/${FC}/include"; fi
if ! test -d "$extlibs_path/${FCKIND}/${FC}/bin";     then mkdir    "$extlibs_path/${FCKIND}/${FC}/bin";     fi
#
cp $extlibs_path/${FCKIND}/${FC}/bin/* $exec_prefix/bin/
cp $extlibs_path/${FCKIND}/${FC}/v3/bin/* $exec_prefix/bin/
#
if [[ "$prefix" != "$srcdir" ]] && [[ "$srcdir" != "." ]] ; then
 if test ! -d "$prefix/driver"     ; then mkdir "$prefix/driver"     ; fi
 if test ! -d "$prefix/ypp"        ; then mkdir "$prefix/ypp"        ; fi
 if test ! -d "$prefix/interfaces" ; then mkdir "$prefix/interfaces" ; fi
 if test ! -d "$prefix/lib"        ; then mkdir "$prefix/lib"        ; fi
 if test ! -d "$prefix/lib/local"  ; then mkdir "$prefix/lib/local"  ; fi
 if test ! -d "$prefix/lib/slatec" ; then mkdir "$prefix/lib/slatec" ; fi
 cd "$srcdir/driver/" ;
 for file in `ls *.h` ; do
  cp "$file" "$prefix/driver" ;
 done ;
 for file in `ls *.c` ; do
  cp "$file" "$prefix/driver" ;
 done ;
 cd "../lib/local" ;
 for file in `ls *.f` ; do
  cp "$file" "$prefix/lib/local" ;
 done ;
 cd "../slatec" ;
 for file in `ls *.f` ; do
  cp "$file" "$prefix/lib/slatec" ;
 done ;
 cd "$prefix" ;
fi
