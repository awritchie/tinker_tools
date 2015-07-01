#! /bin/bash

clear
# compile
make clean
tar cfvz tinker2field.tar *cpp *hpp Makefile
rm -f tinker2field
make
if [ ! -s tinker2field ] ; then exit ; fi
#./tinker2field -h
a1=6 #6268
a2=7
prm=/Users/ritchie/software/tinker/params/webb_amoebapro13.prm



# test xyz + u
#time ./tinker2field -a1 $a1 -a2 $a2 -p $prm -x examples/em_mescn.7000 -u examples/em_mescn.7000u_2 -nt 1

# test xyz + f

# test arc + uind
#time ./tinker2field -a1 $a1 -a2 $a2 -p $prm -x examples/npt.arc -u examples/npt.uind -nt 2

# test arc + frc
#time ./tinker2field -a1 $a1 -a2 $a2 -p $prm -x examples/npt.arc -f examples/npt.frc -nt 3

# test arc + uind + frc
#time ./tinker2field -a1 $a1 -a2 $a2 -p $prm -x examples/npt.arc_2 -u examples/npt.uind -f examples/npt.frc -nt 8

# test arc + uding + frc on small file
#time ./tinker2field -a1 $a1 -a2 $a2 -p $prm -x examples/tmp.arc -u examples/tmp.uind -f examples/tmp.frc -nt 1 #8

# Big test
time ./tinker2field -p ~ritchie/software/tinker/params/webb_amoebapro13.prm  -x examples/em_mescn.arc -u examples/em_mescn.uind -a1 6 -a2 7 -o examples/em_mescn.field

# Big test, same data as ../tinker2h5/*.h5
time ./tinker2field -p ~ritchie/software/tinker/params/webb_amoebapro13.prm  -x examples/mescn.arc -u examples/mescn.uind  -a1 6 -a2 7 -o examples/mescn.field

# test hdf5
time ./tinker2field -a1 $a1 -a2 $a2 -p $prm -x ../tinker2h5/none.h5 -u ../tinker2h5/none.h5  -nt 1 -o ../tinker2h5/none.field
# test hdf5 with gzip
time ./tinker2field -a1 $a1 -a2 $a2 -p $prm -x ../tinker2h5/gzip.h5 -u ../tinker2h5/gzip.h5  -nt 1 -o ../tinker2h5/gzip.field
# test hdf5 with lzf
#time ./tinker2field -a1 $a1 -a2 $a2 -p $prm -x ../tinker2h5/lzf.h5 -u ../tinker2h5/lzf.h5 -f ../tinker2h5/lzf.h5 -nt 1 -o ../tinker2h5/lzf.field