#! /usr/bin/env bash

check() {
if [ ! -s $1 ] ; then echo  "aborting, $1 does not exist" ; exit ; fi
}
rm -f analyze-dip.{o,x}
rm ../bin/analyze-dip
check analyze-dip.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp analyze-dip.f
check analyze-dip.o
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o analyze-dip.x analyze-dip.o libtinker.a ../fftw/lib/libfftw3_threads.a ../fftw/lib/libfftw3.a ; strip analyze-dip.x
check analyze-dip.x
mv analyze-dip.x   ../bin/analyze-dip
