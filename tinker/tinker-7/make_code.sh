#! /usr/bin/env bash

check() {
if [ ! -s $1 ] ; then echo ; echo  "aborting, $1 does not exist" ; exit ; fi
}
rm -f $1.{o,x}
rm ../bin/$1
check $1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp $1.f
check $1.o
echo "Remaking library"
../linux/intel/library.make
echo "Linking"
../linux/intel/link2.make
echo "Copying to ../bin"
../linux/intel/rename.make
if [ -s $1.x ] ; then cp $1.x ../bin/$1 ; fi
