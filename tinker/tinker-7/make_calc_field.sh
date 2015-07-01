#! /usr/bin/env bash

check() {
if [ ! -s $1 ] ; then echo  "aborting, $1 does not exist" ; exit ; fi
}
clear
rm -f calc_field.{o,x}
rm -f ../bin/calc_field
check calc_field.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp calc_field.f
check calc_field.o
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o calc_field.x calc_field.o libtinker.a ../fftw/lib/libfftw3_threads.a ../fftw/lib/libfftw3.a ; strip calc_field.x
check calc_field.x
cp -v calc_field.x   ../bin/calc_field
