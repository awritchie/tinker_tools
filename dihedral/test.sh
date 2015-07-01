#! /bin/bash

arc=/Volumes/Alfheim/ritchie/MD/AMOEBA/hxscn/D2O/4/hxscn.arc
ndx=ndx.txt
out=test.dih

rm dihedral
clear
make
if [ ! -s dihedral ] ; then exit ; fi
./dihedral -x $arc -n $ndx -o $out -e 10
tail $out
echo 
