#!/bin/bash
TARGET=qmc
LIBFFTW3='-lfftw3'
LIBTRNG4='-ltrng4'
LIBBLASLAPACK='-lopenblas -llapack -L/usr/local/opt/lapack/lib -L/usr/local/opt/openblas/lib'

#####
SRCS=(hirschfye_qmc green_tools main)
OBJ=""
for SRC in ${SRCS[@]}; do
  if [ -e ${SRC}.cpp ]; then
    g++ -c -o ${SRC}.o ${SRC}.cpp -std=c++11
    OBJ="$SRC.o "$OBJ
  else
    echo " plz check: $SRC.cpp"
    exit 1
  fi
done

g++ -o $TARGET $OBJ -std=c++11 $LIBFFTW3 $LIBTRNG4 $LIBBLASLAPACK
rm *.o
