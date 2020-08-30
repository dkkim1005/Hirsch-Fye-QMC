#!/bin/bash
TARGET=qmc
LIBFFTW3='-lfftw3'
LIBTRNG4='-ltrng4'
# fill blas lapack link options in 'LIBBLASLAPACK'
LIBBLASLAPACK=DEFAULT
if [ $LIBBLASLAPACK = DEFAULT ]; then
  LIBBLASLAPACK='-llapack -lblas'
fi

#####
SRCS=(hirschfye_qmc green_tools etc main)
OBJ=""
for SRC in ${SRCS[@]}; do
  if [ -e ${SRC}.o ]; then
    continue
  fi
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
