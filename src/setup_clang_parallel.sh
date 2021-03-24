#!/bin/sh

# This file contains compiler settings that enable the code to run in parallel.
#
# WARNING: Although I don't observe any strange behavior now (2021-3-24),
#          experience has taught me that parallel versions of my code are
#          more likely to have bugs.  So if you compile the program this
#          way, please pay attention to the program's output.

export ANSI_C="clang"
export ANSI_CPP="clang++"
export L_COMP="ar rs"

export LFLAGS="-fopenmp"  #-static

export MY_FLAGS="-std=c++11 -O3 -DNDEBUG -ffast-math"
export CFLAGS="-c $MY_FLAGS -fopenmp"
export CPP_PRELINKER_COMMAND="echo"
export COMPILER_TEMP_FILES=""
export LINKER_TEMP_FILES=""
