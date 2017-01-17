#!/bin/bash


# In order to compile successfully, you need to specify the location
# of the header files and libraries using the relevant compiler flags.
# In this particular case, all header files apart from the ones for
# the HYPRE-package are located in ${DEVINCL} and the respective 
# libraries in ${DEVLIB}. The HYPRE library is located in ${LIBHYPRE} 
# with the relevant header files in ${INCLHYPRE}.


# Please also note that you need to extend LD_LIBRARY_PATH
# in order to execute the programme.


export LIBHYPRE=/opt/hypre-2.8.0b/lib
export INCLHYPRE=/opt/hypre-2.8.0b/include

echo "Compiling source files:"
for f in src/*.cpp
do
echo "  -> ${f}"
  mpiCC  -c -O3 -fopenmp -I${DEVINCL} -I${INCLHYPRE} $f  
done
echo "Compilation completed."

echo "Linking objects:"
mpiCC -L${DEVLIB} -L${LIBHYPRE} -O3 -fopenmp -o Impec2D *.o -lblas -llapack -lHYPRE -lconfig++
echo "Linking completed!"

rm *.o                                                                       
echo "Cleaning up completed."
