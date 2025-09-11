all: pikoe1

pikoe1: pikoe1.f90
	gfortran -O3 $< -o $@
	rm -rf *.mod

.PHONY: clean distclean

clean: 
	rm -rf *.mod

distclean:
	rm ./pikoe1 ./sample*/*.outlist ./sample*/*.dat

