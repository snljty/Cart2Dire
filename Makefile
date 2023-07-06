SHELL = cmd
FC := gfortran
FLINKER := $(FC)
DEBUG_LEVEL := 
ifeq ("$(DEBUG_LEVEL)", "")
	OPTS := -s
	OPTLV := -O2
else
	OPTS := -g
	OPTLV := -O0
endif
LAPACKROOT := C:\lapack-3.11.0
LIBPATH := -L $(LAPACKROOT)\lib
LIB := -l lapack -l blas

.PHONY: all
all: Cart2Dire.exe

Cart2Dire.exe: Cart2Dire.obj convert_cell_coordinates.obj read_write_vasp_poscar.obj read_write_qe_pw_inp.obj
	@echo Linking $@ against $^ and LAPACK and BLAS ...
	$(FLINKER) -o $@ $^ -static $(LIBPATH) $(LIB) $(OPTS)

%.obj: %.f90
	@echo Compiling $@ ...
	$(FC) -o $@ -c $< $(OPTLV) $(OPTS)

.PHONY: clean
clean:
	-del /q *.obj 2> NUL

.PHONY: veryclean
veryclean: clean
	-del /q *.exe 2> NUL

