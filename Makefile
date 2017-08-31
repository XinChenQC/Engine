objects = engine.o mod_data.o basisint.o tool_io.o guess.o SCF.o solveKS.o tool_mat.o grid_gen.o gto_eval.o Lebedev.o DFT.o lib_XC.o force.o
FORT90  = gfortran
tagCOMP = -c -free
tagCOMPF77 = -c


all:	Engine
	


basisint.o : basisint.f90 
	$(FORT90) $(tagCOMP) basisint.f90
mod_data.o : mod_data.f90
	$(FORT90) $(tagCOMP) mod_data.f90
engine.o : engine.f90 mod_data.o
	$(FORT90) $(tagCOMP) engine.f90
tool_io.o : tool_io.f90
	$(FORT90) $(tagCOMP) tool_io.f90
guess.o : guess.f90
	$(FORT90) $(tagCOMP) guess.f90
SCF.o   : SCF.f90
	$(FORT90) $(tagCOMP) SCF.f90
solveKS.o : solveKS.f90
	$(FORT90) $(tagCOMP) solveKS.f90
tool_mat.o: tool_mat.f90
	$(FORT90) $(tagCOMP) tool_mat.f90
gto_eval.o: gto_eval.f90
	$(FORT90) $(tagCOMP) gto_eval.f90
grid_gen.o: grid_gen.f90
	$(FORT90) $(tagCOMP) grid_gen.f90
Lebedev.o: Lebedev.F
	$(FORT90) $(tagCOMPF77) Lebedev.F
DFT.o: DFT.f90
	$(FORT90) $(tagCOMPF77) DFT.f90
lib_XC.o : lib_XC.F
	$(FORT90) $(tagCOMPF77) lib_XC.F
force.o:force.f90
	$(FORT90) $(tagCOMPF77) force.f90


Engine : $(objects)
#	$(FORT90) -fno-backtrace -o  Engine  $(objects) libcint.a   -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

	$(FORT90) -o  Engine -debug -fbacktrace $(objects) libcint.a -lblas -llapack 
.PHONY : clean
clean :
	rm  Engine  $(objects) *.mod 
