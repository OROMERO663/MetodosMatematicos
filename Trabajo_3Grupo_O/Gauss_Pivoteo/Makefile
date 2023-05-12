# FORTRAN-90 compiler:
FC = gfortran
# Remove comand:
RM = rm
# Linker:
LD = gfortran
# Compiler Flags:
FFLAGS = -O3 -fbounds-check

.SUFFIXES: .o .f90

# Path to modules folder:
VPATH = ../modules
# List of modules object files (space separated):
OBJS = sistemas_lineales.o main.o

main.exe : $(OBJS)
	$(LD)  -o $@ $(OBJS)

.f90.o:
	$(FC)  -c $^ $(FFLAGS) -o $@

clean:
	$(RM) $(OBJS) *.mod
	
