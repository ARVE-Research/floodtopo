# makefile

FC = gfortran
FCFLAGS = -ffree-form -ffree-line-length-none -ftree-vectorize -Wall

# FCFLAGS = -g  -O0 -ffree-line-length-none -fcheck=all -fno-check-array-temporaries -ffpe-trap=invalid,zero,overflow,underflow -g -fbacktrace -Wall -pedantic

# use the command "nf-config --all" to find the location of your netCDF installation
# and enter the path next to " --prefix    ->" on the line below

# netcdf = /home/public/easybuild/software/netCDF-Fortran/4.6.1-gompi-2023a
netcdf = /usr/local

# should not need to modify anything below this line

# ---------------------------------------------

NC_LIB = $(netcdf)/lib
NC_INC = $(netcdf)/include

CPPFLAGS = -I$(NC_INC)
LDFLAGS  = -L$(NC_LIB)
LIBS     = -lnetcdff

OBJS = utilitymod.o netcdf_createmod.o floodmod.o floodtopo.o

.SUFFIXES: .o .f90 .f .mod

%.o : %.c
	$(CC) $(CFLAGS) -c -o $(*F).o $(CPPFLAGS) $<

%.o : %.f
	$(FC) $(FCFLAGS) -c -o $(*F).o $(CPPFLAGS) $<

%.o : %.f90
	$(FC) $(FCFLAGS) -c -o $(*F).o $(CPPFLAGS) $<

all::	floodtopo

floodtopo:	$(OBJS)
	$(FC) $(FCFLAGS) -o floodtopo $(OBJS) $(LDFLAGS) $(LIBS)

clean::	
	-rm floodtopo *.mod *.o
