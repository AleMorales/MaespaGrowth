##
## Makefile for building a Maespa, could be nicer...
## Martin De Kauwe, 05/12/2013
##

PROG =	maespa.out

SRCS =	growth_module.f90 default_conditions.f90 switches.f90 getmet.f90 maindeclarations.f90 \
        initialize.f90 inout.f90  maespa_growth.f90 maespa.f90 maestcom.f90 metcom.f90 physiol.f90 \
        radn.f90 unstor.f90 utils.f90 watbal.f90

OBJS =	growth_module.o default_conditions.o switches.o getmet.o maindeclarations.o initialize.o \
        inout.o  maespa_growth.o maespa.o maestcom.o metcom.o physiol.o radn.o \
	unstor.o utils.o watbal.o

LIBS =	

INCLS = 

# In Windows we have to use fullpath (even though gfortran is recognized by the command line, make cannot find it!)
F90 = D:\\UserData\\moral005\\Rtools\\gcc-4.6.3\\bin\\gfortran.exe

FFLAGS = -g -fbounds-check -finit-local-zero -Wuninitialized -ftrapv -ffree-form -ffree-line-length-none -O3

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(FFLAGS) -o $@ $(OBJS) $(LIBS) $(INCLS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(FFLAGS) -c $<
	
default_conditions.o: switches.o
getmet.o: maestcom.o metcom.o switches.o
maindeclarations.o: maestcom.o
initialize.o: maindeclarations.o
inout.o: maestcom.o switches.o
growth_module.o: maestcom.o
maespa_growth.o: initialize.o growth_module.o
maespa.o: maestcom.o metcom.o switches.o maindeclarations.o
physiol.o: maestcom.o metcom.o
radn2.o: maestcom.o
unstor.o: maestcom.o
utils.o: maestcom.o
watbal.o: maestcom.o metcom.o
