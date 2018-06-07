# On jupiter, first run module load mvapich2-2.2b-gcc-5.3.0-o4of6w7
CPP := mpicxx
CPPFLAGS := -O3
SRCDIR := ./src
OBJDIR := ./obj
OBJECTS := $(OBJDIR)/util.o $(OBJDIR)/processLC.o $(OBJDIR)/main.o 

# Change trunk loc to match current system...
# Jupiter:
#TRUNK_LOC := /homes/jphollowed/code/lc_interpolation/jupiter
#INCLUDES := -I $(SRCDIR) -I /homes/jphollowed/code/lc_interpolation/genericio
# Cooley:
TRUNK_LOC := /home/hollowed/repos/lc_interpolation/cooley
INCLUDES := -I $(SRCDIR) -I /home/hollowed/repos/lc_interpolation/genericio

LIBS := $(TRUNK_LOC)/mpi/lib/libGenericIOMPI.a

#linking

halo_lightcone_cutout : $(OBJECTS) $(LIBS)
	$(CPP) $(CPPFLAGS) $(OBJECTS) -L$(TRUNK_LOC)/mpi/lib -lGenericIOMPI -L$(TRUNK_LOC)/frontend/lib -lGenericIO -o halo_lc_cutout -fopenmp

#compilation

$(OBJDIR)/util.o: $(SRCDIR)/util.cpp $(SRCDIR)/util.h $(LIBS)
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $(SRCDIR)/util.cpp -o $(OBJDIR)/util.o -fopenmp

$(OBJDIR)/processLC.o: $(SRCDIR)/processLC.cpp $(SRCDIR)/processLC.h $(LIBS)
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $(SRCDIR)/processLC.cpp -o $(OBJDIR)/processLC.o -fopenmp

$(OBJDIR)/main.o: $(SRCDIR)/main.cpp $(LIBS)
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $(SRCDIR)/main.cpp -o $(OBJDIR)/main.o -fopenmp

clean:
	rm $(OBJDIR)/*.o halo_lc_cutout
