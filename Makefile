CFITSIO = $(FITSIOROOT)
CPP = g++
CC = gcc
CFLAGS = -Wall -I$(CFITSIO) $(shell root-config --cflags)
LIBS = -L$(CFITSIO) -lcfitsio -lm $(shell root-config --libs)
GLIBS = 
GLIBS += 
OBJECTS = fits2root.o 
HEADERS = globalConstants.h

ALL : fits2root.exe
	@echo "Listo!"

fits2root.exe : $(OBJECTS)
	$(CPP) $(OBJECTS) -o fits2root.exe $(LIBS) $(GLIBS) $(CFLAGS)

fits2root.o : fits2root.cc $(HEADERS)
	$(CPP) -c fits2root.cc -o fits2root.o $(CFLAGS)

clean:
	rm -f *~ *.o *.exe
