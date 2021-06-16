GXX = g++
OBJS =

ROOTCFLAGS = $(shell root-config --cflags)
INCLUDE = -I ../include
CFLAGS = -Wall -O3 $(ROOTCFLAGS) $(INCLUDE)

ROOTLIBS = $(shell root-config --libs) -lSpectrum
LIBS = $(ROOTLIBS)
LDFLAGS = $(LIBS)

all:
	make Gamma;

Gamma: GammaGamma.o
	$(GXX) $(CFLAGS) -o $@ $^ $(LDFLAGS)
GammaGamma.o: GammaGamma.cpp
	$(GXX) $(CFLAGS) -c $^