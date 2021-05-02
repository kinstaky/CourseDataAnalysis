GXX = g++
OBJS = main.o Simulation.o
SRC = main.cpp Simulation.cpp

ROOTCFLAGS = $(shell root-config --cflags)
INCLUDE = -I ../../include
CFLAGS = -Wall -O3 $(ROOTCFLAGS) $(INCLUDE)

ROOTLIBS = $(shell root-config --libs)
LIBS = $(ROOTLIBS)
LDFLAGS = $(LIBS)


all:
	make Simulation
Simulation: $(OBJS)
	$(GXX) -o $@ $^ $(LDFLAGS)
$(OBJS):%.o:%.cpp
	$(GXX) $(CFLAGS) -c $<

clean:
	rm *.o Simulation
