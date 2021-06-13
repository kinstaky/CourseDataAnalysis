GXX = g++
OBJS = main.o DSSD.o
SRC = main.cpp DSSD.cpp

SOBJS = SimpleNorm.o DSSD.o
SSRC = SimpleNorm.cpp DSSD.cpp

SHOBJS = Shrinker.o Shrink.o

ROOTCFLAGS = $(shell root-config --cflags)
INCLUDE = -I ../include
CFLAGS = -Wall -O3 $(ROOTCFLAGS) $(INCLUDE)

ROOTLIBS = $(shell root-config --libs) -lSpectrum
LIBS = $(ROOTLIBS)
LDFLAGS = $(LIBS)


all:
	make DSSD;
	make Shrinker;
#	make SimpleNorm;

DSSD: $(OBJS)
	$(GXX) -o $@ $^ $(LDFLAGS)
$(OBJS):%.o:%.cpp
	$(GXX) $(CFLAGS) -c $<

SimpleNorm: $(SOBJS)
	$(GXX) -o $@ $^ $(LDFLAGS)
$(SOBJS):%.o:%.cpp
	$(GXX) $(CFLAGS) -c $<


Shrinker: $(SHOBJS)
	$(GXX) -o $@ $^ $(LDFLAGS)
$(SHOBJS):%.o:%.cpp
	$(GXX) $(CFLAGS) -c $<

clean:
	rm *.o DSSD SimpleNorm Shrinker
