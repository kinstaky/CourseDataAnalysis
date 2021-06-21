GXX = g++
OBJS = dssd.o

DSSDMAIN = dssdMain.o
DSSDSRC = dssdMain.cpp

ROOTCFLAGS = $(shell root-config --cflags)
INCLUDE = -I ../include
CFLAGS = -Wall -O3 $(ROOTCFLAGS) $(INCLUDE)

ROOTLIBS = $(shell root-config --libs) -lSpectrum
LIBS = $(ROOTLIBS)
LDFLAGS = $(LIBS)


all:
	make DSSD;
	make AddMWPC;
	make DecayRa;

DSSD: $(DSSDMAIN) $(OBJS)
	$(GXX) -o $@ $^ $(LDFLAGS)
$(DSSDMAIN): $(DSSDSRC)
	$(GXX) -o $@ $(CFLAGS) -c $^
$(OBJS):%.o:%.cpp
	$(GXX) $(CFLAGS) -c $<
AddMWPC: addMWPC.o
	$(GXX) -o $@ $^ $(LDFLAGS)
addMWPC.o: addMWPC.cpp
	$(GXX) $(CFLAGS) -c $^
DecayRa: decayRa.o
	$(GXX) -o $@ $^ $(LDFLAGS)
decayRa.o: decayRa.cpp
	$(GXX) $(CFLAGS) -c $^

clean:
	rm *.o DSSD AddMWPC decayRa
