GXX = g++
OBJS = main.o Adssd.o AdssdBase.o
SRC = main.cpp Adssd.cpp AdssdBase.cpp

ROOTCFLAGS = $(shell root-config --cflags)
INCLUDE = -I ../include
CFLAGS = -Wall -O3 $(ROOTCFLAGS) $(INCLUDE)

ROOTLIBS = $(shell root-config --libs) -lSpectrum
LIBS = $(ROOTLIBS)
LDFLAGS = $(LIBS)


all:
	make Adssd;
	make RandomOffset;

Adssd: $(OBJS)
	$(GXX) -o $@ $^ $(LDFLAGS)
$(OBJS):%.o:%.cpp
	$(GXX) $(CFLAGS) -c $<

RandomOffset: RandomOffset.o
	$(GXX) -o $@ $^ $(LDFLAGS)
RandomOffset.o:%.o:%.cpp
	$(GXX) $(CFLAGS) -c $<

Hit: hit.o
	$(GXX) -o $@ $^ $(LDFLAGS)
hit.o:%.o:%.cpp
	$(GXX) $(CFLAGS) -c $<

clean:
	rm *.o Adssd
