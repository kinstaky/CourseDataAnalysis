GXX = g++
OBJS = main.o Tracking.o trackingBase.o
SRC = main.cpp Tracking.cpp trackingBase.cpp

ROOTCFLAGS = $(shell root-config --cflags)
INCLUDE = -I ../include
CFLAGS = -Wall $(ROOTCFLAGS) $(INCLUDE)

ROOTLIBS = $(shell root-config --libs)
LIBS = $(ROOTLIBS)
LDFLAGS = $(LIBS)


all:
	make Tracking

Tracking: $(OBJS)
	$(GXX) -o $@ $^ $(LDFLAGS)
$(OBJS):%.o:%.cpp
	$(GXX) $(CFLAGS) -c $<

clean:
	rm *.o Tracking
