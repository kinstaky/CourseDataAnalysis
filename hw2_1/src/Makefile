GXX = g++
OBJS = main.o Tracking.o trackingBase.o
OBJSPG = main_pg.o Tracking_pg.o trackingBase_pg.o
SRC = main.cpp Tracking.cpp trackingBase.cpp

ROOTCFLAGS = $(shell root-config --cflags)
INCLUDE = -I ../include
CFLAGS = -Wall -O3 $(ROOTCFLAGS) $(INCLUDE)

ROOTLIBS = $(shell root-config --libs)
LIBS = $(ROOTLIBS)
LDFLAGS = $(LIBS)


all:
	make Tracking

Tracking: $(OBJS)
	$(GXX) -o $@ $^ $(LDFLAGS)
$(OBJS):%.o:%.cpp
	$(GXX) $(CFLAGS) -c $<

profile:
	make Tracking_profile

Tracking_profile: $(OBJSPG)
	$(GXX) -pg -o $@ $^ $(LDFLAGS)
$(OBJSPG):%_pg.o:%.cpp
	$(GXX) $(CFLAGS) -pg -c -lc_p -o $@ $<

clean:
	rm *.o Tracking Tracking_profile
