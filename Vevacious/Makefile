#
# Makefile for Vevacious.
#

SHELL=/bin/sh

CC=g++
LHPCDIR=./LHPC/
LHPCINCLUDEDIR=$(LHPCDIR)/include/
LHPCLIBDIR=$(LHPCDIR)/lib/
CCFLAGS=-Wall -O3 -I./include/ -I$(LHPCINCLUDEDIR)
CCOBJECTFLAGS=-c
CCEXEFLAGS=
LDFLAGS= -L$(LHPCLIBDIR) -lLHPC
SOURCES=source/PotentialMinimizer.cpp \
 source/SarahInterpreter.cpp \
 source/SarahSlhaConverter.cpp \
 source/TadpoleSolver.cpp \
 source/Vevacious.cpp \
 source/VevaciousRunner.cpp \
 source/VevRenamer.cpp
OBJECTS=$(SOURCES:.cpp=.o)
BINDIR=./bin/
EXEFILE=$(BINDIR)/Vevacious.exe

all: $(OBJECTS) $(EXEFILE)

.PHONY: all clean

$(OBJECTS): %.o: %.cpp
	$(CC) $(CCFLAGS) $(CCOBJECTFLAGS) $< -o $@ $(LDFLAGS)
	
$(BINDIR):
	mkdir $@
	
$(EXEFILE): $(BINDIR) $(OBJECTS)
	$(CC) $(CCFLAGS) $(CCEXEFLAGS) -o $@ $(OBJECTS) $(LDFLAGS)

clean:
	\rm -f $(OBJECTS) $(EXEFILE)

deleteobjects:
	\rm -f $(OBJECTS)