#
# Makefile for Vevacious.
#     Author: Ben O'Leary (benjamin.oleary@gmail.com)
#      Copyright 2012 Ben O'Leary
#
#      This file is part of Vevacious, released under the GNU General Public
#      License. Please see the accompanying README.Vevacious.txt file for a
#      full list of files, brief documentation on how to use these classes, and
#      further details on the license.
#

SHELL=/bin/sh

CC=g++
AR=ar
ARFLAGS=rcs
LHPCDIR=./LHPC/
LHPCLIBNAME=LHPC
LHPCINCLUDEDIR=$(LHPCDIR)/include/
LHPCLIBDIR=$(LHPCDIR)/lib/
CCFLAGS=-Wall -O3 -I./include/ -I$(LHPCINCLUDEDIR)
CCOBJECTFLAGS=-c
CCEXEFLAGS=
SOURCES=source/PotentialMinimizer.cpp \
 source/SarahInterpreter.cpp \
 source/SarahSlhaConverter.cpp \
 source/TadpoleSolver.cpp \
 source/VevaciousRunner.cpp \
 source/VevRenamer.cpp
OBJECTS=$(SOURCES:.cpp=.o)
BINDIR=./bin/
LIBNAME=Vevacious
LIBDIR=./lib/
LDFLAGS= -L$(LIBDIR) -l$(LIBNAME)
LIBRARYFILE=$(LIBDIR)/lib$(LIBNAME).a
LDFLAGS= -L$(LIBDIR) -l$(LIBNAME) -L$(LHPCLIBDIR) -l$(LHPCLIBNAME)
SINGLEMAINCPP=source/Vevacious.cpp
SINGLEEXEFILE=$(BINDIR)/Vevacious.exe
BATCHMAINCPP=source/VevaciousBatchRunner.cpp
BATCHEXEFILE=$(BINDIR)/VevaciousBatch.exe
EXES=$(SINGLEEXEFILE) $(BATCHEXEFILE)

all: $(OBJECTS) $(LIBRARYFILE) $(EXES)

.PHONY: all clean

$(OBJECTS): %.o: %.cpp
	$(CC) $(CCFLAGS) $(CCOBJECTFLAGS) $< -o $@
	
$(LIBDIR):
	mkdir $@

$(LIBRARYFILE): $(LIBDIR) $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $(OBJECTS)

$(BINDIR):
	mkdir $@
	
$(SINGLEEXEFILE): $(BINDIR) $(LIBRARYFILE)
	$(CC) $(CCFLAGS) $(CCEXEFLAGS) -o $@ $(SINGLEMAINCPP) $(LDFLAGS)
	
$(BATCHEXEFILE): $(BINDIR) $(LIBRARYFILE)
	$(CC) $(CCFLAGS) $(CCEXEFLAGS) -o $@ $(BATCHMAINCPP) $(LDFLAGS)

clean:
	\rm -f $(OBJECTS) $(LIBRARYFILE) $(EXEFILE)

deleteobjects:
	\rm -f $(OBJECTS)