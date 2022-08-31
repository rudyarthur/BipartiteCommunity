#!/bin/bash
CXX=g++
DIRSRC= ./
CXXFLAGS= -DDEBUG -O5 -Wall -Wno-unused-variable
EXEC=test

all: test_modularity test_barbermodularity test_projectedmodularity test_cpmodularity

test_modularity : $(DIRSRC)test_modularity.o
	$(CXX) -o $@ $^ $(CXXFLAGS)
test_barbermodularity : $(DIRSRC)test_barbermodularity.o
	$(CXX) -o $@ $^ $(CXXFLAGS)
test_projectedmodularity : $(DIRSRC)test_projectedmodularity.o
	$(CXX) -o $@ $^ $(CXXFLAGS)
test_cpmodularity : $(DIRSRC)test_cpmodularity.o
	$(CXX) -o $@ $^ $(CXXFLAGS)
				
##########################################
# Generic rules
##########################################

%.o: %.cpp %.h
	$(CXX) -o  $@ -c $< $(CXXFLAGS)

%.o: %.cpp
	$(CXX) -o $@ -c $< $(CXXFLAGS)

clean:
	rm -f $(DIRSRC)*.o test_modularity test_barbermodularity test_projectedmodularity test_cpmodularity

	
mrproper: clean
	rm -f *~ $(EXEC)
