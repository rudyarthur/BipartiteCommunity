#!/bin/bash
CXX=g++
DIRSRC= ./
CXXFLAGS= -DDEBUG -O5 -Wall -Wno-unused-variable
EXEC=test_network

all: test_network test_weightednetwork test_bipartitenetwork test_projector

test_network : $(DIRSRC)test_network.o 
	$(CXX) -o $@ $^ $(CXXFLAGS)

test_weightednetwork : $(DIRSRC)test_weightednetwork.o 
	$(CXX) -o $@ $^ $(CXXFLAGS)

test_bipartitenetwork : $(DIRSRC)test_bipartitenetwork.o 
	$(CXX) -o $@ $^ $(CXXFLAGS)

test_projector : $(DIRSRC)test_projector.o 
	$(CXX) -o $@ $^ $(CXXFLAGS)
					
##########################################
# Generic rules
##########################################

%.o: %.cpp %.h
	$(CXX) -o  $@ -c $< $(CXXFLAGS)

%.o: %.cpp
	$(CXX) -o $@ -c $< $(CXXFLAGS)

clean:
	rm -f $(DIRSRC)*.o test_network test_weightednetwork test_bipartitenetwork test_projector

	
mrproper: clean
	rm -f *~ $(EXEC)
