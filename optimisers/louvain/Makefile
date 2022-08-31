#!/bin/bash
CXX=g++
DIRSRC= ./
CXXFLAGS= -I../ -DDEBUG -O5 -Wall -Wno-unused-variable
EXEC=test_louvain

all: test_louvain test_invlouvain test_projectedlouvain test_aggregate test_cplouvain test_dualprojection test_bilouvain

test_louvain : $(DIRSRC)test_louvain.o
	$(CXX) -o $@ $^ $(CXXFLAGS)

test_invlouvain : $(DIRSRC)test_invlouvain.o
	$(CXX) -o $@ $^ $(CXXFLAGS)

test_projectedlouvain : $(DIRSRC)test_projectedlouvain.o
	$(CXX) -o $@ $^ $(CXXFLAGS)

test_aggregate : $(DIRSRC)test_aggregate.o
	$(CXX) -o $@ $^ $(CXXFLAGS)
	
test_cplouvain : $(DIRSRC)test_cplouvain.o
	$(CXX) -o $@ $^ $(CXXFLAGS)
	
test_dualprojection : $(DIRSRC)test_dualprojection.o
	$(CXX) -o $@ $^ $(CXXFLAGS)

test_bilouvain : $(DIRSRC)test_bilouvain.o
	$(CXX) -o $@ $^ $(CXXFLAGS)
			
##########################################
# Generic rules
##########################################

%.o: %.cpp %.h
	$(CXX) -o  $@ -c $< $(CXXFLAGS)

%.o: %.cpp
	$(CXX) -o $@ -c $< $(CXXFLAGS)

clean:
	rm -f $(DIRSRC)*.o test_louvain test_invlouvain test_projectedlouvain test_aggregate test_cplouvain test_dualprojection test_bilouvain

	
mrproper: clean
	rm -f *~ $(EXEC)
