#!/bin/bash
CXX=g++
DIRSRC= ./
CXXFLAGS= -DDEBUG -O5 -Wall -Wno-unused-variable
EXEC=test

all: test_synclabelprop test_asynclabelprop

test_synclabelprop : $(DIRSRC)test_synclabelprop.o
	$(CXX) -o $@ $^ $(CXXFLAGS)

test_asynclabelprop : $(DIRSRC)test_asynclabelprop.o
	$(CXX) -o $@ $^ $(CXXFLAGS)
		
##########################################
# Generic rules
##########################################

%.o: %.cpp %.h
	$(CXX) -o  $@ -c $< $(CXXFLAGS)

%.o: %.cpp
	$(CXX) -o $@ -c $< $(CXXFLAGS)

clean:
	rm -f $(DIRSRC)*.o test_synclabelprop test_asynclabelprop

	
mrproper: clean
	rm -f *~ $(EXEC)
