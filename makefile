CXX = $(shell root-config --cxx)
CPPFLAGS = -isystem$(shell root-config --incdir) -I inc
CXXFLAGS = -Wall -Wextra -pedantic -O2 -Wshadow -fPIC $(shell root-config --cflags)
LD = $(shell root-config --ld)
LDFLAGS = $(shell root-config --ldflags) -lGenVector
LDLIBS =  $(shell root-config --glibs) -lMinuit2 -lMathMore

VPATH = inc:src:obj

OBJECTS = obj/WDaughterEllipseCalculator.o obj/hadronicTopSystemChiSquare.o \
	  obj/leptonicTopSystemChiSquare.o obj/lightJetChiSquareMinimumSolver.o \
	  obj/topEventMinimizer.o obj/topSystemChiSquare.o

COMPILE = $(CXX) $(CXXFLAGS) $(CPPFLAGS) -c
LINK = $(LD) $(LDFLAGS)
LINKEND = $(OBJECTS) $(LDLIBS)

# add your executable name here
all : topReconstructionFromLHE converter

topReconstructionFromLHE : topReconstructionFromLHE.o $(OBJECTS)
	$(LINK) -o topReconstructionFromLHE topReconstructionFromLHE.o $(LINKEND)
topReconstructionFromLHE.o: topReconstructionFromLHE.C topReconstructionFromLHE.h topEventMinimizer.h
	$(COMPILE) topReconstructionFromLHE.C

converter : converter.o $(OBJECTS)
	$(LINK) -o converter converter.o $(LINKEND)
converter.o: converter.C converter.h 
	$(COMPILE) converter.C

clean:
	-rm -f topReconstructionFromLHE converter obj/*.o *.o

obj/WDaughterEllipseCalculator.o : WDaughterEllipseCalculator.cxx WDaughterEllipseCalculator.h
	$(COMPILE) src/WDaughterEllipseCalculator.cxx -o obj/WDaughterEllipseCalculator.o

obj/hadronicTopSystemChiSquare.o : hadronicTopSystemChiSquare.cxx hadronicTopSystemChiSquare.h
	$(COMPILE) src/hadronicTopSystemChiSquare.cxx -o obj/hadronicTopSystemChiSquare.o

obj/leptonicTopSystemChiSquare.o : leptonicTopSystemChiSquare.cxx leptonicTopSystemChiSquare.h
	$(COMPILE) src/leptonicTopSystemChiSquare.cxx -o obj/leptonicTopSystemChiSquare.o

obj/lightJetChiSquareMinimumSolver.o : lightJetChiSquareMinimumSolver.cxx lightJetChiSquareMinimumSolver.h
	$(COMPILE) src/lightJetChiSquareMinimumSolver.cxx -o obj/lightJetChiSquareMinimumSolver.o

obj/topEventMinimizer.o : topEventMinimizer.cxx topEventMinimizer.h
	$(COMPILE) src/topEventMinimizer.cxx -o obj/topEventMinimizer.o

obj/topSystemChiSquare.o : topSystemChiSquare.cxx topSystemChiSquare.h
	$(COMPILE) src/topSystemChiSquare.cxx -o obj/topSystemChiSquare.o

.PHONY : clean
