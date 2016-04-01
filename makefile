CXX = $(shell root-config --cxx)
CXXFLAGS = $(shell root-config --cflags) -std=c++11 -Wall -Wextra -pedantic -O2 \
-Wshadow -Wunused-variable -Werror=sign-compare \
-Wunused-but-set-variable -Werror=return-type -Werror=missing-braces \
-Werror=delete-non-virtual-dtor  -fPIC \
$(INCLUDE_RULES)
# needs to be added: -Werror=maybe-uninitialized
INCLUDE_RULES = -I $(shell root-config --incdir) -I . -I inc/
LD = $(shell root-config --ld)
LDFLAGS = $(shell root-config --ldflags)
LDLIBS =  $(shell root-config --glibs) -lMinuit2 -lMathMore -l GenVector

OBJECTS = src/WDaughterEllipseCalculator.o src/lightJetChiSquareMinimumSolver.o \
	  src/topEventMinimizer.o src/topSystemChiSquare.o \
	  src/topReconstructionFromLHE_core.o \
	  src/topReconstructionFromLHE_diagnostics.o \
	  src/topReconstructionFromLHE_methods.o 

COMPILE = $(CXX) $(CXXFLAGS) -c
LINK = $(LD) $(LDFLAGS)
LINKEND = $(OBJECTS) $(LDLIBS)


# add your executable name here
all: topReconstructionFromLHE converter

topReconstructionFromLHE: $(OBJECTS) topReconstructionFromLHE.o
	$(LINK) -o topReconstructionFromLHE topReconstructionFromLHE.o $(LINKEND)

converter: $(OBJECTS) converter.o
	$(LINK) -o converter converter.o $(LINKEND)

clean:
	rm -f topReconstructionFromLHE.o converter.o; \
	rm -f $(OBJECTS)

.PHONY: clean

.SUFFIXES: .cpp .cc .cxx .c

.cxx.o: inc/commonstruct.h
	$(CXX) -c $(CXXFLAGS) -o $@ $<
