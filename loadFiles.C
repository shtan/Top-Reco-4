#include "TROOT.h"

void loadFiles(TString dir, int whichLoop, int maxLoops)
{
    std::cout<<"here"<<std::endl;
  //gROOT->LoadMacro("neutrinoSolutions.cxx+g");
  gROOT->LoadMacro("WDaughterEllipseCalculator.cxx+g");
  gROOT->LoadMacro("topSystemChiSquare.cxx+g");
  gROOT->LoadMacro("leptonicTopSystemChiSquare.cxx+g");
  gROOT->LoadMacro("hadronicTopSystemChiSquare.cxx+g");
  gROOT->LoadMacro("lightJetChiSquareMinimumSolver.cxx+g");
  gROOT->LoadMacro("topEventMinimizer.cxx+g");
  gROOT->LoadMacro("topReconstructionFromLHE.C+g");
  std::cout<<"there"<<std::endl;
  topReconstructionFromLHE t;
  std::cout<<"and"<<std::endl;
  t.Loop(dir,whichLoop,maxLoops);
  std::cout<<"everywhere"<<std::endl;
}
