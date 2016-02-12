//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Dec 11 11:41:46 2015 by ROOT version 6.02/05
// from TTree tree/tree
// found on file: ./output_files/output_0.root
//////////////////////////////////////////////////////////

#ifndef Plotter_h
#define Plotter_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class Plotter
{
  public:
    TTree *fChain; //! pointer to the analyzed TTree or TChain
    Int_t fCurrent; //! current Tree number in a TChain

    // Fixed size dimensions of array or collections stored in the TTree if any.

    // Declaration of leaf types
    Int_t eventNumber;
    Int_t innerMinStatus;
    Int_t outerMinStatus;
    Double_t outerMinEdm;
    Double_t totalChi2;
    Double_t topSystemChi2;
    Double_t topMassChi2;
    Double_t hadronicChi2;
    Double_t nonTopChi2;
    Double_t bJet1PtDelta;
    Double_t bJet1PhiDelta;
    Double_t bJet1EtaDelta;
    Double_t bJet2PtDelta;
    Double_t bJet2PhiDelta;
    Double_t bJet2EtaDelta;
    Double_t W1Daughter1PtDelta;
    Double_t W1Daughter1PhiDelta;
    Double_t W1Daughter1EtaDelta;
    Double_t W2Daughter1PtDelta;
    Double_t W2Daughter1PhiDelta;
    Double_t W2Daughter1EtaDelta;
    Double_t W1Daughter2PtDelta;
    Double_t W1Daughter2PhiDelta;
    Double_t W1Daughter2EtaDelta;
    Double_t W2Daughter2PtDelta;
    Double_t W2Daughter2PhiDelta;
    Double_t W2Daughter2EtaDelta;
    Double_t top1MassDelta;
    Double_t W1MassDelta;
    Double_t top2MassDelta;
    Double_t W2MassDelta;
    Double_t lightJet1PxDelta;
    Double_t lightJet1PyDelta;
    Double_t lightJet2PxDelta;
    Double_t lightJet2PyDelta;
    Double_t top1TruePt;
    Double_t top1TrueEta;
    Double_t top1TruePhi;
    Double_t top1TrueE;
    Double_t top1TrueMass;
    Double_t W1TruePt;
    Double_t W1TrueEta;
    Double_t W1TruePhi;
    Double_t W1TrueE;
    Double_t W1TrueMass;
    Double_t bJet1TruePt;
    Double_t bJet1TrueEta;
    Double_t bJet1TruePhi;
    Double_t bJet1TrueE;
    Double_t W1Daughter1TruePt;
    Double_t W1Daughter1TrueEta;
    Double_t W1Daughter1TruePhi;
    Double_t W1Daughter1TrueE;
    Double_t W1Daughter2TruePt;
    Double_t W1Daughter2TrueEta;
    Double_t W1Daughter2TruePhi;
    Double_t W1Daughter2TrueE;
    Double_t lightJet1TruePx;
    Double_t lightJet1TruePy;
    Double_t lightJet1TruePz;
    Double_t lightJet1TrueE;
    Double_t lightJet1TruePt;
    Double_t top2TruePt;
    Double_t top2TrueEta;
    Double_t top2TruePhi;
    Double_t top2TrueE;
    Double_t top2TrueMass;
    Double_t W2TruePt;
    Double_t W2TrueEta;
    Double_t W2TruePhi;
    Double_t W2TrueE;
    Double_t W2TrueMass;
    Double_t bJet2TruePt;
    Double_t bJet2TrueEta;
    Double_t bJet2TruePhi;
    Double_t bJet2TrueE;
    Double_t W2Daughter1TruePt;
    Double_t W2Daughter1TrueEta;
    Double_t W2Daughter1TruePhi;
    Double_t W2Daughter1TrueE;
    Double_t W2Daughter2TruePt;
    Double_t W2Daughter2TrueEta;
    Double_t W2Daughter2TruePhi;
    Double_t W2Daughter2TrueE;
    Double_t lightJet2TruePx;
    Double_t lightJet2TruePy;
    Double_t lightJet2TruePz;
    Double_t lightJet2TrueE;
    Double_t lightJet2TruePt;
    Double_t top1SmearedPt;
    Double_t top1SmearedEta;
    Double_t top1SmearedPhi;
    Double_t top1SmearedE;
    Double_t top1SmearedMass;
    Double_t W1SmearedPt;
    Double_t W1SmearedEta;
    Double_t W1SmearedPhi;
    Double_t W1SmearedE;
    Double_t W1SmearedMass;
    Double_t bJet1SmearedPt;
    Double_t bJet1SmearedEta;
    Double_t bJet1SmearedPhi;
    Double_t bJet1SmearedE;
    Double_t W1Daughter1SmearedPt;
    Double_t W1Daughter1SmearedEta;
    Double_t W1Daughter1SmearedPhi;
    Double_t W1Daughter1SmearedE;
    Double_t W1Daughter2SmearedPt;
    Double_t W1Daughter2SmearedEta;
    Double_t W1Daughter2SmearedPhi;
    Double_t W1Daughter2SmearedE;
    Double_t lightJet1SmearedPx;
    Double_t lightJet1SmearedPy;
    Double_t lightJet1SmearedPz;
    Double_t lightJet1SmearedE;
    Double_t lightJet1SmearedPt;
    Double_t top2SmearedPt;
    Double_t top2SmearedEta;
    Double_t top2SmearedPhi;
    Double_t top2SmearedE;
    Double_t top2SmearedMass;
    Double_t W2SmearedPt;
    Double_t W2SmearedEta;
    Double_t W2SmearedPhi;
    Double_t W2SmearedE;
    Double_t W2SmearedMass;
    Double_t bJet2SmearedPt;
    Double_t bJet2SmearedEta;
    Double_t bJet2SmearedPhi;
    Double_t bJet2SmearedE;
    Double_t W2Daughter1SmearedPt;
    Double_t W2Daughter1SmearedEta;
    Double_t W2Daughter1SmearedPhi;
    Double_t W2Daughter1SmearedE;
    Double_t W2Daughter2SmearedPt;
    Double_t W2Daughter2SmearedEta;
    Double_t W2Daughter2SmearedPhi;
    Double_t W2Daughter2SmearedE;
    Double_t lightJet2SmearedPx;
    Double_t lightJet2SmearedPy;
    Double_t lightJet2SmearedPz;
    Double_t lightJet2SmearedE;
    Double_t lightJet2SmearedPt;
    Double_t top1BestPt;
    Double_t top1BestEta;
    Double_t top1BestPhi;
    Double_t top1BestE;
    Double_t top1BestMass;
    Double_t W1BestPt;
    Double_t W1BestEta;
    Double_t W1BestPhi;
    Double_t W1BestE;
    Double_t W1BestMass;
    Double_t bJet1BestPt;
    Double_t bJet1BestEta;
    Double_t bJet1BestPhi;
    Double_t bJet1BestE;
    Double_t W1Daughter1BestPt;
    Double_t W1Daughter1BestEta;
    Double_t W1Daughter1BestPhi;
    Double_t W1Daughter1BestE;
    Double_t W1Daughter2BestPt;
    Double_t W1Daughter2BestEta;
    Double_t W1Daughter2BestPhi;
    Double_t W1Daughter2BestE;
    Double_t lightJet1BestPx;
    Double_t lightJet1BestPy;
    Double_t lightJet1BestPz;
    Double_t lightJet1BestE;
    Double_t lightJet1BestPt;
    Double_t top2BestPt;
    Double_t top2BestEta;
    Double_t top2BestPhi;
    Double_t top2BestE;
    Double_t top2BestMass;
    Double_t W2BestPt;
    Double_t W2BestEta;
    Double_t W2BestPhi;
    Double_t W2BestE;
    Double_t W2BestMass;
    Double_t bJet2BestPt;
    Double_t bJet2BestEta;
    Double_t bJet2BestPhi;
    Double_t bJet2BestE;
    Double_t W2Daughter1BestPt;
    Double_t W2Daughter1BestEta;
    Double_t W2Daughter1BestPhi;
    Double_t W2Daughter1BestE;
    Double_t W2Daughter2BestPt;
    Double_t W2Daughter2BestEta;
    Double_t W2Daughter2BestPhi;
    Double_t W2Daughter2BestE;
    Double_t lightJet2BestPx;
    Double_t lightJet2BestPy;
    Double_t lightJet2BestPz;
    Double_t lightJet2BestE;
    Double_t lightJet2BestPt;
    Double_t top1DeltaPtTrueSmeared;
    Double_t top1DeltaRTrueSmeared;
    Double_t top1DeltaMTrueSmeared;
    Double_t W1DeltaPtTrueSmeared;
    Double_t W1DeltaRTrueSmeared;
    Double_t W1DeltaMTrueSmeared;
    Double_t bJet1DeltaPtTrueSmeared;
    Double_t bJet1DeltaRTrueSmeared;
    Double_t W1Daughter1DeltaPtTrueSmeared;
    Double_t W1Daughter1DeltaRTrueSmeared;
    Double_t W1Daughter2DeltaPtTrueSmeared;
    Double_t W1Daughter2DeltaRTrueSmeared;
    Double_t lightJet1DeltaPxTrueSmeared;
    Double_t lightJet1DeltaPyTrueSmeared;
    Double_t lightJet1DeltaRTrueSmeared;
    Double_t top2DeltaPtTrueSmeared;
    Double_t top2DeltaRTrueSmeared;
    Double_t top2DeltaMTrueSmeared;
    Double_t W2DeltaPtTrueSmeared;
    Double_t W2DeltaRTrueSmeared;
    Double_t W2DeltaMTrueSmeared;
    Double_t bJet2DeltaPtTrueSmeared;
    Double_t bJet2DeltaRTrueSmeared;
    Double_t W2Daughter1DeltaPtTrueSmeared;
    Double_t W2Daughter1DeltaRTrueSmeared;
    Double_t W2Daughter2DeltaPtTrueSmeared;
    Double_t W2Daughter2DeltaRTrueSmeared;
    Double_t lightJet2DeltaPxTrueSmeared;
    Double_t lightJet2DeltaPyTrueSmeared;
    Double_t lightJet2DeltaRTrueSmeared;
    Double_t top1DeltaPtTrueBest;
    Double_t top1DeltaRTrueBest;
    Double_t top1DeltaMTrueBest;
    Double_t W1DeltaPtTrueBest;
    Double_t W1DeltaRTrueBest;
    Double_t W1DeltaMTrueBest;
    Double_t bJet1DeltaPtTrueBest;
    Double_t bJet1DeltaRTrueBest;
    Double_t W1Daughter1DeltaPtTrueBest;
    Double_t W1Daughter1DeltaRTrueBest;
    Double_t W1Daughter2DeltaPtTrueBest;
    Double_t W1Daughter2DeltaRTrueBest;
    Double_t lightJet1DeltaPxTrueBest;
    Double_t lightJet1DeltaPyTrueBest;
    Double_t lightJet1DeltaRTrueBest;
    Double_t top2DeltaPtTrueBest;
    Double_t top2DeltaRTrueBest;
    Double_t top2DeltaMTrueBest;
    Double_t W2DeltaPtTrueBest;
    Double_t W2DeltaRTrueBest;
    Double_t W2DeltaMTrueBest;
    Double_t bJet2DeltaPtTrueBest;
    Double_t bJet2DeltaRTrueBest;
    Double_t W2Daughter1DeltaPtTrueBest;
    Double_t W2Daughter1DeltaRTrueBest;
    Double_t W2Daughter2DeltaPtTrueBest;
    Double_t W2Daughter2DeltaRTrueBest;
    Double_t lightJet2DeltaPxTrueBest;
    Double_t lightJet2DeltaPyTrueBest;
    Double_t lightJet2DeltaRTrueBest;
    Double_t top1DeltaPtSmearedBest;
    Double_t top1DeltaRSmearedBest;
    Double_t top1DeltaMSmearedBest;
    Double_t W1DeltaPtSmearedBest;
    Double_t W1DeltaRSmearedBest;
    Double_t W1DeltaMSmearedBest;
    Double_t bJet1DeltaPtSmearedBest;
    Double_t bJet1DeltaRSmearedBest;
    Double_t bJet1DeltaMSmearedBest;
    Double_t W1Daughter1DeltaPtSmearedBest;
    Double_t W1Daughter1DeltaRSmearedBest;
    Double_t W1Daughter1DeltaMSmearedBest;
    Double_t W1Daughter2DeltaPtSmearedBest;
    Double_t W1Daughter2DeltaRSmearedBest;
    Double_t W1Daughter2DeltaMSmearedBest;
    Double_t lightJet1DeltaPxSmearedBest;
    Double_t lightJet1DeltaPySmearedBest;
    Double_t lightJet1DeltaRSmearedBest;
    Double_t top2DeltaPtSmearedBest;
    Double_t top2DeltaRSmearedBest;
    Double_t top2DeltaMSmearedBest;
    Double_t W2DeltaPtSmearedBest;
    Double_t W2DeltaRSmearedBest;
    Double_t W2DeltaMSmearedBest;
    Double_t bJet2DeltaPtSmearedBest;
    Double_t bJet2DeltaRSmearedBest;
    Double_t bJet2DeltaMSmearedBest;
    Double_t W2Daughter1DeltaPtSmearedBest;
    Double_t W2Daughter1DeltaRSmearedBest;
    Double_t W2Daughter1DeltaMSmearedBest;
    Double_t W2Daughter2DeltaPtSmearedBest;
    Double_t W2Daughter2DeltaRSmearedBest;
    Double_t W2Daughter2DeltaMSmearedBest;
    Double_t lightJet2DeltaPxSmearedBest;
    Double_t lightJet2DeltaPySmearedBest;
    Double_t lightJet2DeltaRSmearedBest;

    // List of branches
    TBranch *b_eventNumber;                   //!
    TBranch *b_innerMinStatus;                //!
    TBranch *b_outerMinStatus;                //!
    TBranch *b_outerMinEdm;                   //!
    TBranch *b_totalChi2;                     //!
    TBranch *b_topSystemChi2;                 //!
    TBranch *b_topMassChi2;                   //!
    TBranch *b_hadronicChi2;                  //!
    TBranch *b_nonTopChi2;                    //!
    TBranch *b_bJet1PtDelta;                  //!
    TBranch *b_bJet1PhiDelta;                 //!
    TBranch *b_bJet1EtaDelta;                 //!
    TBranch *b_bJet2PtDelta;                  //!
    TBranch *b_bJet2PhiDelta;                 //!
    TBranch *b_bJet2EtaDelta;                 //!
    TBranch *b_W1Daughter1PtDelta;            //!
    TBranch *b_W1Daughter1PhiDelta;           //!
    TBranch *b_W1Daughter1EtaDelta;           //!
    TBranch *b_W2Daughter1PtDelta;            //!
    TBranch *b_W2Daughter1PhiDelta;           //!
    TBranch *b_W2Daughter1EtaDelta;           //!
    TBranch *b_W1Daughter2PtDelta;            //!
    TBranch *b_W1Daughter2PhiDelta;           //!
    TBranch *b_W1Daughter2EtaDelta;           //!
    TBranch *b_W2Daughter2PtDelta;            //!
    TBranch *b_W2Daughter2PhiDelta;           //!
    TBranch *b_W2Daughter2EtaDelta;           //!
    TBranch *b_top1MassDelta;                 //!
    TBranch *b_W1MassDelta;                   //!
    TBranch *b_top2MassDelta;                 //!
    TBranch *b_W2MassDelta;                   //!
    TBranch *b_lightJet1PxDelta;              //!
    TBranch *b_lightJet1PyDelta;              //!
    TBranch *b_lightJet2PxDelta;              //!
    TBranch *b_lightJet2PyDelta;              //!
    TBranch *b_top1TruePt;                    //!
    TBranch *b_top1TrueEta;                   //!
    TBranch *b_top1TruePhi;                   //!
    TBranch *b_top1TrueE;                     //!
    TBranch *b_top1TrueMass;                  //!
    TBranch *b_W1TruePt;                      //!
    TBranch *b_W1TrueEta;                     //!
    TBranch *b_W1TruePhi;                     //!
    TBranch *b_W1TrueE;                       //!
    TBranch *b_W1TrueMass;                    //!
    TBranch *b_bJet1TruePt;                   //!
    TBranch *b_bJet1TrueEta;                  //!
    TBranch *b_bJet1TruePhi;                  //!
    TBranch *b_bJet1TrueE;                    //!
    TBranch *b_W1Daughter1TruePt;             //!
    TBranch *b_W1Daughter1TrueEta;            //!
    TBranch *b_W1Daughter1TruePhi;            //!
    TBranch *b_W1Daughter1TrueE;              //!
    TBranch *b_W1Daughter2TruePt;             //!
    TBranch *b_W1Daughter2TrueEta;            //!
    TBranch *b_W1Daughter2TruePhi;            //!
    TBranch *b_W1Daughter2TrueE;              //!
    TBranch *b_lightJet1TruePx;               //!
    TBranch *b_lightJet1TruePy;               //!
    TBranch *b_lightJet1TruePz;               //!
    TBranch *b_lightJet1TrueE;                //!
    TBranch *b_lightJet1TruePt;               //!
    TBranch *b_top2TruePt;                    //!
    TBranch *b_top2TrueEta;                   //!
    TBranch *b_top2TruePhi;                   //!
    TBranch *b_top2TrueE;                     //!
    TBranch *b_top2TrueMass;                  //!
    TBranch *b_W2TruePt;                      //!
    TBranch *b_W2TrueEta;                     //!
    TBranch *b_W2TruePhi;                     //!
    TBranch *b_W2TrueE;                       //!
    TBranch *b_W2TrueMass;                    //!
    TBranch *b_bJet2TruePt;                   //!
    TBranch *b_bJet2TrueEta;                  //!
    TBranch *b_bJet2TruePhi;                  //!
    TBranch *b_bJet2TrueE;                    //!
    TBranch *b_W2Daughter1TruePt;             //!
    TBranch *b_W2Daughter1TrueEta;            //!
    TBranch *b_W2Daughter1TruePhi;            //!
    TBranch *b_W2Daughter1TrueE;              //!
    TBranch *b_W2Daughter2TruePt;             //!
    TBranch *b_W2Daughter2TrueEta;            //!
    TBranch *b_W2Daughter2TruePhi;            //!
    TBranch *b_W2Daughter2TrueE;              //!
    TBranch *b_lightJet2TruePx;               //!
    TBranch *b_lightJet2TruePy;               //!
    TBranch *b_lightJet2TruePz;               //!
    TBranch *b_lightJet2TrueE;                //!
    TBranch *b_lightJet2TruePt;               //!
    TBranch *b_top1SmearedPt;                 //!
    TBranch *b_top1SmearedEta;                //!
    TBranch *b_top1SmearedPhi;                //!
    TBranch *b_top1SmearedE;                  //!
    TBranch *b_top1SmearedMass;               //!
    TBranch *b_W1SmearedPt;                   //!
    TBranch *b_W1SmearedEta;                  //!
    TBranch *b_W1SmearedPhi;                  //!
    TBranch *b_W1SmearedE;                    //!
    TBranch *b_W1SmearedMass;                 //!
    TBranch *b_bJet1SmearedPt;                //!
    TBranch *b_bJet1SmearedEta;               //!
    TBranch *b_bJet1SmearedPhi;               //!
    TBranch *b_bJet1SmearedE;                 //!
    TBranch *b_W1Daughter1SmearedPt;          //!
    TBranch *b_W1Daughter1SmearedEta;         //!
    TBranch *b_W1Daughter1SmearedPhi;         //!
    TBranch *b_W1Daughter1SmearedE;           //!
    TBranch *b_W1Daughter2SmearedPt;          //!
    TBranch *b_W1Daughter2SmearedEta;         //!
    TBranch *b_W1Daughter2SmearedPhi;         //!
    TBranch *b_W1Daughter2SmearedE;           //!
    TBranch *b_lightJet1SmearedPx;            //!
    TBranch *b_lightJet1SmearedPy;            //!
    TBranch *b_lightJet1SmearedPz;            //!
    TBranch *b_lightJet1SmearedE;             //!
    TBranch *b_lightJet1SmearedPt;            //!
    TBranch *b_top2SmearedPt;                 //!
    TBranch *b_top2SmearedEta;                //!
    TBranch *b_top2SmearedPhi;                //!
    TBranch *b_top2SmearedE;                  //!
    TBranch *b_top2SmearedMass;               //!
    TBranch *b_W2SmearedPt;                   //!
    TBranch *b_W2SmearedEta;                  //!
    TBranch *b_W2SmearedPhi;                  //!
    TBranch *b_W2SmearedE;                    //!
    TBranch *b_W2SmearedMass;                 //!
    TBranch *b_bJet2SmearedPt;                //!
    TBranch *b_bJet2SmearedEta;               //!
    TBranch *b_bJet2SmearedPhi;               //!
    TBranch *b_bJet2SmearedE;                 //!
    TBranch *b_W2Daughter1SmearedPt;          //!
    TBranch *b_W2Daughter1SmearedEta;         //!
    TBranch *b_W2Daughter1SmearedPhi;         //!
    TBranch *b_W2Daughter1SmearedE;           //!
    TBranch *b_W2Daughter2SmearedPt;          //!
    TBranch *b_W2Daughter2SmearedEta;         //!
    TBranch *b_W2Daughter2SmearedPhi;         //!
    TBranch *b_W2Daughter2SmearedE;           //!
    TBranch *b_lightJet2SmearedPx;            //!
    TBranch *b_lightJet2SmearedPy;            //!
    TBranch *b_lightJet2SmearedPz;            //!
    TBranch *b_lightJet2SmearedE;             //!
    TBranch *b_lightJet2SmearedPt;            //!
    TBranch *b_top1BestPt;                    //!
    TBranch *b_top1BestEta;                   //!
    TBranch *b_top1BestPhi;                   //!
    TBranch *b_top1BestE;                     //!
    TBranch *b_top1BestMass;                  //!
    TBranch *b_W1BestPt;                      //!
    TBranch *b_W1BestEta;                     //!
    TBranch *b_W1BestPhi;                     //!
    TBranch *b_W1BestE;                       //!
    TBranch *b_W1BestMass;                    //!
    TBranch *b_bJet1BestPt;                   //!
    TBranch *b_bJet1BestEta;                  //!
    TBranch *b_bJet1BestPhi;                  //!
    TBranch *b_bJet1BestE;                    //!
    TBranch *b_W1Daughter1BestPt;             //!
    TBranch *b_W1Daughter1BestEta;            //!
    TBranch *b_W1Daughter1BestPhi;            //!
    TBranch *b_W1Daughter1BestE;              //!
    TBranch *b_W1Daughter2BestPt;             //!
    TBranch *b_W1Daughter2BestEta;            //!
    TBranch *b_W1Daughter2BestPhi;            //!
    TBranch *b_W1Daughter2BestE;              //!
    TBranch *b_lightJet1BestPx;               //!
    TBranch *b_lightJet1BestPy;               //!
    TBranch *b_lightJet1BestPz;               //!
    TBranch *b_lightJet1BestE;                //!
    TBranch *b_lightJet1BestPt;               //!
    TBranch *b_top2BestPt;                    //!
    TBranch *b_top2BestEta;                   //!
    TBranch *b_top2BestPhi;                   //!
    TBranch *b_top2BestE;                     //!
    TBranch *b_top2BestMass;                  //!
    TBranch *b_W2BestPt;                      //!
    TBranch *b_W2BestEta;                     //!
    TBranch *b_W2BestPhi;                     //!
    TBranch *b_W2BestE;                       //!
    TBranch *b_W2BestMass;                    //!
    TBranch *b_bJet2BestPt;                   //!
    TBranch *b_bJet2BestEta;                  //!
    TBranch *b_bJet2BestPhi;                  //!
    TBranch *b_bJet2BestE;                    //!
    TBranch *b_W2Daughter1BestPt;             //!
    TBranch *b_W2Daughter1BestEta;            //!
    TBranch *b_W2Daughter1BestPhi;            //!
    TBranch *b_W2Daughter1BestE;              //!
    TBranch *b_W2Daughter2BestPt;             //!
    TBranch *b_W2Daughter2BestEta;            //!
    TBranch *b_W2Daughter2BestPhi;            //!
    TBranch *b_W2Daughter2BestE;              //!
    TBranch *b_lightJet2BestPx;               //!
    TBranch *b_lightJet2BestPy;               //!
    TBranch *b_lightJet2BestPz;               //!
    TBranch *b_lightJet2BestE;                //!
    TBranch *b_lightJet2BestPt;               //!
    TBranch *b_top1DeltaPtTrueSmeared;        //!
    TBranch *b_top1DeltaRTrueSmeared;         //!
    TBranch *b_top1DeltaMTrueSmeared;         //!
    TBranch *b_W1DeltaPtTrueSmeared;          //!
    TBranch *b_W1DeltaRTrueSmeared;           //!
    TBranch *b_W1DeltaMTrueSmeared;           //!
    TBranch *b_bJet1DeltaPtTrueSmeared;       //!
    TBranch *b_bJet1DeltaRTrueSmeared;        //!
    TBranch *b_W1Daughter1DeltaPtTrueSmeared; //!
    TBranch *b_W1Daughter1DeltaRTrueSmeared;  //!
    TBranch *b_W1Daughter2DeltaPtTrueSmeared; //!
    TBranch *b_W1Daughter2DeltaRTrueSmeared;  //!
    TBranch *b_lightJet1DeltaPxTrueSmeared;   //!
    TBranch *b_lightJet1DeltaPyTrueSmeared;   //!
    TBranch *b_lightJet1DeltaRTrueSmeared;    //!
    TBranch *b_top2DeltaPtTrueSmeared;        //!
    TBranch *b_top2DeltaRTrueSmeared;         //!
    TBranch *b_top2DeltaMTrueSmeared;         //!
    TBranch *b_W2DeltaPtTrueSmeared;          //!
    TBranch *b_W2DeltaRTrueSmeared;           //!
    TBranch *b_W2DeltaMTrueSmeared;           //!
    TBranch *b_bJet2DeltaPtTrueSmeared;       //!
    TBranch *b_bJet2DeltaRTrueSmeared;        //!
    TBranch *b_W2Daughter1DeltaPtTrueSmeared; //!
    TBranch *b_W2Daughter1DeltaRTrueSmeared;  //!
    TBranch *b_W2Daughter2DeltaPtTrueSmeared; //!
    TBranch *b_W2Daughter2DeltaRTrueSmeared;  //!
    TBranch *b_lightJet2DeltaPxTrueSmeared;   //!
    TBranch *b_lightJet2DeltaPyTrueSmeared;   //!
    TBranch *b_lightJet2DeltaRTrueSmeared;    //!
    TBranch *b_top1DeltaPtTrueBest;           //!
    TBranch *b_top1DeltaRTrueBest;            //!
    TBranch *b_top1DeltaMTrueBest;            //!
    TBranch *b_W1DeltaPtTrueBest;             //!
    TBranch *b_W1DeltaRTrueBest;              //!
    TBranch *b_W1DeltaMTrueBest;              //!
    TBranch *b_bJet1DeltaPtTrueBest;          //!
    TBranch *b_bJet1DeltaRTrueBest;           //!
    TBranch *b_W1Daughter1DeltaPtTrueBest;    //!
    TBranch *b_W1Daughter1DeltaRTrueBest;     //!
    TBranch *b_W1Daughter2DeltaPtTrueBest;    //!
    TBranch *b_W1Daughter2DeltaRTrueBest;     //!
    TBranch *b_lightJet1DeltaPxTrueBest;      //!
    TBranch *b_lightJet1DeltaPyTrueBest;      //!
    TBranch *b_lightJet1DeltaRTrueBest;       //!
    TBranch *b_top2DeltaPtTrueBest;           //!
    TBranch *b_top2DeltaRTrueBest;            //!
    TBranch *b_top2DeltaMTrueBest;            //!
    TBranch *b_W2DeltaPtTrueBest;             //!
    TBranch *b_W2DeltaRTrueBest;              //!
    TBranch *b_W2DeltaMTrueBest;              //!
    TBranch *b_bJet2DeltaPtTrueBest;          //!
    TBranch *b_bJet2DeltaRTrueBest;           //!
    TBranch *b_W2Daughter1DeltaPtTrueBest;    //!
    TBranch *b_W2Daughter1DeltaRTrueBest;     //!
    TBranch *b_W2Daughter2DeltaPtTrueBest;    //!
    TBranch *b_W2Daughter2DeltaRTrueBest;     //!
    TBranch *b_lightJet2DeltaPxTrueBest;      //!
    TBranch *b_lightJet2DeltaPyTrueBest;      //!
    TBranch *b_lightJet2DeltaRTrueBest;       //!
    TBranch *b_top1DeltaPtSmearedBest;        //!
    TBranch *b_top1DeltaRSmearedBest;         //!
    TBranch *b_top1DeltaMSmearedBest;         //!
    TBranch *b_W1DeltaPtSmearedBest;          //!
    TBranch *b_W1DeltaRSmearedBest;           //!
    TBranch *b_W1DeltaMSmearedBest;           //!
    TBranch *b_bJet1DeltaPtSmearedBest;       //!
    TBranch *b_bJet1DeltaRSmearedBest;        //!
    TBranch *b_bJet1DeltaMSmearedBest;        //!
    TBranch *b_W1Daughter1DeltaPtSmearedBest; //!
    TBranch *b_W1Daughter1DeltaRSmearedBest;  //!
    TBranch *b_W1Daughter1DeltaMSmearedBest;  //!
    TBranch *b_W1Daughter2DeltaPtSmearedBest; //!
    TBranch *b_W1Daughter2DeltaRSmearedBest;  //!
    TBranch *b_W1Daughter2DeltaMSmearedBest;  //!
    TBranch *b_lightJet1DeltaPxSmearedBest;   //!
    TBranch *b_lightJet1DeltaPySmearedBest;   //!
    TBranch *b_lightJet1DeltaRSmearedBest;    //!
    TBranch *b_top2DeltaPtSmearedBest;        //!
    TBranch *b_top2DeltaRSmearedBest;         //!
    TBranch *b_top2DeltaMSmearedBest;         //!
    TBranch *b_W2DeltaPtSmearedBest;          //!
    TBranch *b_W2DeltaRSmearedBest;           //!
    TBranch *b_W2DeltaMSmearedBest;           //!
    TBranch *b_bJet2DeltaPtSmearedBest;       //!
    TBranch *b_bJet2DeltaRSmearedBest;        //!
    TBranch *b_bJet2DeltaMSmearedBest;        //!
    TBranch *b_W2Daughter1DeltaPtSmearedBest; //!
    TBranch *b_W2Daughter1DeltaRSmearedBest;  //!
    TBranch *b_W2Daughter1DeltaMSmearedBest;  //!
    TBranch *b_W2Daughter2DeltaPtSmearedBest; //!
    TBranch *b_W2Daughter2DeltaRSmearedBest;  //!
    TBranch *b_W2Daughter2DeltaMSmearedBest;  //!
    TBranch *b_lightJet2DeltaPxSmearedBest;   //!
    TBranch *b_lightJet2DeltaPySmearedBest;   //!
    TBranch *b_lightJet2DeltaRSmearedBest;    //!

    Plotter(TTree *tree = 0);
    virtual ~Plotter();
    virtual Int_t Cut(Long64_t entry);
    virtual Int_t GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void Init(TTree *tree);
    virtual void Loop();
    virtual Bool_t Notify();
    virtual void Show(Long64_t entry = -1);
};

#endif

#ifdef Plotter_cxx
Plotter::Plotter(TTree *tree) : fChain(0)
{
    // if parameter tree is not specified (or zero), connect the file
    // used to generate this class and read the Tree.
    if (tree == 0) {
        TFile *f = (TFile *)gROOT->GetListOfFiles()->FindObject(
            "./output_files/output_0.root");
        if (!f || !f->IsOpen()) {
            f = new TFile("./output_files/output_0.root");
        }
        f->GetObject("tree", tree);
    }
    Init(tree);
}

Plotter::~Plotter()
{
    if (!fChain)
        return;
    delete fChain->GetCurrentFile();
}

Int_t Plotter::GetEntry(Long64_t entry)
{
    // Read contents of entry.
    if (!fChain)
        return 0;
    return fChain->GetEntry(entry);
}
Long64_t Plotter::LoadTree(Long64_t entry)
{
    // Set the environment to read one entry
    if (!fChain)
        return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0)
        return centry;
    if (fChain->GetTreeNumber() != fCurrent) {
        fCurrent = fChain->GetTreeNumber();
        Notify();
    }
    return centry;
}

void Plotter::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set branch addresses and branch pointers
    if (!tree)
        return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
    fChain->SetBranchAddress("innerMinStatus", &innerMinStatus,
                             &b_innerMinStatus);
    fChain->SetBranchAddress("outerMinStatus", &outerMinStatus,
                             &b_outerMinStatus);
    fChain->SetBranchAddress("outerMinEdm", &outerMinEdm, &b_outerMinEdm);
    fChain->SetBranchAddress("totalChi2", &totalChi2, &b_totalChi2);
    fChain->SetBranchAddress("topSystemChi2", &topSystemChi2, &b_topSystemChi2);
    fChain->SetBranchAddress("topMassChi2", &topMassChi2, &b_topMassChi2);
    fChain->SetBranchAddress("hadronicChi2", &hadronicChi2, &b_hadronicChi2);
    fChain->SetBranchAddress("nonTopChi2", &nonTopChi2, &b_nonTopChi2);
    fChain->SetBranchAddress("bJet1PtDelta", &bJet1PtDelta, &b_bJet1PtDelta);
    fChain->SetBranchAddress("bJet1PhiDelta", &bJet1PhiDelta, &b_bJet1PhiDelta);
    fChain->SetBranchAddress("bJet1EtaDelta", &bJet1EtaDelta, &b_bJet1EtaDelta);
    fChain->SetBranchAddress("bJet2PtDelta", &bJet2PtDelta, &b_bJet2PtDelta);
    fChain->SetBranchAddress("bJet2PhiDelta", &bJet2PhiDelta, &b_bJet2PhiDelta);
    fChain->SetBranchAddress("bJet2EtaDelta", &bJet2EtaDelta, &b_bJet2EtaDelta);
    fChain->SetBranchAddress("W1Daughter1PtDelta", &W1Daughter1PtDelta,
                             &b_W1Daughter1PtDelta);
    fChain->SetBranchAddress("W1Daughter1PhiDelta", &W1Daughter1PhiDelta,
                             &b_W1Daughter1PhiDelta);
    fChain->SetBranchAddress("W1Daughter1EtaDelta", &W1Daughter1EtaDelta,
                             &b_W1Daughter1EtaDelta);
    fChain->SetBranchAddress("W2Daughter1PtDelta", &W2Daughter1PtDelta,
                             &b_W2Daughter1PtDelta);
    fChain->SetBranchAddress("W2Daughter1PhiDelta", &W2Daughter1PhiDelta,
                             &b_W2Daughter1PhiDelta);
    fChain->SetBranchAddress("W2Daughter1EtaDelta", &W2Daughter1EtaDelta,
                             &b_W2Daughter1EtaDelta);
    fChain->SetBranchAddress("W1Daughter2PtDelta", &W1Daughter2PtDelta,
                             &b_W1Daughter2PtDelta);
    fChain->SetBranchAddress("W1Daughter2PhiDelta", &W1Daughter2PhiDelta,
                             &b_W1Daughter2PhiDelta);
    fChain->SetBranchAddress("W1Daughter2EtaDelta", &W1Daughter2EtaDelta,
                             &b_W1Daughter2EtaDelta);
    fChain->SetBranchAddress("W2Daughter2PtDelta", &W2Daughter2PtDelta,
                             &b_W2Daughter2PtDelta);
    fChain->SetBranchAddress("W2Daughter2PhiDelta", &W2Daughter2PhiDelta,
                             &b_W2Daughter2PhiDelta);
    fChain->SetBranchAddress("W2Daughter2EtaDelta", &W2Daughter2EtaDelta,
                             &b_W2Daughter2EtaDelta);
    fChain->SetBranchAddress("top1MassDelta", &top1MassDelta, &b_top1MassDelta);
    fChain->SetBranchAddress("W1MassDelta", &W1MassDelta, &b_W1MassDelta);
    fChain->SetBranchAddress("top2MassDelta", &top2MassDelta, &b_top2MassDelta);
    fChain->SetBranchAddress("W2MassDelta", &W2MassDelta, &b_W2MassDelta);
    fChain->SetBranchAddress("lightJet1PxDelta", &lightJet1PxDelta,
                             &b_lightJet1PxDelta);
    fChain->SetBranchAddress("lightJet1PyDelta", &lightJet1PyDelta,
                             &b_lightJet1PyDelta);
    fChain->SetBranchAddress("lightJet2PxDelta", &lightJet2PxDelta,
                             &b_lightJet2PxDelta);
    fChain->SetBranchAddress("lightJet2PyDelta", &lightJet2PyDelta,
                             &b_lightJet2PyDelta);
    fChain->SetBranchAddress("top1TruePt", &top1TruePt, &b_top1TruePt);
    fChain->SetBranchAddress("top1TrueEta", &top1TrueEta, &b_top1TrueEta);
    fChain->SetBranchAddress("top1TruePhi", &top1TruePhi, &b_top1TruePhi);
    fChain->SetBranchAddress("top1TrueE", &top1TrueE, &b_top1TrueE);
    fChain->SetBranchAddress("top1TrueMass", &top1TrueMass, &b_top1TrueMass);
    fChain->SetBranchAddress("W1TruePt", &W1TruePt, &b_W1TruePt);
    fChain->SetBranchAddress("W1TrueEta", &W1TrueEta, &b_W1TrueEta);
    fChain->SetBranchAddress("W1TruePhi", &W1TruePhi, &b_W1TruePhi);
    fChain->SetBranchAddress("W1TrueE", &W1TrueE, &b_W1TrueE);
    fChain->SetBranchAddress("W1TrueMass", &W1TrueMass, &b_W1TrueMass);
    fChain->SetBranchAddress("bJet1TruePt", &bJet1TruePt, &b_bJet1TruePt);
    fChain->SetBranchAddress("bJet1TrueEta", &bJet1TrueEta, &b_bJet1TrueEta);
    fChain->SetBranchAddress("bJet1TruePhi", &bJet1TruePhi, &b_bJet1TruePhi);
    fChain->SetBranchAddress("bJet1TrueE", &bJet1TrueE, &b_bJet1TrueE);
    fChain->SetBranchAddress("W1Daughter1TruePt", &W1Daughter1TruePt,
                             &b_W1Daughter1TruePt);
    fChain->SetBranchAddress("W1Daughter1TrueEta", &W1Daughter1TrueEta,
                             &b_W1Daughter1TrueEta);
    fChain->SetBranchAddress("W1Daughter1TruePhi", &W1Daughter1TruePhi,
                             &b_W1Daughter1TruePhi);
    fChain->SetBranchAddress("W1Daughter1TrueE", &W1Daughter1TrueE,
                             &b_W1Daughter1TrueE);
    fChain->SetBranchAddress("W1Daughter2TruePt", &W1Daughter2TruePt,
                             &b_W1Daughter2TruePt);
    fChain->SetBranchAddress("W1Daughter2TrueEta", &W1Daughter2TrueEta,
                             &b_W1Daughter2TrueEta);
    fChain->SetBranchAddress("W1Daughter2TruePhi", &W1Daughter2TruePhi,
                             &b_W1Daughter2TruePhi);
    fChain->SetBranchAddress("W1Daughter2TrueE", &W1Daughter2TrueE,
                             &b_W1Daughter2TrueE);
    fChain->SetBranchAddress("lightJet1TruePx", &lightJet1TruePx,
                             &b_lightJet1TruePx);
    fChain->SetBranchAddress("lightJet1TruePy", &lightJet1TruePy,
                             &b_lightJet1TruePy);
    fChain->SetBranchAddress("lightJet1TruePz", &lightJet1TruePz,
                             &b_lightJet1TruePz);
    fChain->SetBranchAddress("lightJet1TrueE", &lightJet1TrueE,
                             &b_lightJet1TrueE);
    fChain->SetBranchAddress("lightJet1TruePt", &lightJet1TruePt,
                             &b_lightJet1TruePt);
    fChain->SetBranchAddress("top2TruePt", &top2TruePt, &b_top2TruePt);
    fChain->SetBranchAddress("top2TrueEta", &top2TrueEta, &b_top2TrueEta);
    fChain->SetBranchAddress("top2TruePhi", &top2TruePhi, &b_top2TruePhi);
    fChain->SetBranchAddress("top2TrueE", &top2TrueE, &b_top2TrueE);
    fChain->SetBranchAddress("top2TrueMass", &top2TrueMass, &b_top2TrueMass);
    fChain->SetBranchAddress("W2TruePt", &W2TruePt, &b_W2TruePt);
    fChain->SetBranchAddress("W2TrueEta", &W2TrueEta, &b_W2TrueEta);
    fChain->SetBranchAddress("W2TruePhi", &W2TruePhi, &b_W2TruePhi);
    fChain->SetBranchAddress("W2TrueE", &W2TrueE, &b_W2TrueE);
    fChain->SetBranchAddress("W2TrueMass", &W2TrueMass, &b_W2TrueMass);
    fChain->SetBranchAddress("bJet2TruePt", &bJet2TruePt, &b_bJet2TruePt);
    fChain->SetBranchAddress("bJet2TrueEta", &bJet2TrueEta, &b_bJet2TrueEta);
    fChain->SetBranchAddress("bJet2TruePhi", &bJet2TruePhi, &b_bJet2TruePhi);
    fChain->SetBranchAddress("bJet2TrueE", &bJet2TrueE, &b_bJet2TrueE);
    fChain->SetBranchAddress("W2Daughter1TruePt", &W2Daughter1TruePt,
                             &b_W2Daughter1TruePt);
    fChain->SetBranchAddress("W2Daughter1TrueEta", &W2Daughter1TrueEta,
                             &b_W2Daughter1TrueEta);
    fChain->SetBranchAddress("W2Daughter1TruePhi", &W2Daughter1TruePhi,
                             &b_W2Daughter1TruePhi);
    fChain->SetBranchAddress("W2Daughter1TrueE", &W2Daughter1TrueE,
                             &b_W2Daughter1TrueE);
    fChain->SetBranchAddress("W2Daughter2TruePt", &W2Daughter2TruePt,
                             &b_W2Daughter2TruePt);
    fChain->SetBranchAddress("W2Daughter2TrueEta", &W2Daughter2TrueEta,
                             &b_W2Daughter2TrueEta);
    fChain->SetBranchAddress("W2Daughter2TruePhi", &W2Daughter2TruePhi,
                             &b_W2Daughter2TruePhi);
    fChain->SetBranchAddress("W2Daughter2TrueE", &W2Daughter2TrueE,
                             &b_W2Daughter2TrueE);
    fChain->SetBranchAddress("lightJet2TruePx", &lightJet2TruePx,
                             &b_lightJet2TruePx);
    fChain->SetBranchAddress("lightJet2TruePy", &lightJet2TruePy,
                             &b_lightJet2TruePy);
    fChain->SetBranchAddress("lightJet2TruePz", &lightJet2TruePz,
                             &b_lightJet2TruePz);
    fChain->SetBranchAddress("lightJet2TrueE", &lightJet2TrueE,
                             &b_lightJet2TrueE);
    fChain->SetBranchAddress("lightJet2TruePt", &lightJet2TruePt,
                             &b_lightJet2TruePt);
    fChain->SetBranchAddress("top1SmearedPt", &top1SmearedPt, &b_top1SmearedPt);
    fChain->SetBranchAddress("top1SmearedEta", &top1SmearedEta,
                             &b_top1SmearedEta);
    fChain->SetBranchAddress("top1SmearedPhi", &top1SmearedPhi,
                             &b_top1SmearedPhi);
    fChain->SetBranchAddress("top1SmearedE", &top1SmearedE, &b_top1SmearedE);
    fChain->SetBranchAddress("top1SmearedMass", &top1SmearedMass,
                             &b_top1SmearedMass);
    fChain->SetBranchAddress("W1SmearedPt", &W1SmearedPt, &b_W1SmearedPt);
    fChain->SetBranchAddress("W1SmearedEta", &W1SmearedEta, &b_W1SmearedEta);
    fChain->SetBranchAddress("W1SmearedPhi", &W1SmearedPhi, &b_W1SmearedPhi);
    fChain->SetBranchAddress("W1SmearedE", &W1SmearedE, &b_W1SmearedE);
    fChain->SetBranchAddress("W1SmearedMass", &W1SmearedMass, &b_W1SmearedMass);
    fChain->SetBranchAddress("bJet1SmearedPt", &bJet1SmearedPt,
                             &b_bJet1SmearedPt);
    fChain->SetBranchAddress("bJet1SmearedEta", &bJet1SmearedEta,
                             &b_bJet1SmearedEta);
    fChain->SetBranchAddress("bJet1SmearedPhi", &bJet1SmearedPhi,
                             &b_bJet1SmearedPhi);
    fChain->SetBranchAddress("bJet1SmearedE", &bJet1SmearedE, &b_bJet1SmearedE);
    fChain->SetBranchAddress("W1Daughter1SmearedPt", &W1Daughter1SmearedPt,
                             &b_W1Daughter1SmearedPt);
    fChain->SetBranchAddress("W1Daughter1SmearedEta", &W1Daughter1SmearedEta,
                             &b_W1Daughter1SmearedEta);
    fChain->SetBranchAddress("W1Daughter1SmearedPhi", &W1Daughter1SmearedPhi,
                             &b_W1Daughter1SmearedPhi);
    fChain->SetBranchAddress("W1Daughter1SmearedE", &W1Daughter1SmearedE,
                             &b_W1Daughter1SmearedE);
    fChain->SetBranchAddress("W1Daughter2SmearedPt", &W1Daughter2SmearedPt,
                             &b_W1Daughter2SmearedPt);
    fChain->SetBranchAddress("W1Daughter2SmearedEta", &W1Daughter2SmearedEta,
                             &b_W1Daughter2SmearedEta);
    fChain->SetBranchAddress("W1Daughter2SmearedPhi", &W1Daughter2SmearedPhi,
                             &b_W1Daughter2SmearedPhi);
    fChain->SetBranchAddress("W1Daughter2SmearedE", &W1Daughter2SmearedE,
                             &b_W1Daughter2SmearedE);
    fChain->SetBranchAddress("lightJet1SmearedPx", &lightJet1SmearedPx,
                             &b_lightJet1SmearedPx);
    fChain->SetBranchAddress("lightJet1SmearedPy", &lightJet1SmearedPy,
                             &b_lightJet1SmearedPy);
    fChain->SetBranchAddress("lightJet1SmearedPz", &lightJet1SmearedPz,
                             &b_lightJet1SmearedPz);
    fChain->SetBranchAddress("lightJet1SmearedE", &lightJet1SmearedE,
                             &b_lightJet1SmearedE);
    fChain->SetBranchAddress("lightJet1SmearedPt", &lightJet1SmearedPt,
                             &b_lightJet1SmearedPt);
    fChain->SetBranchAddress("top2SmearedPt", &top2SmearedPt, &b_top2SmearedPt);
    fChain->SetBranchAddress("top2SmearedEta", &top2SmearedEta,
                             &b_top2SmearedEta);
    fChain->SetBranchAddress("top2SmearedPhi", &top2SmearedPhi,
                             &b_top2SmearedPhi);
    fChain->SetBranchAddress("top2SmearedE", &top2SmearedE, &b_top2SmearedE);
    fChain->SetBranchAddress("top2SmearedMass", &top2SmearedMass,
                             &b_top2SmearedMass);
    fChain->SetBranchAddress("W2SmearedPt", &W2SmearedPt, &b_W2SmearedPt);
    fChain->SetBranchAddress("W2SmearedEta", &W2SmearedEta, &b_W2SmearedEta);
    fChain->SetBranchAddress("W2SmearedPhi", &W2SmearedPhi, &b_W2SmearedPhi);
    fChain->SetBranchAddress("W2SmearedE", &W2SmearedE, &b_W2SmearedE);
    fChain->SetBranchAddress("W2SmearedMass", &W2SmearedMass, &b_W2SmearedMass);
    fChain->SetBranchAddress("bJet2SmearedPt", &bJet2SmearedPt,
                             &b_bJet2SmearedPt);
    fChain->SetBranchAddress("bJet2SmearedEta", &bJet2SmearedEta,
                             &b_bJet2SmearedEta);
    fChain->SetBranchAddress("bJet2SmearedPhi", &bJet2SmearedPhi,
                             &b_bJet2SmearedPhi);
    fChain->SetBranchAddress("bJet2SmearedE", &bJet2SmearedE, &b_bJet2SmearedE);
    fChain->SetBranchAddress("W2Daughter1SmearedPt", &W2Daughter1SmearedPt,
                             &b_W2Daughter1SmearedPt);
    fChain->SetBranchAddress("W2Daughter1SmearedEta", &W2Daughter1SmearedEta,
                             &b_W2Daughter1SmearedEta);
    fChain->SetBranchAddress("W2Daughter1SmearedPhi", &W2Daughter1SmearedPhi,
                             &b_W2Daughter1SmearedPhi);
    fChain->SetBranchAddress("W2Daughter1SmearedE", &W2Daughter1SmearedE,
                             &b_W2Daughter1SmearedE);
    fChain->SetBranchAddress("W2Daughter2SmearedPt", &W2Daughter2SmearedPt,
                             &b_W2Daughter2SmearedPt);
    fChain->SetBranchAddress("W2Daughter2SmearedEta", &W2Daughter2SmearedEta,
                             &b_W2Daughter2SmearedEta);
    fChain->SetBranchAddress("W2Daughter2SmearedPhi", &W2Daughter2SmearedPhi,
                             &b_W2Daughter2SmearedPhi);
    fChain->SetBranchAddress("W2Daughter2SmearedE", &W2Daughter2SmearedE,
                             &b_W2Daughter2SmearedE);
    fChain->SetBranchAddress("lightJet2SmearedPx", &lightJet2SmearedPx,
                             &b_lightJet2SmearedPx);
    fChain->SetBranchAddress("lightJet2SmearedPy", &lightJet2SmearedPy,
                             &b_lightJet2SmearedPy);
    fChain->SetBranchAddress("lightJet2SmearedPz", &lightJet2SmearedPz,
                             &b_lightJet2SmearedPz);
    fChain->SetBranchAddress("lightJet2SmearedE", &lightJet2SmearedE,
                             &b_lightJet2SmearedE);
    fChain->SetBranchAddress("lightJet2SmearedPt", &lightJet2SmearedPt,
                             &b_lightJet2SmearedPt);
    fChain->SetBranchAddress("top1BestPt", &top1BestPt, &b_top1BestPt);
    fChain->SetBranchAddress("top1BestEta", &top1BestEta, &b_top1BestEta);
    fChain->SetBranchAddress("top1BestPhi", &top1BestPhi, &b_top1BestPhi);
    fChain->SetBranchAddress("top1BestE", &top1BestE, &b_top1BestE);
    fChain->SetBranchAddress("top1BestMass", &top1BestMass, &b_top1BestMass);
    fChain->SetBranchAddress("W1BestPt", &W1BestPt, &b_W1BestPt);
    fChain->SetBranchAddress("W1BestEta", &W1BestEta, &b_W1BestEta);
    fChain->SetBranchAddress("W1BestPhi", &W1BestPhi, &b_W1BestPhi);
    fChain->SetBranchAddress("W1BestE", &W1BestE, &b_W1BestE);
    fChain->SetBranchAddress("W1BestMass", &W1BestMass, &b_W1BestMass);
    fChain->SetBranchAddress("bJet1BestPt", &bJet1BestPt, &b_bJet1BestPt);
    fChain->SetBranchAddress("bJet1BestEta", &bJet1BestEta, &b_bJet1BestEta);
    fChain->SetBranchAddress("bJet1BestPhi", &bJet1BestPhi, &b_bJet1BestPhi);
    fChain->SetBranchAddress("bJet1BestE", &bJet1BestE, &b_bJet1BestE);
    fChain->SetBranchAddress("W1Daughter1BestPt", &W1Daughter1BestPt,
                             &b_W1Daughter1BestPt);
    fChain->SetBranchAddress("W1Daughter1BestEta", &W1Daughter1BestEta,
                             &b_W1Daughter1BestEta);
    fChain->SetBranchAddress("W1Daughter1BestPhi", &W1Daughter1BestPhi,
                             &b_W1Daughter1BestPhi);
    fChain->SetBranchAddress("W1Daughter1BestE", &W1Daughter1BestE,
                             &b_W1Daughter1BestE);
    fChain->SetBranchAddress("W1Daughter2BestPt", &W1Daughter2BestPt,
                             &b_W1Daughter2BestPt);
    fChain->SetBranchAddress("W1Daughter2BestEta", &W1Daughter2BestEta,
                             &b_W1Daughter2BestEta);
    fChain->SetBranchAddress("W1Daughter2BestPhi", &W1Daughter2BestPhi,
                             &b_W1Daughter2BestPhi);
    fChain->SetBranchAddress("W1Daughter2BestE", &W1Daughter2BestE,
                             &b_W1Daughter2BestE);
    fChain->SetBranchAddress("lightJet1BestPx", &lightJet1BestPx,
                             &b_lightJet1BestPx);
    fChain->SetBranchAddress("lightJet1BestPy", &lightJet1BestPy,
                             &b_lightJet1BestPy);
    fChain->SetBranchAddress("lightJet1BestPz", &lightJet1BestPz,
                             &b_lightJet1BestPz);
    fChain->SetBranchAddress("lightJet1BestE", &lightJet1BestE,
                             &b_lightJet1BestE);
    fChain->SetBranchAddress("lightJet1BestPt", &lightJet1BestPt,
                             &b_lightJet1BestPt);
    fChain->SetBranchAddress("top2BestPt", &top2BestPt, &b_top2BestPt);
    fChain->SetBranchAddress("top2BestEta", &top2BestEta, &b_top2BestEta);
    fChain->SetBranchAddress("top2BestPhi", &top2BestPhi, &b_top2BestPhi);
    fChain->SetBranchAddress("top2BestE", &top2BestE, &b_top2BestE);
    fChain->SetBranchAddress("top2BestMass", &top2BestMass, &b_top2BestMass);
    fChain->SetBranchAddress("W2BestPt", &W2BestPt, &b_W2BestPt);
    fChain->SetBranchAddress("W2BestEta", &W2BestEta, &b_W2BestEta);
    fChain->SetBranchAddress("W2BestPhi", &W2BestPhi, &b_W2BestPhi);
    fChain->SetBranchAddress("W2BestE", &W2BestE, &b_W2BestE);
    fChain->SetBranchAddress("W2BestMass", &W2BestMass, &b_W2BestMass);
    fChain->SetBranchAddress("bJet2BestPt", &bJet2BestPt, &b_bJet2BestPt);
    fChain->SetBranchAddress("bJet2BestEta", &bJet2BestEta, &b_bJet2BestEta);
    fChain->SetBranchAddress("bJet2BestPhi", &bJet2BestPhi, &b_bJet2BestPhi);
    fChain->SetBranchAddress("bJet2BestE", &bJet2BestE, &b_bJet2BestE);
    fChain->SetBranchAddress("W2Daughter1BestPt", &W2Daughter1BestPt,
                             &b_W2Daughter1BestPt);
    fChain->SetBranchAddress("W2Daughter1BestEta", &W2Daughter1BestEta,
                             &b_W2Daughter1BestEta);
    fChain->SetBranchAddress("W2Daughter1BestPhi", &W2Daughter1BestPhi,
                             &b_W2Daughter1BestPhi);
    fChain->SetBranchAddress("W2Daughter1BestE", &W2Daughter1BestE,
                             &b_W2Daughter1BestE);
    fChain->SetBranchAddress("W2Daughter2BestPt", &W2Daughter2BestPt,
                             &b_W2Daughter2BestPt);
    fChain->SetBranchAddress("W2Daughter2BestEta", &W2Daughter2BestEta,
                             &b_W2Daughter2BestEta);
    fChain->SetBranchAddress("W2Daughter2BestPhi", &W2Daughter2BestPhi,
                             &b_W2Daughter2BestPhi);
    fChain->SetBranchAddress("W2Daughter2BestE", &W2Daughter2BestE,
                             &b_W2Daughter2BestE);
    fChain->SetBranchAddress("lightJet2BestPx", &lightJet2BestPx,
                             &b_lightJet2BestPx);
    fChain->SetBranchAddress("lightJet2BestPy", &lightJet2BestPy,
                             &b_lightJet2BestPy);
    fChain->SetBranchAddress("lightJet2BestPz", &lightJet2BestPz,
                             &b_lightJet2BestPz);
    fChain->SetBranchAddress("lightJet2BestE", &lightJet2BestE,
                             &b_lightJet2BestE);
    fChain->SetBranchAddress("lightJet2BestPt", &lightJet2BestPt,
                             &b_lightJet2BestPt);
    fChain->SetBranchAddress("top1DeltaPtTrueSmeared", &top1DeltaPtTrueSmeared,
                             &b_top1DeltaPtTrueSmeared);
    fChain->SetBranchAddress("top1DeltaRTrueSmeared", &top1DeltaRTrueSmeared,
                             &b_top1DeltaRTrueSmeared);
    fChain->SetBranchAddress("top1DeltaMTrueSmeared", &top1DeltaMTrueSmeared,
                             &b_top1DeltaMTrueSmeared);
    fChain->SetBranchAddress("W1DeltaPtTrueSmeared", &W1DeltaPtTrueSmeared,
                             &b_W1DeltaPtTrueSmeared);
    fChain->SetBranchAddress("W1DeltaRTrueSmeared", &W1DeltaRTrueSmeared,
                             &b_W1DeltaRTrueSmeared);
    fChain->SetBranchAddress("W1DeltaMTrueSmeared", &W1DeltaMTrueSmeared,
                             &b_W1DeltaMTrueSmeared);
    fChain->SetBranchAddress("bJet1DeltaPtTrueSmeared",
                             &bJet1DeltaPtTrueSmeared,
                             &b_bJet1DeltaPtTrueSmeared);
    fChain->SetBranchAddress("bJet1DeltaRTrueSmeared", &bJet1DeltaRTrueSmeared,
                             &b_bJet1DeltaRTrueSmeared);
    fChain->SetBranchAddress("W1Daughter1DeltaPtTrueSmeared",
                             &W1Daughter1DeltaPtTrueSmeared,
                             &b_W1Daughter1DeltaPtTrueSmeared);
    fChain->SetBranchAddress("W1Daughter1DeltaRTrueSmeared",
                             &W1Daughter1DeltaRTrueSmeared,
                             &b_W1Daughter1DeltaRTrueSmeared);
    fChain->SetBranchAddress("W1Daughter2DeltaPtTrueSmeared",
                             &W1Daughter2DeltaPtTrueSmeared,
                             &b_W1Daughter2DeltaPtTrueSmeared);
    fChain->SetBranchAddress("W1Daughter2DeltaRTrueSmeared",
                             &W1Daughter2DeltaRTrueSmeared,
                             &b_W1Daughter2DeltaRTrueSmeared);
    fChain->SetBranchAddress("lightJet1DeltaPxTrueSmeared",
                             &lightJet1DeltaPxTrueSmeared,
                             &b_lightJet1DeltaPxTrueSmeared);
    fChain->SetBranchAddress("lightJet1DeltaPyTrueSmeared",
                             &lightJet1DeltaPyTrueSmeared,
                             &b_lightJet1DeltaPyTrueSmeared);
    fChain->SetBranchAddress("lightJet1DeltaRTrueSmeared",
                             &lightJet1DeltaRTrueSmeared,
                             &b_lightJet1DeltaRTrueSmeared);
    fChain->SetBranchAddress("top2DeltaPtTrueSmeared", &top2DeltaPtTrueSmeared,
                             &b_top2DeltaPtTrueSmeared);
    fChain->SetBranchAddress("top2DeltaRTrueSmeared", &top2DeltaRTrueSmeared,
                             &b_top2DeltaRTrueSmeared);
    fChain->SetBranchAddress("top2DeltaMTrueSmeared", &top2DeltaMTrueSmeared,
                             &b_top2DeltaMTrueSmeared);
    fChain->SetBranchAddress("W2DeltaPtTrueSmeared", &W2DeltaPtTrueSmeared,
                             &b_W2DeltaPtTrueSmeared);
    fChain->SetBranchAddress("W2DeltaRTrueSmeared", &W2DeltaRTrueSmeared,
                             &b_W2DeltaRTrueSmeared);
    fChain->SetBranchAddress("W2DeltaMTrueSmeared", &W2DeltaMTrueSmeared,
                             &b_W2DeltaMTrueSmeared);
    fChain->SetBranchAddress("bJet2DeltaPtTrueSmeared",
                             &bJet2DeltaPtTrueSmeared,
                             &b_bJet2DeltaPtTrueSmeared);
    fChain->SetBranchAddress("bJet2DeltaRTrueSmeared", &bJet2DeltaRTrueSmeared,
                             &b_bJet2DeltaRTrueSmeared);
    fChain->SetBranchAddress("W2Daughter1DeltaPtTrueSmeared",
                             &W2Daughter1DeltaPtTrueSmeared,
                             &b_W2Daughter1DeltaPtTrueSmeared);
    fChain->SetBranchAddress("W2Daughter1DeltaRTrueSmeared",
                             &W2Daughter1DeltaRTrueSmeared,
                             &b_W2Daughter1DeltaRTrueSmeared);
    fChain->SetBranchAddress("W2Daughter2DeltaPtTrueSmeared",
                             &W2Daughter2DeltaPtTrueSmeared,
                             &b_W2Daughter2DeltaPtTrueSmeared);
    fChain->SetBranchAddress("W2Daughter2DeltaRTrueSmeared",
                             &W2Daughter2DeltaRTrueSmeared,
                             &b_W2Daughter2DeltaRTrueSmeared);
    fChain->SetBranchAddress("lightJet2DeltaPxTrueSmeared",
                             &lightJet2DeltaPxTrueSmeared,
                             &b_lightJet2DeltaPxTrueSmeared);
    fChain->SetBranchAddress("lightJet2DeltaPyTrueSmeared",
                             &lightJet2DeltaPyTrueSmeared,
                             &b_lightJet2DeltaPyTrueSmeared);
    fChain->SetBranchAddress("lightJet2DeltaRTrueSmeared",
                             &lightJet2DeltaRTrueSmeared,
                             &b_lightJet2DeltaRTrueSmeared);
    fChain->SetBranchAddress("top1DeltaPtTrueBest", &top1DeltaPtTrueBest,
                             &b_top1DeltaPtTrueBest);
    fChain->SetBranchAddress("top1DeltaRTrueBest", &top1DeltaRTrueBest,
                             &b_top1DeltaRTrueBest);
    fChain->SetBranchAddress("top1DeltaMTrueBest", &top1DeltaMTrueBest,
                             &b_top1DeltaMTrueBest);
    fChain->SetBranchAddress("W1DeltaPtTrueBest", &W1DeltaPtTrueBest,
                             &b_W1DeltaPtTrueBest);
    fChain->SetBranchAddress("W1DeltaRTrueBest", &W1DeltaRTrueBest,
                             &b_W1DeltaRTrueBest);
    fChain->SetBranchAddress("W1DeltaMTrueBest", &W1DeltaMTrueBest,
                             &b_W1DeltaMTrueBest);
    fChain->SetBranchAddress("bJet1DeltaPtTrueBest", &bJet1DeltaPtTrueBest,
                             &b_bJet1DeltaPtTrueBest);
    fChain->SetBranchAddress("bJet1DeltaRTrueBest", &bJet1DeltaRTrueBest,
                             &b_bJet1DeltaRTrueBest);
    fChain->SetBranchAddress("W1Daughter1DeltaPtTrueBest",
                             &W1Daughter1DeltaPtTrueBest,
                             &b_W1Daughter1DeltaPtTrueBest);
    fChain->SetBranchAddress("W1Daughter1DeltaRTrueBest",
                             &W1Daughter1DeltaRTrueBest,
                             &b_W1Daughter1DeltaRTrueBest);
    fChain->SetBranchAddress("W1Daughter2DeltaPtTrueBest",
                             &W1Daughter2DeltaPtTrueBest,
                             &b_W1Daughter2DeltaPtTrueBest);
    fChain->SetBranchAddress("W1Daughter2DeltaRTrueBest",
                             &W1Daughter2DeltaRTrueBest,
                             &b_W1Daughter2DeltaRTrueBest);
    fChain->SetBranchAddress("lightJet1DeltaPxTrueBest",
                             &lightJet1DeltaPxTrueBest,
                             &b_lightJet1DeltaPxTrueBest);
    fChain->SetBranchAddress("lightJet1DeltaPyTrueBest",
                             &lightJet1DeltaPyTrueBest,
                             &b_lightJet1DeltaPyTrueBest);
    fChain->SetBranchAddress("lightJet1DeltaRTrueBest",
                             &lightJet1DeltaRTrueBest,
                             &b_lightJet1DeltaRTrueBest);
    fChain->SetBranchAddress("top2DeltaPtTrueBest", &top2DeltaPtTrueBest,
                             &b_top2DeltaPtTrueBest);
    fChain->SetBranchAddress("top2DeltaRTrueBest", &top2DeltaRTrueBest,
                             &b_top2DeltaRTrueBest);
    fChain->SetBranchAddress("top2DeltaMTrueBest", &top2DeltaMTrueBest,
                             &b_top2DeltaMTrueBest);
    fChain->SetBranchAddress("W2DeltaPtTrueBest", &W2DeltaPtTrueBest,
                             &b_W2DeltaPtTrueBest);
    fChain->SetBranchAddress("W2DeltaRTrueBest", &W2DeltaRTrueBest,
                             &b_W2DeltaRTrueBest);
    fChain->SetBranchAddress("W2DeltaMTrueBest", &W2DeltaMTrueBest,
                             &b_W2DeltaMTrueBest);
    fChain->SetBranchAddress("bJet2DeltaPtTrueBest", &bJet2DeltaPtTrueBest,
                             &b_bJet2DeltaPtTrueBest);
    fChain->SetBranchAddress("bJet2DeltaRTrueBest", &bJet2DeltaRTrueBest,
                             &b_bJet2DeltaRTrueBest);
    fChain->SetBranchAddress("W2Daughter1DeltaPtTrueBest",
                             &W2Daughter1DeltaPtTrueBest,
                             &b_W2Daughter1DeltaPtTrueBest);
    fChain->SetBranchAddress("W2Daughter1DeltaRTrueBest",
                             &W2Daughter1DeltaRTrueBest,
                             &b_W2Daughter1DeltaRTrueBest);
    fChain->SetBranchAddress("W2Daughter2DeltaPtTrueBest",
                             &W2Daughter2DeltaPtTrueBest,
                             &b_W2Daughter2DeltaPtTrueBest);
    fChain->SetBranchAddress("W2Daughter2DeltaRTrueBest",
                             &W2Daughter2DeltaRTrueBest,
                             &b_W2Daughter2DeltaRTrueBest);
    fChain->SetBranchAddress("lightJet2DeltaPxTrueBest",
                             &lightJet2DeltaPxTrueBest,
                             &b_lightJet2DeltaPxTrueBest);
    fChain->SetBranchAddress("lightJet2DeltaPyTrueBest",
                             &lightJet2DeltaPyTrueBest,
                             &b_lightJet2DeltaPyTrueBest);
    fChain->SetBranchAddress("lightJet2DeltaRTrueBest",
                             &lightJet2DeltaRTrueBest,
                             &b_lightJet2DeltaRTrueBest);
    fChain->SetBranchAddress("top1DeltaPtSmearedBest", &top1DeltaPtSmearedBest,
                             &b_top1DeltaPtSmearedBest);
    fChain->SetBranchAddress("top1DeltaRSmearedBest", &top1DeltaRSmearedBest,
                             &b_top1DeltaRSmearedBest);
    fChain->SetBranchAddress("top1DeltaMSmearedBest", &top1DeltaMSmearedBest,
                             &b_top1DeltaMSmearedBest);
    fChain->SetBranchAddress("W1DeltaPtSmearedBest", &W1DeltaPtSmearedBest,
                             &b_W1DeltaPtSmearedBest);
    fChain->SetBranchAddress("W1DeltaRSmearedBest", &W1DeltaRSmearedBest,
                             &b_W1DeltaRSmearedBest);
    fChain->SetBranchAddress("W1DeltaMSmearedBest", &W1DeltaMSmearedBest,
                             &b_W1DeltaMSmearedBest);
    fChain->SetBranchAddress("bJet1DeltaPtSmearedBest",
                             &bJet1DeltaPtSmearedBest,
                             &b_bJet1DeltaPtSmearedBest);
    fChain->SetBranchAddress("bJet1DeltaRSmearedBest", &bJet1DeltaRSmearedBest,
                             &b_bJet1DeltaRSmearedBest);
    fChain->SetBranchAddress("bJet1DeltaMSmearedBest", &bJet1DeltaMSmearedBest,
                             &b_bJet1DeltaMSmearedBest);
    fChain->SetBranchAddress("W1Daughter1DeltaPtSmearedBest",
                             &W1Daughter1DeltaPtSmearedBest,
                             &b_W1Daughter1DeltaPtSmearedBest);
    fChain->SetBranchAddress("W1Daughter1DeltaRSmearedBest",
                             &W1Daughter1DeltaRSmearedBest,
                             &b_W1Daughter1DeltaRSmearedBest);
    fChain->SetBranchAddress("W1Daughter1DeltaMSmearedBest",
                             &W1Daughter1DeltaMSmearedBest,
                             &b_W1Daughter1DeltaMSmearedBest);
    fChain->SetBranchAddress("W1Daughter2DeltaPtSmearedBest",
                             &W1Daughter2DeltaPtSmearedBest,
                             &b_W1Daughter2DeltaPtSmearedBest);
    fChain->SetBranchAddress("W1Daughter2DeltaRSmearedBest",
                             &W1Daughter2DeltaRSmearedBest,
                             &b_W1Daughter2DeltaRSmearedBest);
    fChain->SetBranchAddress("W1Daughter2DeltaMSmearedBest",
                             &W1Daughter2DeltaMSmearedBest,
                             &b_W1Daughter2DeltaMSmearedBest);
    fChain->SetBranchAddress("lightJet1DeltaPxSmearedBest",
                             &lightJet1DeltaPxSmearedBest,
                             &b_lightJet1DeltaPxSmearedBest);
    fChain->SetBranchAddress("lightJet1DeltaPySmearedBest",
                             &lightJet1DeltaPySmearedBest,
                             &b_lightJet1DeltaPySmearedBest);
    fChain->SetBranchAddress("lightJet1DeltaRSmearedBest",
                             &lightJet1DeltaRSmearedBest,
                             &b_lightJet1DeltaRSmearedBest);
    fChain->SetBranchAddress("top2DeltaPtSmearedBest", &top2DeltaPtSmearedBest,
                             &b_top2DeltaPtSmearedBest);
    fChain->SetBranchAddress("top2DeltaRSmearedBest", &top2DeltaRSmearedBest,
                             &b_top2DeltaRSmearedBest);
    fChain->SetBranchAddress("top2DeltaMSmearedBest", &top2DeltaMSmearedBest,
                             &b_top2DeltaMSmearedBest);
    fChain->SetBranchAddress("W2DeltaPtSmearedBest", &W2DeltaPtSmearedBest,
                             &b_W2DeltaPtSmearedBest);
    fChain->SetBranchAddress("W2DeltaRSmearedBest", &W2DeltaRSmearedBest,
                             &b_W2DeltaRSmearedBest);
    fChain->SetBranchAddress("W2DeltaMSmearedBest", &W2DeltaMSmearedBest,
                             &b_W2DeltaMSmearedBest);
    fChain->SetBranchAddress("bJet2DeltaPtSmearedBest",
                             &bJet2DeltaPtSmearedBest,
                             &b_bJet2DeltaPtSmearedBest);
    fChain->SetBranchAddress("bJet2DeltaRSmearedBest", &bJet2DeltaRSmearedBest,
                             &b_bJet2DeltaRSmearedBest);
    fChain->SetBranchAddress("bJet2DeltaMSmearedBest", &bJet2DeltaMSmearedBest,
                             &b_bJet2DeltaMSmearedBest);
    fChain->SetBranchAddress("W2Daughter1DeltaPtSmearedBest",
                             &W2Daughter1DeltaPtSmearedBest,
                             &b_W2Daughter1DeltaPtSmearedBest);
    fChain->SetBranchAddress("W2Daughter1DeltaRSmearedBest",
                             &W2Daughter1DeltaRSmearedBest,
                             &b_W2Daughter1DeltaRSmearedBest);
    fChain->SetBranchAddress("W2Daughter1DeltaMSmearedBest",
                             &W2Daughter1DeltaMSmearedBest,
                             &b_W2Daughter1DeltaMSmearedBest);
    fChain->SetBranchAddress("W2Daughter2DeltaPtSmearedBest",
                             &W2Daughter2DeltaPtSmearedBest,
                             &b_W2Daughter2DeltaPtSmearedBest);
    fChain->SetBranchAddress("W2Daughter2DeltaRSmearedBest",
                             &W2Daughter2DeltaRSmearedBest,
                             &b_W2Daughter2DeltaRSmearedBest);
    fChain->SetBranchAddress("W2Daughter2DeltaMSmearedBest",
                             &W2Daughter2DeltaMSmearedBest,
                             &b_W2Daughter2DeltaMSmearedBest);
    fChain->SetBranchAddress("lightJet2DeltaPxSmearedBest",
                             &lightJet2DeltaPxSmearedBest,
                             &b_lightJet2DeltaPxSmearedBest);
    fChain->SetBranchAddress("lightJet2DeltaPySmearedBest",
                             &lightJet2DeltaPySmearedBest,
                             &b_lightJet2DeltaPySmearedBest);
    fChain->SetBranchAddress("lightJet2DeltaRSmearedBest",
                             &lightJet2DeltaRSmearedBest,
                             &b_lightJet2DeltaRSmearedBest);
    Notify();
}

Bool_t Plotter::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

void Plotter::Show(Long64_t entry)
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain)
        return;
    fChain->Show(entry);
}
Int_t Plotter::Cut(Long64_t entry)
{
    // This function may be called from Loop.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    return 1;
}
#endif // #ifdef Plotter_cxx
