//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr  8 15:44:02 2014 by ROOT version 5.34/04
// from TTree Physics/Physics
// found on file: lhe.root
//////////////////////////////////////////////////////////

#ifndef topReconstructionFromLHE_h
#define topReconstructionFromLHE_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include "Math/GenVector/LorentzVector.h"
#include "TLorentzVector.h"

using namespace std;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>
    XYZTLorentzVector;

// Fixed size dimensions of array or collections stored in the TTree if any.

/*struct handleEvent {
    TLorentzVector trueParticles [] = {
}*/

class topReconstructionFromLHE
{
  public:
    TTree *fChain; //! pointer to the analyzed TTree or TChain
    Int_t fCurrent; //! current Tree number in a TChain

    TFile *inFilePlot;
    TTree *inTreePlot;

    // Declaration of list of variables
    static vector<string> particleNames;
    static vector<string> names;
    static vector<string> chinames;

    // Declaration and definition of structure
    struct handleEvent {

        map<string, XYZTLorentzVector *> trueParticles, smearedParticles,
            bestParticles, trueParticlesLH, smearedParticlesLH, bestParticlesLH;
        map<string, double *> chiSquareds;
        bool leptonFlag;

        handleEvent()
        {

            for (vector<string>::const_iterator t = particleNames.begin();
                 t < particleNames.end(); t++) {
                string name = *t;
                trueParticles[name] = new XYZTLorentzVector();
                smearedParticles[name] = new XYZTLorentzVector();
                bestParticles[name] = new XYZTLorentzVector();
            }

            for (vector<string>::const_iterator t = names.begin();
                 t < names.end(); t++) {
                string name = *t;
                trueParticlesLH[name] = new XYZTLorentzVector();
                smearedParticlesLH[name] = new XYZTLorentzVector();
                bestParticlesLH[name] = new XYZTLorentzVector();
            }

            for (vector<string>::const_iterator t = chinames.begin();
                 t < chinames.end(); t++) {
                string name = *t;
                chiSquareds[name] = new double;
            }
        }
    };

    // Declare maps to store histograms and canvases for plots of resolution and
    // plots of chi2
    typedef map<string, map<string, map<string, TH1D *>>> hmap3;
    typedef map<string, map<string, TH1D *>> hmap2;
    typedef map<string, TH1D *> hmap1;
    hmap3 histdif;

    typedef map<string, TCanvas *> cmap1;
    typedef map<string, cmap1> cmap2;
    cmap2 canvasdif;

    hmap1 histchi;
    cmap1 canvaschi;

    // Declaration of leaf types
    Int_t n_particles;
    vector<int> *PID;
    vector<double> *P_X;
    vector<double> *P_Y;
    vector<double> *P_Z;
    vector<double> *E;
    vector<double> *M;
    vector<int> *status;
    //   vector<int>     *particleID;
    //   vector<int>     *parent1ID;
    //   vector<int>     *parent2ID;

    // List of branches
    TBranch *b_n_particles; //!
    TBranch *b_PID;         //!
    TBranch *b_P_X;         //!
    TBranch *b_P_Y;         //!
    TBranch *b_P_Z;         //!
    TBranch *b_E;           //!
    TBranch *b_M;           //!
    TBranch *b_status;      //!
                            //   TBranch        *b_particleID;   //!
                            //   TBranch        *b_parent1ID;   //!
                            //   TBranch        *b_parent2ID;   //!

    topReconstructionFromLHE(TTree *tree = 0);
    virtual ~topReconstructionFromLHE();
    virtual Int_t Cut(Long64_t entry);
    virtual Int_t GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void Init(TTree *tree);
    virtual void Loop(TString, int, int);
    virtual Bool_t Notify();
    virtual void Show(Long64_t entry = -1);
    void printVector(XYZTLorentzVector &);
    void moveStatsBox(TH1D *);

    double deltaR(XYZTLorentzVector &, XYZTLorentzVector &);
    double deltaPhi(XYZTLorentzVector &, XYZTLorentzVector &);

    // store minimizer results in a tree;
    TFile *outFile;
    TTree *outTree;

    // Initialising functions
    void initOutput(TString, int);
    void initPlotting(TString);

    // Variables for plotting-only function
    void Plot(TString);
    TFile *outFilePlot;

    // Branches

    int eventNumber;

    int innerMinStatus;
    int outerMinStatus;

    double totalChi2;
    double topSystemChi2;
    double topMassChi2;
    double hadronicChi2;
    double nonTopChi2;

    double outerMinEdm;

    bool smearingSwitchedLightJetOrdering;

  private:
    typedef void (*FnPtr)(int, double, double, double, double);
    map<string, FnPtr> getBestMap;

    void DeclareMaps();
    void DeclareHists();
    void FillHists(handleEvent);
    void FillLH(handleEvent);
    void DeclareCanvases();
    void PlotHists();
    void DeclareOutBranches(handleEvent);
    void DeclareInBranchesForPlotting(handleEvent &);

    XYZTLorentzVector testvec;

    // handleEvent evh_outside;

    vector<string> varTypes = {"Pt", "Eta", "Phi", "Px", "Py", "M"};
    vector<string> difTypes = {"smearedTrue", "bestTrue"};

    vector<vector<string>> nameMap = {
        {"0", "Leptonic_Bottom", "bottom"},
        {"0", "Hadronic_Bottom", "antiBottom"},
        {"0", "Leptonic_Top", "top"},
        {"0", "Hadronic_Top", "antiTop"},
        {"0", "Leptonic_W", "Wplus"},
        {"0", "Hadronic_W", "Wminus"},
        {"0", "Lepton_or_AntiLepton", "antiLepton"},
        {"0", "Neutrino_or_AntiNeutrino", "neutrino"},
        {"0", "Quark_from_W", "qFromW"},
        {"0", "Antiquark_from_W", "qbarFromW"},
        {"0", "Higgs", "higgs"},
        {"0", "B_from_H", "bFromH"},
        {"0", "Bbar_from_H", "bbarFromH"},
        {"1", "Leptonic_Bottom", "antiBottom"},
        {"1", "Hadronic_Bottom", "bottom"},
        {"1", "Leptonic_Top", "antiTop"},
        {"1", "Hadronic_Top", "top"},
        {"1", "Leptonic_W", "Wminus"},
        {"1", "Hadronic_W", "Wplus"},
        {"1", "Lepton_or_AntiLepton", "lepton"},
        {"1", "Neutrino_or_AntiNeutrino", "antiNeutrino"},
        {"1", "Quark_from_W", "qFromW"},
        {"1", "Antiquark_from_W", "qbarFromW"},
        {"1", "Higgs", "higgs"},
        {"1", "B_from_H", "bFromH"},
        {"1", "Bbar_from_H", "bbarFromH"}};
};

#endif

#ifdef topReconstructionFromLHE_cxx
topReconstructionFromLHE::topReconstructionFromLHE(TTree *tree) : fChain(0)
{
    // if parameter tree is not specified (or zero), connect the file
    // used to generate this class and read the Tree.
    if (tree == 0) {
        //      TFile *f =
        //      (TFile*)gROOT->GetListOfFiles()->FindObject("skimmedntuple.root");
        TFile *f = (TFile *)gROOT->GetListOfFiles()->FindObject(
            "fixedsinglemuntuple.root");
        if (!f || !f->IsOpen()) {
            //         f = new TFile("skimmedntuple.root");
            f = new TFile("fixedsinglemuntuple.root");
        }
        f->GetObject("physics", tree);
    }
    Init(tree);

    string particleNameArray[15] = {
        "top",    "antiTop",   "bottom",       "antiBottom", "Wplus",
        "Wminus", "lepton",    "antiNeutrino", "antiLepton", "neutrino",
        "qFromW", "qbarFromW", "higgs",        "bFromH",     "bbarFromH"};
    particleNames.clear();
    for (int i = 0; i < 15; i++) {
        particleNames.push_back(particleNameArray[i]);
    }
    // particleNames(particleNameArray, particleNameArray+13);

    string namesArray[13] = {"Leptonic_Bottom",
                             "Hadronic_Bottom",
                             "Leptonic_Top",
                             "Hadronic_Top",
                             "Leptonic_W",
                             "Hadronic_W",
                             "Lepton_or_AntiLepton",
                             "Neutrino_or_AntiNeutrino",
                             "Quark_from_W",
                             "Antiquark_from_W",
                             "Higgs",
                             "B_from_H",
                             "Bbar_from_H"};
    names.clear();
    for (int i = 0; i < 13; i++) {
        names.push_back(namesArray[i]);
    }

    //    string chinamesArray [12] = {"total", "topSystem", "topMass",
    //    "hadronic", "nonTop", "leptonicBottom", "leptonicWMass",
    //    "hadronicWMass", "hadronicBottom", "lepton", "qFromW", "qbarFromW" };
    string chinamesArray[5] = {"total", "topSystem", "topMass", "hadronic",
                               "nonTop"};
    chinames.clear();
    for (int i = 0; i < 5; i++) {
        chinames.push_back(chinamesArray[i]);
    }
}

topReconstructionFromLHE::~topReconstructionFromLHE()
{
    if (!fChain)
        return;
    delete fChain->GetCurrentFile();
}

double topReconstructionFromLHE::deltaR(XYZTLorentzVector &p2,
                                        XYZTLorentzVector &p1)
{
    double thisDeltaEta = p1.Eta() - p2.Eta();
    double thisDeltaPhi = p1.Phi() - p2.Phi();
    while (thisDeltaPhi > 3.14159265359) {
        thisDeltaPhi -= 2 * 3.14159265359;
    }
    while (thisDeltaPhi < -3.14159265359) {
        thisDeltaPhi += 2 * 3.14159265359;
    }
    return sqrt(thisDeltaPhi * thisDeltaPhi + thisDeltaEta * thisDeltaEta);
}

double topReconstructionFromLHE::deltaPhi(XYZTLorentzVector &p2,
                                          XYZTLorentzVector &p1)
{
    double thisDeltaPhi = p1.Phi() - p2.Phi();
    while (thisDeltaPhi > 3.14159265359) {
        thisDeltaPhi -= 2 * 3.14159265359;
    }
    while (thisDeltaPhi < -3.14159265359) {
        thisDeltaPhi += 2 * 3.14159265359;
    }
    return thisDeltaPhi;
}

void topReconstructionFromLHE::initPlotting(TString dir)
{
    TString fileName = dir;
    if (dir != "")
        fileName += "/";
    fileName += "output_";
    fileName += ".root";

    float leftbound = -150.;
    float rightbound = 150.;
    float canvsize1 = 700;
    float canvsize2 = 700;
    int numbins = 100;

    outFile = new TFile(fileName, "RECREATE", "tree");
    // outTree = new TTree("tree","tree");
}

void topReconstructionFromLHE::initOutput(TString dir, int whichLoop)
{
    TString fileName = dir;
    if (dir != "")
        fileName += "/";
    fileName += "output_";
    fileName += whichLoop;
    fileName += ".root";

    float leftbound = -150.;
    float rightbound = 150.;
    float canvsize1 = 700;
    float canvsize2 = 700;
    int numbins = 100;

    outFile = new TFile(fileName, "RECREATE", "tree");
    outTree = new TTree("tree", "tree");

    // outTree->Branch( "smearingSwitchedLightJetOrdering" ,
    // &smearingSwitchedLightJetOrdering ) ;
}

Int_t topReconstructionFromLHE::GetEntry(Long64_t entry)
{
    // Read contents of entry.
    if (!fChain)
        return 0;
    return fChain->GetEntry(entry);
}
Long64_t topReconstructionFromLHE::LoadTree(Long64_t entry)
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

void topReconstructionFromLHE::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set object pointer
    PID = 0;
    P_X = 0;
    P_Y = 0;
    P_Z = 0;
    E = 0;
    M = 0;
    status = 0;
    //   particleID = 0;
    //   parent1ID = 0;
    //   parent2ID = 0;
    // Set branch addresses and branch pointers
    if (!tree)
        return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("n_particles", &n_particles, &b_n_particles);
    fChain->SetBranchAddress("PID", &PID, &b_PID);
    fChain->SetBranchAddress("P_X", &P_X, &b_P_X);
    fChain->SetBranchAddress("P_Y", &P_Y, &b_P_Y);
    fChain->SetBranchAddress("P_Z", &P_Z, &b_P_Z);
    fChain->SetBranchAddress("E", &E, &b_E);
    fChain->SetBranchAddress("M", &M, &b_M);
    fChain->SetBranchAddress("status", &status, &b_status);
    //   fChain->SetBranchAddress("particleID", &particleID, &b_particleID);
    //   fChain->SetBranchAddress("parent1ID", &parent1ID, &b_parent1ID);
    //   fChain->SetBranchAddress("parent2ID", &parent2ID, &b_parent2ID);
    Notify();
}

Bool_t topReconstructionFromLHE::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

void topReconstructionFromLHE::Show(Long64_t entry)
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain)
        return;
    fChain->Show(entry);
}
Int_t topReconstructionFromLHE::Cut(Long64_t entry)
{
    // This function may be called from Loop.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    return 1;
}
#endif // #ifdef topReconstructionFromLHE_cxx
