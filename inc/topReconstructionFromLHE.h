//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr  8 15:44:02 2014 by ROOT version 5.34/04
// from TTree Physics/Physics
// found on file: lhe.root
//////////////////////////////////////////////////////////

// #ifndef topReconstructionFromLHE_h
// #define topReconstructionFromLHE_h

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
#include "topEventMinimizer.h"

using namespace std;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>
    XYZTLorentzVector;

// Fixed size dimensions of array or collections stored in the TTree if any.

/*struct handleEvent {
    TLorentzVector trueParticles [] = {
}*/

// Declaration and definition of structure
class handleEvent {
  public:
    // ROOT needs double pointer (**) to update containers, i.e., vectors
    map<string, XYZTLorentzVector *> trueParticles, smearedParticles,
        bestParticles, trueParticlesLH, smearedParticlesLH, bestParticlesLH;
    map<string, double> chiSquareds;
    bool leptonFlag;

    handleEvent(vector<string> &particleNames, vector<string> &names,
                vector<string> &chinames)
    {
        for (auto t = particleNames.begin(); t != particleNames.end(); ++t) {
            const string name = *t;
            trueParticles[name] = new XYZTLorentzVector;
            smearedParticles[name] = new XYZTLorentzVector;
            bestParticles[name] = new XYZTLorentzVector;
        }

        for (auto t = names.begin(); t < names.end(); ++t) {
            const string name = *t;
            trueParticlesLH[name] = new XYZTLorentzVector;
            smearedParticlesLH[name] = new XYZTLorentzVector;
            bestParticlesLH[name] = new XYZTLorentzVector;
        }

        for (auto t = chinames.begin(); t < chinames.end(); ++t) {
            const string name = *t;
            chiSquareds[name] = double();
        }
    }
    
    // clears a map of pointers
    template<class LV_pointer_map>
    void Clear_map(LV_pointer_map &m)
    {
        for (auto it = m.begin(); it != m.end(); ++it)
            delete (*it).second;
    }
    
    ~handleEvent()
    {
        // clears maps of pointers
        Clear_map(trueParticles);
        Clear_map(smearedParticles);
        Clear_map(bestParticles);
        Clear_map(trueParticlesLH);
        Clear_map(smearedParticlesLH);
        Clear_map(bestParticlesLH);
    }
};

class topReconstructionFromLHE
{
  public:
    // verbosity
    int debug_verbosity;
    
    // parameters
    double mTop, mW;
    double sigmaMTop, sigmaMW, sigmaPtLep, sigmaPhiLep, sigmaEtaLep,
           sigmaPhiJet, sigmaEtaJet;
    
    TTree *fChain; //! pointer to the analyzed TTree or TChain
    Int_t fCurrent; //! current Tree number in a TChain

    TFile *inFilePlot;
    TTree *inTreePlot;

    // Declaration of list of variables
    vector<string> particleNames;
    vector<string> names;
    vector<string> chinames;

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

    topReconstructionFromLHE(TTree *tree = NULL);
    ~topReconstructionFromLHE();
    void Set_def_parameters();
    Int_t Cut(Long64_t entry);
    Int_t GetEntry(Long64_t entry);
    Long64_t LoadTree(Long64_t entry);
    void Init(TTree *tree);
    void Loop(TString, const int, const int, const int);
    Bool_t Notify();
    void Show(Long64_t entry = -1);
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

    // For printing-only function
    void Print();

    // Branches
    int eventNumber;

    int innerMinStatus;
    int outerMinStatus;

    double outerMinEdm;
    double rel_error;

//     bool smearingSwitchedLightJetOrdering;

  private:
    typedef void (*FnPtr)(int, double, double, double, double);
    map<string, FnPtr> getBestMap;

    void DeclareMaps();
    void DeclareHists();
    void FillHists(handleEvent &);
    void FillHists_(const XYZTLorentzVector *, const XYZTLorentzVector *,
                    const string &, const string &,
                    hmap3 &);
    void FillLH(handleEvent &);
    void DeclareCanvases();
    void PlotHists();
    void DeclareOutBranches(handleEvent &);
    void DeclareInBranchesForPlotting(handleEvent &);
    void PrintMass(string, handleEvent&);
    void PrintPt(string, handleEvent&);
    void PrintPhi(string, handleEvent&);
    void PrintEta(string, handleEvent&);

    void Loop_diagnostics(handleEvent &);
    void Loop_load_eventh_SM(handleEvent &);
    void Loop_load_eventh_enu_SM(handleEvent &);
    void Loop_load_event_tt_SM(handleEvent &, topEventMinimizer &);
    void Loop_fill_results_SM(topEventMinimizer &, handleEvent &);
    void Print_smear_bs_SM(handleEvent &);
    void Print_smear_tt_SM(handleEvent &);
    
    string Get_pname(const string &, const string &);
    
    double Calc_rel_error(handleEvent &);
    double Calc_rel_error_(handleEvent &, const string);

    XYZTLorentzVector testvec;

    // handleEvent evh_outside;

    vector<string> varTypes = {"Pt", "Eta", "Phi", "M", "Px", "Py", "Pz", "E",
                               "Pt_", "Eta_", "Phi_", "M_", "Px_", "Py_", "Pz_", "E_"};
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

// #endif
