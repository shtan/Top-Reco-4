#ifndef topReconstructionFromLHE_methods_cxx
#define topReconstructionFromLHE_methods_cxx

#include "topReconstructionFromLHE.h"

topReconstructionFromLHE::topReconstructionFromLHE(TTree *tree) : fChain(0)
{
    debug = false;
    
    // if parameter tree is not specified (or zero), connect the file
    // used to generate this class and read the Tree.
    if (tree == NULL) {
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

    const int part_name_size = 15;
    string particleNameArray[part_name_size] = {
        "top",    "antiTop",   "bottom",       "antiBottom", "Wplus",
        "Wminus", "lepton",    "antiNeutrino", "antiLepton", "neutrino",
        "qFromW", "qbarFromW", "higgs",        "bFromH",     "bbarFromH"};
    particleNames.reserve(part_name_size);
    for (int i = 0; i < part_name_size; ++i) {
        particleNames.push_back(particleNameArray[i]);
    }
    // particleNames(particleNameArray, particleNameArray+13);

    const int name_size = 13;
    string namesArray[name_size] = {"Leptonic_Bottom",
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
    names.reserve(name_size);
    for (int i = 0; i < name_size; ++i) {
        names.push_back(namesArray[i]);
    }

    //    string chinamesArray [12] = {"total", "topSystem", "topMass",
    //    "hadronic", "nonTop", "leptonicBottom", "leptonicWMass",
    //    "hadronicWMass", "hadronicBottom", "lepton", "qFromW", "qbarFromW" };
    const int chi_size = 5;
    string chinamesArray[chi_size] = {"total", "topSystem", "topMass",
                                      "hadronic", "nonTop"};
    chinames.reserve(chi_size);
    for (int i = 0; i < chi_size; ++i) {
        chinames.push_back(chinamesArray[i]);
    }
    
    Set_def_parameters();
}

void topReconstructionFromLHE::Set_def_parameters()
{
    mTop = 173.;
    mW = 80.4;
    sigmaMTop = 2.0;
    sigmaMW = 2.085;
    sigmaPtLep = 0;
    sigmaPhiLep = 0;
    sigmaEtaLep = 0;
    sigmaPhiJet = 0.01;
    sigmaEtaJet = 0.01;
//         sigmaPtJet = 0.1;
//         sigmaPhiJet = 0.;
//         sigmaEtaJet = 0.;
//         sigmaPtJet = 0.1;
//         sigmaPtLep = 0.1;
//         sigmaPhiLep = 0.01;
//         sigmaEtaLep = 0.01;
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

//     float leftbound = -150.;
//     float rightbound = 150.;
//     float canvsize1 = 700;
//     float canvsize2 = 700;
//     int numbins = 100;

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

//     float leftbound = -150.;
//     float rightbound = 150.;
//     float canvsize1 = 700;
//     float canvsize2 = 700;
//     int numbins = 100;

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
    PID = NULL;
    P_X = NULL;
    P_Y = NULL;
    P_Z = NULL;
    E = NULL;
    M = NULL;
    status = NULL;
    //   particleID = 0;
    //   parent1ID = 0;
    //   parent2ID = 0;
    // Set branch addresses and branch pointers
    if (!tree)
        return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    // specific tree structure implemented by Shao Min
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
#endif // #ifdef topReconstructionFromLHE_methods_cxx
