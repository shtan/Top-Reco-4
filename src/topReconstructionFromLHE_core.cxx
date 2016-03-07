#ifndef topReconstructionFromLHE_core_cxx
#define topReconstructionFromLHE_core_cxx

#include "topReconstructionFromLHE.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TString.h>
#include <TH1F.h>
#include <TCanvas.h>
//#include "WDaughterEllipseCalculator.h"
#include <TPaveStats.h>
#include "Math/GenVector/LorentzVector.h"

using namespace std;
using namespace ROOT::Math;

void topReconstructionFromLHE::printVector(XYZTLorentzVector &v)
{
    cout << "px = " << v.Px() << endl;
    cout << "py = " << v.Py() << endl;
    cout << "pz = " << v.Pz() << endl;
    cout << "E  = " << v.E() << endl;
}

void topReconstructionFromLHE::DeclareOutBranches(handleEvent &evh)
{
    for (auto t = particleNames.begin(); t != particleNames.end(); ++t) {
        outTree->Branch((*t + "_smeared").c_str(), evh.smearedParticles[*t]);
        outTree->Branch((*t + "_true").c_str(), evh.trueParticles[*t]);
        outTree->Branch((*t + "_best").c_str(), evh.bestParticles[*t]);
    }
    outTree->Branch("leptonFlag", &evh.leptonFlag);

    for (auto t = names.begin(); t != names.end(); ++t) {
        outTree->Branch(("LH_" + *t + "_smeared").c_str(),
                        evh.smearedParticlesLH[*t]);
        outTree->Branch(("LH_" + *t + "_true").c_str(),
                        evh.trueParticlesLH[*t]);
        outTree->Branch(("LH_" + *t + "_best").c_str(),
                        evh.bestParticlesLH[*t]);
    }
    outTree->Branch("eventNumber", &eventNumber);
    outTree->Branch("innerMinStatus", &innerMinStatus);
    outTree->Branch("outerMinStatus", &outerMinStatus);
    outTree->Branch("outerMinEdm", &outerMinEdm);
    outTree->Branch("rel_error", &rel_error);

    for (auto t = chinames.begin(); t != chinames.end(); ++t) {
        outTree->Branch(("chiSquared_" + *t).c_str(), &(evh.chiSquareds[*t]));
    }

    // outTree->Branch( "evh", &evh );
}

void topReconstructionFromLHE::DeclareInBranchesForPlotting(handleEvent &evh)
{
    for (auto t = particleNames.begin(); t != particleNames.end(); ++t) {
        inTreePlot->SetBranchAddress((*t + "_smeared").c_str(),
                                     &(evh.smearedParticles[*t]));
        inTreePlot->SetBranchAddress((*t + "_true").c_str(),
                                     &(evh.trueParticles[*t]));
        inTreePlot->SetBranchAddress((*t + "_best").c_str(),
                                     &(evh.bestParticles[*t]));
    }
    inTreePlot->SetBranchAddress("leptonFlag", &evh.leptonFlag);

    for (auto t = names.begin(); t != names.end(); ++t) {
        inTreePlot->SetBranchAddress(("LH_" + *t + "_smeared").c_str(),
                                     &(evh.smearedParticlesLH[*t]));
        inTreePlot->SetBranchAddress(("LH_" + *t + "_true").c_str(),
                                     &(evh.trueParticlesLH[*t]));
        inTreePlot->SetBranchAddress(("LH_" + *t + "_best").c_str(),
                                     &(evh.bestParticlesLH[*t]));
    }
    inTreePlot->SetBranchAddress("eventNumber", &eventNumber);
    inTreePlot->SetBranchAddress("innerMinStatus", &innerMinStatus);
    inTreePlot->SetBranchAddress("outerMinStatus", &outerMinStatus);
    inTreePlot->SetBranchAddress("outerMinEdm", &outerMinEdm);
    inTreePlot->SetBranchAddress("rel_error", &rel_error);

    for (auto t = chinames.begin(); t != chinames.end(); ++t) {
        inTreePlot->SetBranchAddress(("chiSquared_" + *t).c_str(),
                                     &(evh.chiSquareds[*t]));
    }
}

void topReconstructionFromLHE::DeclareMaps()
{
    // typedef void (*FnPtr)(int, double, double, double, double);
    // getBestMap["top"] = topEventMinimizer::getTop;
}

void topReconstructionFromLHE::DeclareHists()
{
    for (auto name : names) {
        for (auto vtype : varTypes) {
            double lbound = -1;
            double rbound = 1;
            if (vtype == "Pt" || vtype == "Px" || vtype == "Py" ||
                vtype == "Pz") {
                lbound = -150;
                rbound = 150;
            }
            if (vtype == "Eta") {
                lbound = -5;
                rbound = 5;
            }
            if (vtype == "Phi") {
                lbound = -5;
                rbound = 5;
            }
            if (vtype == "M" || vtype == "E") {
                lbound = -150;
                rbound = 150;
            }

            for (auto diff : difTypes) {
                const string hname = name + "_" + vtype + "_" + diff;
                histdif[diff][vtype][name] =
                    new TH1D(hname.c_str(), hname.c_str(), 100, lbound, rbound);
            }
        }
    }

    for (auto name : chinames) {
        const string hname = "ChiSquared_" + name;
        histchi[name] = new TH1D(hname.c_str(), hname.c_str(), 100, 0, 10);
    }
}

void topReconstructionFromLHE::DeclareCanvases()
{
    for (auto name : names) {
        for (auto vtype : varTypes) {
            const string cname = "c_" + name + "_" + vtype;
            canvasdif[vtype][name] =
                new TCanvas(cname.c_str(), cname.c_str(), 700, 700);
        }
    }

    for (auto chiname : chinames) {
        const string cname = "c_ChiSquared_" + chiname;
        canvaschi[chiname] =
            new TCanvas(cname.c_str(), cname.c_str(), 700, 700);
    }
}

string topReconstructionFromLHE::Get_pname(const string &lep_flag,
                                           const string &name)
{
    string pname = "noname";
    for (auto it = nameMap.begin(); it != nameMap.end(); ++it) {
        if (lep_flag == it->at(0) and name == it->at(1))
            pname = it->at(2);
    }

    return pname;
}

void topReconstructionFromLHE::FillLH(handleEvent &evh)
{
    string leptonFlagStr = "0";
    if (evh.leptonFlag)
        leptonFlagStr = "1";

    for (auto name : names) {
        const string pname = Get_pname(leptonFlagStr, name);
        *evh.smearedParticlesLH[name] = *evh.smearedParticles[pname];
        *evh.bestParticlesLH[name] = *evh.bestParticles[pname];
        *evh.trueParticlesLH[name] = *evh.trueParticles[pname];
    }
}

void topReconstructionFromLHE::FillHists(handleEvent &evh)
{
    string leptonFlagStr = "0";
    if (evh.leptonFlag)
        leptonFlagStr = "1";

    for (auto name : names) {
        const string pname = Get_pname(leptonFlagStr, name);

        FillHists_(evh.trueParticles[pname], evh.smearedParticles[pname],
                   "smearedTrue", name, histdif);

        FillHists_(evh.trueParticles[pname], evh.bestParticles[pname],
                   "bestTrue", name, histdif);
    }

    for (auto h : histchi)
        h.second->Fill(evh.chiSquareds[h.first]);
}

// diff two loretz vectors for diff; first and third level of hmap3
void topReconstructionFromLHE::FillHists_(const XYZTLorentzVector *v1,
                                          const XYZTLorentzVector *v2,
                                          const string &DT, const string &name,
                                          hmap3 &hm)
{
    hm[DT]["Pt"][name]->Fill(v1->Pt() - v2->Pt());
    hm[DT]["Pt_"][name]->Fill((v1->Pt() - v2->Pt()) / v1->Pt());

    hm[DT]["Px"][name]->Fill(v1->Px() - v2->Px());
    hm[DT]["Px_"][name]->Fill((v1->Px() - v2->Px()) / v1->Px());
    hm[DT]["Py"][name]->Fill(v1->Py() - v2->Py());
    hm[DT]["Py_"][name]->Fill((v1->Py() - v2->Py()) / v1->Py());
    hm[DT]["Pz"][name]->Fill(v1->Pz() - v2->Pz());
    hm[DT]["Pz_"][name]->Fill((v1->Pz() - v2->Pz()) / v1->Pz());

    hm[DT]["M"][name]->Fill(v1->M() - v2->M());
    hm[DT]["M_"][name]->Fill((v1->M() - v2->M()) / v1->M());

    hm[DT]["E"][name]->Fill(v1->E() - v2->E());
    hm[DT]["E_"][name]->Fill((v1->E() - v2->E()) / v1->E());

    hm[DT]["Eta"][name]->Fill(v1->Eta() - v2->Eta());
    hm[DT]["Eta_"][name]->Fill((v1->Eta() - v2->Eta()) / v1->Eta());
    hm[DT]["Phi"][name]->Fill(v1->Phi() - v2->Phi());
    hm[DT]["Phi_"][name]->Fill((v1->Phi() - v2->Phi()) / v1->Phi());
}

void topReconstructionFromLHE::PlotHists()
{
    cout << "Drawing histograms..." << endl;
    int mc1 = 5;
    int mc2 = 4;
    //     int mc3 = 3;
    //     int mc0 = 0;

    int iter = 0;
    for (auto name : names) {
        for (auto vtype : varTypes) {
            string unit = "";
            if (vtype == "Pt" or vtype == "Px" or vtype == "Py" or
                vtype == "Pz" or vtype == "M" or vtype == "E") {
                unit = "(GeV)";
            }
            TH1D *h_smeared = histdif["smearedTrue"][vtype][name];
            TH1D *h_true = histdif["bestTrue"][vtype][name];
            TCanvas *canv = canvasdif[vtype][name];

            h_smeared->SetFillColor(mc1);
            h_smeared->SetLineColor(mc1);
            h_true->SetLineColor(mc2);

            if (vtype == "Pt" or vtype == "Px" or vtype == "Py" or
                vtype == "Phi" or vtype == "Pz" or vtype == "M" or
                vtype == "E" or vtype == "Eta") {
                h_smeared->GetXaxis()->SetTitle(
                    (vtype + " Residual " + unit).c_str());
            } else {
                h_smeared->GetXaxis()->SetTitle(
                    (vtype + " Resolution " + unit).c_str());
            }
            h_smeared->GetYaxis()->SetTitle("Events");
            h_smeared->SetTitle((name + "_" + vtype).c_str());
            h_smeared->SetMaximum(
                max(h_smeared->GetMaximum(), h_true->GetMaximum()) + 1);

            canv->cd();

            h_smeared->Draw("HIST");
            h_true->Draw("SAMES");

            moveStatsBox(h_smeared);

            canv->Write();
            // canv->Print((name + "_" + vtype + ".pdf").c_str());

            if (iter == 0) {
                canv->Print("./pdfplots2/plots.pdf[");
                canv->Print("./pdfplots2/plots.pdf");
            } else if (iter ==
                       (int)(names.size()) * (int)(varTypes.size()) - 1) {
                canv->Print("./pdfplots2/plots.pdf");
                canv->Print("./pdfplots2/plots.pdf]");
            } else {
                canv->Print("./pdfplots2/plots.pdf");
            }
            canv->ls();

            ++iter;
        }
    }

    int iterchi = 0;
    for (auto name : chinames) {
        histchi[name]->SetTitle(("chiSquared_" + name).c_str());
        canvaschi[name]->cd();
        histchi[name]->Draw("HIST");
        canvaschi[name]->Write();

        if (iterchi == 0) {
            canvaschi[name]->Print("./pdfplots2/chiplots.pdf[");
            canvaschi[name]->Print("./pdfplots2/chiplots.pdf");
        } else if (iterchi == (int)(chinames.size()) - 1) {
            canvaschi[name]->Print("./pdfplots2/chiplots.pdf");
            canvaschi[name]->Print("./pdfplots2/chiplots.pdf]");
        } else {
            canvaschi[name]->Print("./pdfplots2/chiplots.pdf");
        }
        ++iterchi;

        canvaschi[name]->ls();
    }
}

void topReconstructionFromLHE::PrintPt(string whichParticle, handleEvent &evh)
{
    cout << whichParticle << endl;
    cout << "Pt true = " << evh.trueParticlesLH[whichParticle]->Pt() << endl;
    cout << "Pt smeared = " << evh.smearedParticlesLH[whichParticle]->Pt()
         << endl;
    cout << "Pt best = " << evh.bestParticlesLH[whichParticle]->Pt() << endl;
    cout << "Pt smeared - true = "
         << evh.smearedParticlesLH[whichParticle]->Pt() -
                evh.trueParticlesLH[whichParticle]->Pt()
         << endl;
    cout << "Pt best - true = "
         << evh.bestParticlesLH[whichParticle]->Pt() -
                evh.trueParticlesLH[whichParticle]->Pt()
         << endl;
}

void topReconstructionFromLHE::PrintPhi(string whichParticle, handleEvent &evh)
{
    cout << whichParticle << endl;
    cout << "Phi true = " << evh.trueParticlesLH[whichParticle]->Phi() << endl;
    cout << "Phi smeared = " << evh.smearedParticlesLH[whichParticle]->Phi()
         << endl;
    cout << "Phi best = " << evh.bestParticlesLH[whichParticle]->Phi() << endl;
    cout << "Phi smeared - true = "
         << evh.smearedParticlesLH[whichParticle]->Phi() -
                evh.trueParticlesLH[whichParticle]->Phi()
         << endl;
    cout << "Phi best - true = "
         << evh.bestParticlesLH[whichParticle]->Phi() -
                evh.trueParticlesLH[whichParticle]->Phi()
         << endl;
}

void topReconstructionFromLHE::PrintEta(string whichParticle, handleEvent &evh)
{
    cout << whichParticle << endl;
    cout << "Eta true = " << evh.trueParticlesLH[whichParticle]->Eta() << endl;
    cout << "Eta smeared = " << evh.smearedParticlesLH[whichParticle]->Eta()
         << endl;
    cout << "Eta best = " << evh.bestParticlesLH[whichParticle]->Eta() << endl;
    cout << "Eta smeared - true = "
         << evh.smearedParticlesLH[whichParticle]->Eta() -
                evh.trueParticlesLH[whichParticle]->Eta()
         << endl;
    cout << "Eta best - true = "
         << evh.bestParticlesLH[whichParticle]->Eta() -
                evh.trueParticlesLH[whichParticle]->Eta()
         << endl;
}

void topReconstructionFromLHE::PrintMass(string whichParticle, handleEvent &evh)
{
    cout << whichParticle << endl;
    cout << "Mass true = " << evh.trueParticlesLH[whichParticle]->M() << endl;
    cout << "Mass smeared = " << evh.smearedParticlesLH[whichParticle]->M()
         << endl;
    cout << "Mass best = " << evh.bestParticlesLH[whichParticle]->M() << endl;
    cout << "Mass smeared - true = "
         << evh.smearedParticlesLH[whichParticle]->M() -
                evh.trueParticlesLH[whichParticle]->M()
         << endl;
    cout << "Mass best - true = "
         << evh.bestParticlesLH[whichParticle]->M() -
                evh.trueParticlesLH[whichParticle]->M()
         << endl;
}

void topReconstructionFromLHE::Plot(TString dir)
{
    inFilePlot = new TFile("output_files/output_0.root");
    inTreePlot = (TTree *)inFilePlot->Get("tree");

    outFilePlot = new TFile(dir + "/output_0_plots.root", "RECREATE");

    initPlotting(dir);

    handleEvent evh(particleNames, names, chinames);
    DeclareInBranchesForPlotting(evh);
    DeclareHists();

    int numEvents = inTreePlot->GetEntries();

    for (int i = 0; i < numEvents; ++i) {
        // cout<<"blah"<<endl;
        inTreePlot->GetEntry(i);
        // cout<<"balh"<<endl;
        if ((innerMinStatus == 0 or innerMinStatus == 1) and
            (outerMinStatus == 0 or outerMinStatus == 1)) {
            // if (*(evh.chiSquareds["total"]) < 1){
            // if ( (
            // (evh.bestParticlesLH["Leptonic_W"]->M()-evh.trueParticlesLH["Leptonic_W"]->M())
            // > -15 ) and (
            // (evh.bestParticlesLH["Leptonic_W"]->M()-evh.trueParticlesLH["Leptonic_W"]->M())
            // < 15 ) ){

            FillHists(evh);
            //}
            //}
        }
        // cout<<"after fillhists"<<endl;
    }
    cout << "Opening outFile" << endl;
    outFilePlot->cd();
    DeclareCanvases();
    PlotHists();
    outFilePlot->Write();
    outFilePlot->Close();

    cout << "Done!" << endl;
}

void topReconstructionFromLHE::Print()
{
    inFilePlot = new TFile("output_files/output_0.root");
    inTreePlot = (TTree *)inFilePlot->Get("tree");

    handleEvent evh(particleNames, names, chinames);
    DeclareInBranchesForPlotting(evh);
    DeclareHists();

    int numEvents = inTreePlot->GetEntries();

    for (int i = 0; i < numEvents; i++) {
        inTreePlot->GetEntry(i);
        if ((innerMinStatus == 0 or innerMinStatus == 1) and
            (outerMinStatus == 0 or outerMinStatus == 1)) {
            cout << endl;
            cout << "eventNumber = " << eventNumber << endl;
            PrintMass("Leptonic_W", evh);
            PrintPt("Leptonic_Bottom", evh);
            PrintPhi("Leptonic_Bottom", evh);
            PrintEta("Leptonic_Bottom", evh);
            PrintPt("Neutrino_or_AntiNeutrino", evh);
            PrintPhi("Neutrino_or_AntiNeutrino", evh);
            PrintEta("Neutrino_or_AntiNeutrino", evh);
            PrintMass("Leptonic_Top", evh);
            PrintMass("Hadronic_W", evh);
            PrintPt("Hadronic_Bottom", evh);
            PrintPt("Quark_from_W", evh);
            PrintPt("Antiquark_from_W", evh);

            for (vector<string>::const_iterator t = chinames.begin();
                 t < chinames.end(); t++) {
                cout << *t << " " << evh.chiSquareds[*t] << endl;
            }

            /*     cout<<"Px total = "<<
               evh.bestParticlesLH["Leptonic_Bottom"]->Px() <<endl;
                 cout<<evh.bestParticlesLH["Hadronic_Bottom"]->Px() <<endl;
                 cout<<evh.bestParticlesLH["Lepton_or_Antilepton"]->Px() <<endl;
                 cout<<evh.bestParticlesLH["Neutrino_or_AntiNeutrino"]->Px()
               <<endl;
                 cout<<evh.bestParticlesLH["Quark_from_W"]->Px() <<endl;
                 cout<<evh.bestParticlesLH["Antiquark_from_W"]->Px() <<endl;
                 cout<<evh.bestParticlesLH["B_from_H"]->Px() <<endl;
                 cout<<evh.bestParticlesLH["Bbar_from_H"]->Px() << endl; */

            cout << "Px total = "
                 << evh.bestParticlesLH["Leptonic_Bottom"]->Px() +
                        evh.bestParticlesLH["Hadronic_Bottom"]->Px() +
                        evh.bestParticlesLH["Lepton_or_AntiLepton"]->Px() +
                        evh.bestParticlesLH["Neutrino_or_AntiNeutrino"]->Px() +
                        evh.bestParticlesLH["Quark_from_W"]->Px() +
                        evh.bestParticlesLH["Antiquark_from_W"]->Px() +
                        evh.bestParticlesLH["B_from_H"]->Px() +
                        evh.bestParticlesLH["Bbar_from_H"]->Px()
                 << endl;
            cout << "Py total = "
                 << evh.bestParticlesLH["Leptonic_Bottom"]->Py() +
                        evh.bestParticlesLH["Hadronic_Bottom"]->Py() +
                        evh.bestParticlesLH["Lepton_or_AntiLepton"]->Py() +
                        evh.bestParticlesLH["Neutrino_or_AntiNeutrino"]->Py() +
                        evh.bestParticlesLH["Quark_from_W"]->Py() +
                        evh.bestParticlesLH["Antiquark_from_W"]->Py() +
                        evh.bestParticlesLH["B_from_H"]->Py() +
                        evh.bestParticlesLH["Bbar_from_H"]->Py()
                 << endl;
            cout << "Top Px total = "
                 << evh.bestParticlesLH["Leptonic_Bottom"]->Px() +
                        evh.bestParticlesLH["Hadronic_Bottom"]->Px() +
                        evh.bestParticlesLH["Lepton_or_AntiLepton"]->Px() +
                        evh.bestParticlesLH["Neutrino_or_AntiNeutrino"]->Px() +
                        evh.bestParticlesLH["Quark_from_W"]->Px() +
                        evh.bestParticlesLH["Antiquark_from_W"]->Px()
                 << endl;
            cout << "Nontop Px total = "
                 << evh.bestParticlesLH["B_from_H"]->Px() +
                        evh.bestParticlesLH["Bbar_from_H"]->Px()
                 << endl;

            cout << "Top Py total = "
                 << evh.bestParticlesLH["Leptonic_Bottom"]->Py() +
                        evh.bestParticlesLH["Hadronic_Bottom"]->Py() +
                        evh.bestParticlesLH["Lepton_or_AntiLepton"]->Py() +
                        evh.bestParticlesLH["Neutrino_or_AntiNeutrino"]->Py() +
                        evh.bestParticlesLH["Quark_from_W"]->Py() +
                        evh.bestParticlesLH["Antiquark_from_W"]->Py()
                 << endl;
            cout << "Nontop Py total = "
                 << evh.bestParticlesLH["B_from_H"]->Py() +
                        evh.bestParticlesLH["Bbar_from_H"]->Py()
                 << endl;
        }
        // cout<<"after fillhists"<<endl;
    }
}

void topReconstructionFromLHE::Loop(TString dir, const int whichLoop,
                                    const int maxLoops, const int max_evts)
{
    // int whichLoop = 1; int maxLoops = 1;
    //   In a ROOT session, you can do:
    //      Root > .L topReconstructionFromLHE.C
    //      Root > topReconstructionFromLHE t
    //      Root > t.GetEntry(12); // Fill t data members with entry number 12
    //      Root > t.Show();       // Show values of entry 12
    //      Root > t.Show(16);     // Read and show values of entry 16
    //      Root > t.Loop();       // Loop on all entries
    //

    //     This is the loop skeleton where:
    //    jentry is the global entry number in the chain
    //    ientry is the entry number in the current Tree
    //  Note that the argument to GetEntry must be:
    //    jentry for TChain::GetEntry
    //    ientry for TTree::GetEntry and TBranch::GetEntry
    //
    //       To read only selected branches, Insert statements like:
    // METHOD1:
    //    fChain->SetBranchStatus("*",0);  // disable all branches
    //    fChain->SetBranchStatus("branchname",1);  // activate branchname
    // METHOD2: replace line
    //    fChain->GetEntry(jentry);       //read all branches
    // by  b_branchname->GetEntry(ientry); //read only this branch
    if (fChain == NULL)
        return;
    TRandom3 rand;
    rand.SetSeed(10);

    initOutput(dir, whichLoop);

    // gStyle->SetOptStat(101110);
    // TCanvas c1("c1","c1",500,500);
    // TH1D h_chi2("h_chi2","",100,0,100);

    const Long64_t nentries = fChain->GetEntriesFast();
    // Long64_t nentries = intree->GetEntries();
    // int jStart = whichLoop*(nentries/maxLoops) +
    // ((whichLoop>(maxLoops-nentries%maxLoops))?(whichLoop+nentries%maxLoops-maxLoops):0);
    const int jStart = 0;
    int jFinish = max_evts;
    if (max_evts < 0)
        jFinish = jStart + (nentries + whichLoop) / maxLoops;

    // Declare structure (used as a container) to store all the particles
    handleEvent evh(particleNames, names, chinames);

    DeclareHists();
    DeclareOutBranches(evh);

    if (debug_verbosity >= 2)
        cout << "Loop start, finish, n_entries: " << jStart << ", " << jFinish
             << ", " << nentries << endl;
    for (Long64_t jentry = jStart; jentry < jFinish; ++jentry) {
        if (debug_verbosity >= 2)
            cout << "BEGINNING BRANCH NUMBER " << jentry << endl;
        Long64_t ientry = LoadTree(jentry);
        // intree->GetEntry(jentry);
        if (ientry < 0)
            break;
        // if (jentry > 0) break;
        fChain->GetEntry(jentry); // Long64_t nbytes = 0, nb = 0;
        //         nbytes += nb;
        // if (Cut(ientry) < 0) continue;

        //      bool leptonFlag = 0; //1 if lepton and antineutrino (i.e. tbar
        //      branch is the leptonic one), 0 otherwise.
        evh.leptonFlag = false; // 1 if lepton and antineutrino (i.e. tbar
                                // branch is the leptonic one), 0 otherwise.

        // double METx(0.), METy(0.);
        //         int iJet = 0;

        if (debug_verbosity >= 2) {
            cout << "Before loading event.\n";
            Print_smear_bs_SM(evh);
        }

        // Load evh, follows Shao Min's tree definitions
        Loop_load_eventh_SM(evh);

        if (debug_verbosity >= 2) {
            cout << "After loading event.\n";
            Print_smear_bs_SM(evh);
            cout << "Setting light partons" << endl;
        }
        vector<XYZTLorentzVector> smearedOtherLightPartons;
        for (int ii = 0; ii < ((int)PID->size() - 21); ++ii) {
            XYZTLorentzVector temp;
            temp.SetPxPyPzE(P_X->at(ii + 21), P_Y->at(ii + 21),
                            P_Z->at(ii + 21), E->at(ii + 21));
            smearedOtherLightPartons.push_back(temp);
        }
        if (debug_verbosity >= 2) {
            cout << "smearedOtherLightPartons.size() = "
                 << smearedOtherLightPartons.size() << endl
                 << "Setting leptons and neutrinos" << endl;
        }
        Loop_load_eventh_enu_SM(evh);
        if (debug_verbosity >= 2)
            Print_smear_tt_SM(evh);

        // Set smeared higgs by adding bFromH and bbarFromH
        *(evh.smearedParticles["higgs"]) = *(evh.smearedParticles["bFromH"]) +
                                           *(evh.smearedParticles["bbarFromH"]);
        if (debug_verbosity >= 2)
            cout << "Setting non-top objects" << endl;
        vector<XYZTLorentzVector> nonTopObjects;
        vector<double> nonTopObjectPtWidths;
        vector<double> nonTopObjectEtaWidths;
        vector<double> nonTopObjectPhiWidths;

        nonTopObjects.push_back(*(evh.smearedParticles["bFromH"]));
        nonTopObjectPtWidths.push_back(
            sqrt(evh.smearedParticles["bFromH"]->Pt()));
        nonTopObjectEtaWidths.push_back(0.01);
        nonTopObjectPhiWidths.push_back(0.01);

        nonTopObjects.push_back(*(evh.smearedParticles["bbarFromH"]));
        nonTopObjectPtWidths.push_back(
            sqrt(evh.smearedParticles["bbarFromH"]->Pt()));
        nonTopObjectEtaWidths.push_back(0.01);
        nonTopObjectPhiWidths.push_back(0.01);
        // nonTopObjects.push_back(*(evh.smearedParticles["bFromH"]) +
        //     *(evh.smearedParticles["bbarFromH"]));
        // nonTopObjectPtWidths.push_back(30);
        // nonTopObjectEtaWidths.push_back(0.01);
        // nonTopObjectPhiWidths.push_back(0.01);

        for (auto particle : smearedOtherLightPartons) {
            nonTopObjects.push_back(particle);
            nonTopObjectPtWidths.push_back(sqrt(particle.pt()));
            nonTopObjectEtaWidths.push_back(0.01);
            nonTopObjectPhiWidths.push_back(0.01);
        }
        if (debug_verbosity >= 2) {
            cout << "nonTopObjects.size() = " << nonTopObjects.size() << endl
                 << "Setup top event minimizer." << endl;
        }

        // FIXME these quantities are not well known in real data
        const double totalTopPx = evh.smearedParticles["top"]->Px() +
                                  evh.smearedParticles["antiTop"]->Px();
        const double totalTopPy = evh.smearedParticles["top"]->Py() +
                                  evh.smearedParticles["antiTop"]->Py();
        const double totalTopPz = evh.smearedParticles["top"]->Pz() +
                                  evh.smearedParticles["antiTop"]->Pz();

        topEventMinimizer ev(nonTopObjects, nonTopObjectPtWidths,
                             nonTopObjectEtaWidths, nonTopObjectPhiWidths, mTop,
                             sigmaMTop, mW, sigmaMW, totalTopPx, totalTopPy,
                             totalTopPz);

        if (debug_verbosity >= 2)
            cout << "Before add tops:" << endl;
        Loop_load_event_tt_SM(evh, ev);

        // ev.printTopConstituents();
        ev.initializeDeltas();
        // ev.calcTopMassRanges();
        if (debug_verbosity >= 2)
            cout << "After add tops & initialise deltas" << endl;

        // ev.findStartingValues(50);
        // TString thisPlotName="ellipses_";
        // thisPlotName+=jentry;
        // ev.plotEllipses(thisPlotName);

        //         ev.minimizeNonTopChiSquare();

        ev.minimizeTotalChiSquare();
        if (debug_verbosity >= 2)
            cout << "After minimizeTotalChiSquare()." << endl;

        // h_chi2.Fill(ev.getChiSquare());

        // Fill branches
        eventNumber = jentry;

        // old stuff: can delete later
        innerMinStatus = ev.innerMinStatus;
        outerMinStatus = ev.outerMinStatus;
        outerMinEdm = ev.outerMin_Edm;

        // NEW Get best values
        if (debug_verbosity >= 2)
            cout << "Getting best values" << endl;
        Loop_fill_results_SM(ev, evh);

        // Set best higgs by adding bFromH and bbarFromH
        *evh.bestParticles["higgs"] =
            *(evh.bestParticles["bFromH"]) + *(evh.bestParticles["bbarFromH"]);

        // Fill Hists
        if (debug_verbosity >= 2)
            cout << "Filling hists" << endl;
        if ((innerMinStatus == 0 or innerMinStatus == 1) and
            (outerMinStatus == 0 or outerMinStatus == 1)) {
            FillHists(evh);
        }
        FillLH(evh);

        if (debug_verbosity >= 1)
            Loop_diagnostics(evh);

        outTree->Fill();

    } // end loop over events

    cout << "Opening outFile" << endl;
    outFile->cd();
    DeclareCanvases();
    PlotHists();
    outFile->Write();
    outFile->Close();

    cout << "Done!" << endl;
}

void topReconstructionFromLHE::Loop_load_eventh_SM(handleEvent &evh)
{
    // smearedLepton.SetPxPyPzE(P_X->at(0),P_Y->at(0),P_Z->at(0),E->at(0));
    evh.smearedParticles["bottom"]->SetPxPyPzE(P_X->at(1), P_Y->at(1),
                                               P_Z->at(1), E->at(1));
    evh.smearedParticles["antiBottom"]->SetPxPyPzE(P_X->at(2), P_Y->at(2),
                                                   P_Z->at(2), E->at(2));
    evh.smearedParticles["qFromW"]->SetPxPyPzE(P_X->at(3), P_Y->at(3),
                                               P_Z->at(3), E->at(3));
    evh.smearedParticles["qbarFromW"]->SetPxPyPzE(P_X->at(4), P_Y->at(4),
                                                  P_Z->at(4), E->at(4));
    evh.smearedParticles["bFromH"]->SetPxPyPzE(P_X->at(5), P_Y->at(5),
                                               P_Z->at(5), E->at(5));
    evh.smearedParticles["bbarFromH"]->SetPxPyPzE(P_X->at(6), P_Y->at(6),
                                                  P_Z->at(6), E->at(6));

    // lepton.SetPxPyPzE(P_X->at(7),P_Y->at(7),P_Z->at(7),E->at(7));
    // neutrino.SetPxPyPzE(P_X->at(18),P_Y->at(18),P_Z->at(18),E->at(18));
    evh.trueParticles["bottom"]->SetPxPyPzE(P_X->at(8), P_Y->at(8), P_Z->at(8),
                                            E->at(8));
    evh.trueParticles["antiBottom"]->SetPxPyPzE(P_X->at(9), P_Y->at(9),
                                                P_Z->at(9), E->at(9));
    evh.trueParticles["qFromW"]->SetPxPyPzE(P_X->at(10), P_Y->at(10),
                                            P_Z->at(10), E->at(10));
    evh.trueParticles["qbarFromW"]->SetPxPyPzE(P_X->at(11), P_Y->at(11),
                                               P_Z->at(11), E->at(11));
    evh.trueParticles["bFromH"]->SetPxPyPzE(P_X->at(12), P_Y->at(12),
                                            P_Z->at(12), E->at(12));
    evh.trueParticles["bbarFromH"]->SetPxPyPzE(P_X->at(13), P_Y->at(13),
                                               P_Z->at(13), E->at(13));

    evh.trueParticles["higgs"]->SetPxPyPzE(P_X->at(19), P_Y->at(19),
                                           P_Z->at(19), E->at(19));
    evh.trueParticles["top"]->SetPxPyPzE(P_X->at(14), P_Y->at(14), P_Z->at(14),
                                         E->at(14));
    evh.trueParticles["antiTop"]->SetPxPyPzE(P_X->at(15), P_Y->at(15),
                                             P_Z->at(15), E->at(15));
    evh.trueParticles["Wplus"]->SetPxPyPzE(P_X->at(16), P_Y->at(16),
                                           P_Z->at(16), E->at(16));
    evh.trueParticles["Wminus"]->SetPxPyPzE(P_X->at(17), P_Y->at(17),
                                            P_Z->at(17), E->at(17));
}

void topReconstructionFromLHE::Loop_load_eventh_enu_SM(handleEvent &evh)
{
    if (PID->at(0) == 13) {
        evh.smearedParticles["lepton"]->SetPxPyPzE(P_X->at(0), P_Y->at(0),
                                                   P_Z->at(0), E->at(0));
        // define smeared anti neutrino as the MET; Pz set to zero; E set
        // based on zero mass and Px and Py
        evh.smearedParticles["antiNeutrino"]->SetPxPyPzE(
            P_X->at(20), P_Y->at(20), 0,
            sqrt(pow(P_X->at(20), 2) + pow(P_Y->at(20), 2)));

        evh.trueParticles["lepton"]->SetPxPyPzE(P_X->at(7), P_Y->at(7),
                                                P_Z->at(7), E->at(7));
        evh.trueParticles["antiNeutrino"]->SetPxPyPzE(P_X->at(18), P_Y->at(18),
                                                      P_Z->at(18), E->at(18));

        *evh.smearedParticles["Wminus"] =
            *(evh.smearedParticles["lepton"]) +
            *(evh.smearedParticles["antiNeutrino"]);
        *evh.smearedParticles["Wplus"] = *(evh.smearedParticles["qFromW"]) +
                                         *(evh.smearedParticles["qbarFromW"]);
        *evh.smearedParticles["top"] = *(evh.smearedParticles["Wplus"]) +
                                       *(evh.smearedParticles["bottom"]);
        *evh.smearedParticles["antiTop"] =
            *(evh.smearedParticles["Wminus"]) +
            *(evh.smearedParticles["antiBottom"]);

        evh.leptonFlag = true;
    } else if (PID->at(0) == -13) {
        evh.smearedParticles["antiLepton"]->SetPxPyPzE(P_X->at(0), P_Y->at(0),
                                                       P_Z->at(0), E->at(0));
        // define smeared neutrino as the MET; Pz set to zero; E set based
        // on zero mass and Px and Py
        evh.smearedParticles["neutrino"]->SetPxPyPzE(
            P_X->at(20), P_Y->at(20), 0,
            sqrt(pow(P_X->at(20), 2) + pow(P_Y->at(20), 2)));

        evh.trueParticles["antiLepton"]->SetPxPyPzE(P_X->at(7), P_Y->at(7),
                                                    P_Z->at(7), E->at(7));
        evh.trueParticles["neutrino"]->SetPxPyPzE(P_X->at(18), P_Y->at(18),
                                                  P_Z->at(18), E->at(18));

        *evh.smearedParticles["Wplus"] = *(evh.smearedParticles["antiLepton"]) +
                                         *(evh.smearedParticles["neutrino"]);
        *evh.smearedParticles["Wminus"] = *(evh.smearedParticles["qFromW"]) +
                                          *(evh.smearedParticles["qbarFromW"]);
        *evh.smearedParticles["top"] = *(evh.smearedParticles["Wplus"]) +
                                       *(evh.smearedParticles["bottom"]);
        *evh.smearedParticles["antiTop"] =
            *(evh.smearedParticles["Wminus"]) +
            *(evh.smearedParticles["antiBottom"]);

        evh.leptonFlag = false;
    }
}

void topReconstructionFromLHE::Loop_load_event_tt_SM(handleEvent &evh,
                                                     topEventMinimizer &ev)
{
    if (evh.leptonFlag == false) {
        const auto b = evh.smearedParticles["bottom"];
        const auto al = evh.smearedParticles["antiLepton"];
        ev.addLeptonicTop(b->Px(), b->Py(), b->Pz(), b->E(), sqrt(b->Pt()),
                          sigmaEtaJet, sigmaPhiJet, al->Px(), al->Py(),
                          al->Pz(), al->E(), sigmaPtLep, sigmaEtaLep,
                          sigmaPhiLep, mTop, sigmaMTop, mW, sigmaMW);

        const auto bbar = evh.smearedParticles["antiBottom"];
        const auto Wq1 = evh.smearedParticles["qFromW"];
        const auto Wq2 = evh.smearedParticles["qbarFromW"];
        ev.addHadronicTop(bbar->Px(), bbar->Py(), bbar->Pz(), bbar->E(),
                          sqrt(bbar->Pt()), sigmaEtaJet, sigmaPhiJet, Wq1->Px(),
                          Wq1->Py(), Wq1->Pz(), Wq1->E(), sqrt(Wq1->Pt()),
                          sigmaEtaJet, sigmaPhiJet, Wq2->Px(), Wq2->Py(),
                          Wq2->Pz(), Wq2->E(), sqrt(Wq2->Pt()), sigmaEtaJet,
                          sigmaPhiJet, mTop, sigmaMTop, mW, sigmaMW);
    }

    if (evh.leptonFlag == true) {
        const auto b = evh.smearedParticles["bottom"];
        const auto Wq1 = evh.smearedParticles["qFromW"];
        const auto Wq2 = evh.smearedParticles["qbarFromW"];
        ev.addHadronicTop(b->Px(), b->Py(), b->Pz(), b->E(), sqrt(b->Pt()),
                          sigmaEtaJet, sigmaPhiJet, Wq1->Px(), Wq1->Py(),
                          Wq1->Pz(), Wq1->E(), sqrt(Wq1->Pt()), sigmaEtaJet,
                          sigmaPhiJet, Wq2->Px(), Wq2->Py(), Wq2->Pz(),
                          Wq2->E(), sqrt(Wq2->Pt()), sigmaEtaJet, sigmaPhiJet,
                          mTop, sigmaMTop, mW, sigmaMW);

        const auto bbar = evh.smearedParticles["antiBottom"];
        const auto l = evh.smearedParticles["lepton"];
        ev.addLeptonicTop(bbar->Px(), bbar->Py(), bbar->Pz(), bbar->E(),
                          sqrt(bbar->Pt()), sigmaEtaJet, sigmaPhiJet, l->Px(),
                          l->Py(), l->Pz(), l->E(), sigmaPtLep, sigmaEtaLep,
                          sigmaPhiLep, mTop, sigmaMTop, mW, sigmaMW);
    }
}

void topReconstructionFromLHE::Loop_fill_results_SM(topEventMinimizer &ev,
                                                    handleEvent &evh)
{
    // Minimizer stuff
    evh.chiSquareds["nonTop"] = ev.getNonTopChiSquare();
    evh.chiSquareds["hadronic"] = ev.getHadronicChiSquare();
    evh.chiSquareds["topMass"] = ev.getTopMassChiSquare();
    evh.chiSquareds["topSystem"] = ev.getTopChiSquare();
    evh.chiSquareds["total"] = ev.getChiSquare();
    evh.chiSquareds["qbarFromW"] = ev.getHadronicChiSquare();

    evh.chiSquareds["top_topMass"] = ev.getOneTopMassChiSquare(0);
    evh.chiSquareds["antiTop_topMass"] = ev.getOneTopMassChiSquare(1);
    evh.chiSquareds["bottom"] = ev.getOneBChiSquare(0);
    evh.chiSquareds["antiBottom"] = ev.getOneBChiSquare(1);
    evh.chiSquareds["Wplus_Wmass"] = ev.getOneWMassChiSquare(0);
    evh.chiSquareds["Wminus_Wmass"] = ev.getOneWMassChiSquare(1);

    if (evh.leptonFlag == 0) {
        evh.chiSquareds["leptonicTopMass"] = ev.getOneTopMassChiSquare(0);
        evh.chiSquareds["hadronicTopMass"] = ev.getOneTopMassChiSquare(1);
        evh.chiSquareds["leptonicBottom"] = ev.getOneBChiSquare(0);
        evh.chiSquareds["hadronicBottom"] = ev.getOneBChiSquare(1);
        evh.chiSquareds["leptonicWMass"] = ev.getOneWMassChiSquare(0);
        evh.chiSquareds["hadronicWMass"] = ev.getOneWMassChiSquare(1);
        evh.chiSquareds["qFromW"] = ev.getOneWDaughter1ChiSquare(1);
        evh.chiSquareds["lepton"] = ev.getOneWDaughter1ChiSquare(0);
        evh.chiSquareds["lepton_or_antilepton"] =
            ev.getOneWDaughter1ChiSquare(0);

    } else if (evh.leptonFlag == 1) {
        evh.chiSquareds["leptonicTopMass"] = ev.getOneTopMassChiSquare(1);
        evh.chiSquareds["hadronicTopMass"] = ev.getOneTopMassChiSquare(0);
        evh.chiSquareds["leptonicBottom"] = ev.getOneBChiSquare(1);
        evh.chiSquareds["hadronicBottom"] = ev.getOneBChiSquare(0);
        evh.chiSquareds["leptonicWMass"] = ev.getOneWMassChiSquare(1);
        evh.chiSquareds["hadronicWMass"] = ev.getOneWMassChiSquare(0);
        evh.chiSquareds["qFromW"] = ev.getOneWDaughter1ChiSquare(0);
        evh.chiSquareds["antiLepton"] = ev.getOneWDaughter1ChiSquare(1);
        evh.chiSquareds["lepton_or_antilepton"] =
            ev.getOneWDaughter1ChiSquare(1);
    }

    // for (vector<string>::const_iterator t = chinames.begin(); t <
    // chinames.end(); t++){
    //    cout<<*t<<" " << *(evh.chiSquareds[*t])<<endl;
    //}

    *evh.bestParticles["top"] = ev.getConverter("getTop", 0);
    *evh.bestParticles["bottom"] = ev.getConverter("getBJet", 0);
    *evh.bestParticles["Wplus"] = ev.getConverter("getW", 0);
    // FIXME
    if (evh.leptonFlag == 0) {
        *evh.bestParticles["antiLepton"] = ev.getConverter("getWDaughter1", 0);
        *evh.bestParticles["neutrino"] = ev.getConverter("getWDaughter2", 0);
    } else {
        *evh.bestParticles["qFromW"] = ev.getConverter("getWDaughter1", 0);
        *evh.bestParticles["qbarFromW"] = ev.getConverter("getWDaughter2", 0);
    }
    *evh.bestParticles["bFromH"] = ev.getConverter("getNonTopObject4", 0);

    *evh.bestParticles["antiTop"] = ev.getConverter("getTop", 1);
    *evh.bestParticles["antiBottom"] = ev.getConverter("getBJet", 1);
    *evh.bestParticles["Wminus"] = ev.getConverter("getW", 1);
    // FIXME
    if (evh.leptonFlag == 0) {
        *evh.bestParticles["qFromW"] = ev.getConverter("getWDaughter1", 1);
        *evh.bestParticles["qbarFromW"] = ev.getConverter("getWDaughter2", 1);
    } else {
        *evh.bestParticles["lepton"] = ev.getConverter("getWDaughter1", 1);
        *evh.bestParticles["antiNeutrino"] =
            ev.getConverter("getWDaughter2", 1);
    }
    *evh.bestParticles["bbarFromH"] = ev.getConverter("getNonTopObject4", 1);
}

void topReconstructionFromLHE::Loop_diagnostics(handleEvent &evh)
{
    XYZTLorentzVector gen_all =
        *evh.smearedParticles["bottom"] + *evh.smearedParticles["antiBottom"] +
        *evh.smearedParticles["qFromW"] + *evh.smearedParticles["qbarFromW"] +
        *evh.smearedParticles["bFromH"] + *evh.smearedParticles["bbarFromH"] +
        *evh.smearedParticles["lepton"] + *evh.smearedParticles["antiNeutrino"];
    XYZTLorentzVector best_all =
        *evh.bestParticles["top"] + *evh.bestParticles["antiTop"] +
        *evh.bestParticles["bFromH"] + *evh.bestParticles["bbarFromH"];
    XYZTLorentzVector best_tt =
        *evh.bestParticles["top"] + *evh.bestParticles["antiTop"];

    cout << "Event, pT(gen_all), pT(best_all), pT(tt) = " << eventNumber << ", "
         << gen_all.Pt() << ", " << best_all.Pt() << ", " << best_tt.Pt()
         << endl;
}

void topReconstructionFromLHE::Print_smear_bs_SM(handleEvent &evh)
{
    cout << "Smeared bottom (E, px, py, pz) = ("
         << evh.smearedParticles["bottom"]->E() << ", "
         << evh.smearedParticles["bottom"]->Px() << ", "
         << evh.smearedParticles["bottom"]->Py() << ", "
         << evh.smearedParticles["bottom"]->Pz() << ")\n";

    cout << "Smeared bfromH (E, px, py, pz) = ("
         << evh.smearedParticles["bFromH"]->E() << ", "
         << evh.smearedParticles["bFromH"]->Px() << ", "
         << evh.smearedParticles["bFromH"]->Py() << ", "
         << evh.smearedParticles["bFromH"]->Pz() << ")\n";

    cout << "Smeared bbarFromH (E, px, py, pz) = ("
         << evh.smearedParticles["bbarFromH"]->E() << ", "
         << evh.smearedParticles["bbarFromH"]->Px() << ", "
         << evh.smearedParticles["bbarFromH"]->Py() << ", "
         << evh.smearedParticles["bbarFromH"]->Pz() << ")\n";
}

void topReconstructionFromLHE::Print_smear_tt_SM(handleEvent &evh)
{
    cout << "Smeared top (E, px, py, pz) = ("
         << evh.smearedParticles["top"]->E() << ", "
         << evh.smearedParticles["top"]->Px() << ", "
         << evh.smearedParticles["top"]->Py() << ", "
         << evh.smearedParticles["top"]->Pz() << ")\n";

    cout << "Smeared antiTop (E, px, py, pz) = ("
         << evh.smearedParticles["antiTop"]->E() << ", "
         << evh.smearedParticles["antiTop"]->Px() << ", "
         << evh.smearedParticles["antiTop"]->Py() << ", "
         << evh.smearedParticles["antiTop"]->Pz() << ")\n";
}

// void topReconstructionFromLHE::getBestObjects(){}

void topReconstructionFromLHE::moveStatsBox(TH1D *hist)
{
    cout << "in stats box func" << endl;
    gPad->Update();
    cout << "after gpad update" << endl;
    TPaveStats *s = (TPaveStats *)hist->FindObject("stats");
    cout << "1" << endl;
    /*(    if (s == NULL) {
            //return;
            cout <<"null pointer"<<endl;
            return;
        }*/
    // float x1 = s->GetX1NDC();
    // float x2 = s->GetX2NDC();
    float y1 = s->GetY1NDC();
    float y2 = s->GetY2NDC();
    cout << "2" << endl;
    s->SetY1NDC(y1 - (y2 - y1));
    s->SetY2NDC(y2 - (y2 - y1));
    cout << "3" << endl;
}

#endif
