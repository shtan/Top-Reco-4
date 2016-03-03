#define topReconstructionFromLHE_cxx
#include "converter.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TString.h>
#include <TH1F.h>
#include <TCanvas.h>
//#include "WDaughterEllipseCalculator.h"
//#include "topEventMinimizer.h"
#include <TPaveStats.h>
#include <TLorentzVector.h>

using namespace std;
using namespace ROOT::Math;

int main()
{

    converter t;
    // t.Loop("output_files",0,1);
    t.Loop();

    return 0;
}

void converter::Filler(int i, int j)
{
    PID.at(i) = Pid.at(j);
    P_X.at(i) = Px.at(j);
    P_Y.at(i) = Py.at(j);
    P_Z.at(i) = Pz.at(j);
    E.at(i) = En.at(j);
    M.at(i) = Ma.at(j);

    return;
}

void converter::Smearer(TLorentzVector &jet)
{
    double pt = jet.Pt();
    double phi = jet.Phi();
    double eta = jet.Eta();
    double e = jet.Energy();
    double newpt = pt + sqrt(pt) * rand.Gaus();
    double newphi = phi + 0.01 * rand.Gaus();
    jet.SetPtEtaPhiE(newpt, eta, newphi, e);
    return;
}

void converter::SmearedFiller(int i, TLorentzVector &jet, int pid)
{
    PID.at(i) = pid;
    P_X.at(i) = jet.Px();
    P_Y.at(i) = jet.Py();
    P_Z.at(i) = jet.Pz();
    E.at(i) = jet.Energy();
    M.at(i) = jet.M();
    return;
}

// void converter::Loop(TString dir, int whichLoop, int maxLoops)
void converter::Loop()
// int main()
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
    // std::cout<<"before fchain"<<std::endl;
    // if (fChain == 0) return;
    //   TRandom3 rand;
    rand.SetSeed(10);

    /*   vector<double> Px, Py, Pz, E, M;
       vector<int> PID;
       int numpar;
    */
    TFile *infile = TFile::Open("./singlemu_onshell2_6nov.root");
    TTree *intree = (TTree *)infile->Get("LHEF");

    // intree->SetMakeClass(1);

    // vector<double> Px, Py, Pz, En, Ma;
    // vector<int> Pid, mother1, mother2;
    // int numpar;
    int numparr = -99;

    // double Pxd, Pyd, Pzd, Mad;
    double_t End[200];
    double *Pxd[20];
    double *Pyd[20];
    double *Pzd[20];
    double *Mad[20];
    // End = -999;
    int *Pidd[20];
    int *mother1d[20];
    int *mother2d[20];
    *(End)[0] = 14;
    *(End)[1] = 16;
    *(End)[10] = 57;
    *(Pxd)[0] = 15;

    intree->SetMakeClass(1);

    intree->SetBranchAddress("Particle.Px", Pxd);
    intree->SetBranchAddress("Particle.Py", Pyd);
    intree->SetBranchAddress("Particle.Pz", Pzd);
    intree->SetBranchAddress("Particle.E", End);
    intree->SetBranchAddress("Particle.M", Mad);
    intree->SetBranchAddress("Particle.PID", Pidd);
    intree->SetBranchAddress("Particle_size", &numparr);
    intree->SetBranchAddress("Particle.Mother1", mother1d);
    intree->SetBranchAddress("Particle.Mother2", mother2d);
    /*
       intree->SetBranchAddress("Particle.Px", &Px);
       intree->SetBranchAddress("Particle.Py", &Py);
       intree->SetBranchAddress("Particle.Pz", &Pz);
       intree->SetBranchAddress("Particle.E", &En);
       intree->SetBranchAddress("Particle.M", &Ma);
       intree->SetBranchAddress("Particle.PID", &Pid);
       intree->SetBranchAddress("Particle.@size",&numparr);
       intree->SetBranchAddress("Particle.Mother1", &mother1);
       intree->SetBranchAddress("Particle.Mother2", &mother2);
    */
    TFile *outfile =
        new TFile("./fixedsinglemuntuple.root", "RECREATE", "tree");
    TTree *outtree = new TTree("outtree", "outtree");

    // vector<double> P_X, P_Y, P_Z, E, M;
    // vector<int> PID;

    outtree->Branch("P_X", &P_X);
    outtree->Branch("P_Y", &P_Y);
    outtree->Branch("P_Z", &P_Z);
    outtree->Branch("E", &E);
    outtree->Branch("M", &M);
    outtree->Branch("PID", &PID);

    for (int i = 0; i < 20; i++) {
        P_X.push_back(0);
        P_Y.push_back(0);
        P_Z.push_back(0);
        E.push_back(0);
        M.push_back(0);
        PID.push_back(-9999);
    }

    // initOutput(dir,whichLoop);

    // Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nentries = intree->GetEntries();

    // Long64_t nbytes = 0, nb = 0;

    vector<XYZTLorentzVector> nonTopObjects;
    vector<double> nonTopObjectPtWidths;
    vector<double> nonTopObjectEtaWidths;
    vector<double> nonTopObjectPhiWidths;

    // gStyle->SetOptStat(101110);
    // TCanvas c1("c1","c1",500,500);
    // TH1D h_chi2("h_chi2","",100,0,100);
    std::cout << "hi!" << std::endl;

    // int jStart = whichLoop*(nentries/maxLoops) +
    // ((whichLoop>(maxLoops-nentries%maxLoops))?(whichLoop+nentries%maxLoops-maxLoops):0);
    // int jStart = 17;
    // int jFinish = jStart + (nentries+whichLoop)/maxLoops;
    std::cout << "number of entries = " << nentries << std::endl;
    int jStart = 0;
    // int jFinish = nentries;
    int jFinish = 5;

    std::cout << jStart << std::endl;
    std::cout << "nentries = " << nentries << std::endl;
    std::cout << jFinish << std::endl;

    // TH1F* bottomPt = new TH1F("bottomPt","bottomPt",100,-300,300);

    for (Long64_t jentry = jStart; jentry < jFinish; jentry++) {
        cout << "BEGINNING BRANCH NUMBER " << jentry << endl;
        // Long64_t ientry = LoadTree(jentry);
        intree->GetEntry(jentry);

        cout << "numparr = " << numparr << endl;

        cout << *(End)[0] << endl;
        cout << *(End)[1] << endl;
        cout << *(End)[10] << endl;
        cout << *(Pxd)[0] << endl;

        numpar = (int)(Px.size());

        if (numpar != 15) {
            cout << "numpar != 15" << endl;
            cout << numpar << endl;
            continue;
        }
        // if (ientry < 0) break;
        // if (jentry > 0) break;
        // nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;
        //

        TLorentzVector bJetFromTop, bbarJetFromTop, bJetFromH, bbarJetFromH,
            qFromW, qbarFromW, muon;
        double leptonPx, leptonPy;
        int qFromWPID, qbarFromWPID;

        int counter = 0;

        for (int ivec = 0; ivec < numpar; ivec++) {
            if (Pid.at(ivec) == 5) {
                if (Pid.at(mother1.at(ivec)) == 6 or
                    Pid.at(mother2.at(ivec)) == 6) {
                    Filler(8, ivec);
                    counter++;
                    bJetFromTop.SetPxPyPzE(Px.at(ivec), Py.at(ivec),
                                           Pz.at(ivec), E.at(ivec));
                }
                if (Pid.at(mother1.at(ivec)) == 25 or
                    Pid.at(mother2.at(ivec)) == 25) {
                    Filler(12, ivec);
                    counter++;
                    bJetFromH.SetPxPyPzE(Px.at(ivec), Py.at(ivec), Pz.at(ivec),
                                         E.at(ivec));
                }
            } else if (Pid.at(ivec) == -5) {
                if (Pid.at(mother1.at(ivec)) == -6 or
                    Pid.at(mother2.at(ivec)) == -6) {
                    Filler(9, ivec);
                    counter++;
                    bbarJetFromTop.SetPxPyPzE(Px.at(ivec), Py.at(ivec),
                                              Pz.at(ivec), E.at(ivec));
                }
                if (Pid.at(mother1.at(ivec)) == 25 or
                    Pid.at(mother2.at(ivec)) == 25) {
                    Filler(13, ivec);
                    counter++;
                    bbarJetFromH.SetPxPyPzE(Px.at(ivec), Py.at(ivec),
                                            Pz.at(ivec), E.at(ivec));
                }
            } else if (Pid.at(ivec) == 13) {
                Filler(7, ivec);
                counter++;
                muon.SetPxPyPzE(Px.at(ivec), Py.at(ivec), Pz.at(ivec),
                                E.at(ivec));
                SmearedFiller(0, muon, 13);
                leptonPx = Px.at(ivec);
                leptonPy = Py.at(ivec);
            } else if (Pid.at(ivec) == 1 or Pid.at(ivec) == 2 or
                       Pid.at(ivec) == 3 or Pid.at(ivec) == 4) {
                Filler(10, ivec);
                counter++;
                qFromW.SetPxPyPzE(Px.at(ivec), Py.at(ivec), Pz.at(ivec),
                                  E.at(ivec));
                qFromWPID = Pid.at(ivec);
            } else if (Pid.at(ivec) == -1 or Pid.at(ivec) == -2 or
                       Pid.at(ivec) == -3 or Pid.at(ivec) == -4) {
                Filler(11, ivec);
                counter++;
                qbarFromW.SetPxPyPzE(Px.at(ivec), Py.at(ivec), Pz.at(ivec),
                                     E.at(ivec));
                qbarFromWPID = Pid.at(ivec);
            } else if (Pid.at(ivec) == -14) {
                Filler(18, ivec);
                counter++;
            } else if (Pid.at(ivec) == 6) {
                Filler(14, ivec);
                counter++;
            } else if (Pid.at(ivec) == -6) {
                Filler(15, ivec);
            } else if (Pid.at(ivec) == 24) {
                Filler(16, ivec);
                counter++;
            } else if (Pid.at(ivec) == -24) {
                Filler(17, ivec);
                counter++;
            } else if (Pid.at(ivec) == 25) {
                Filler(19, ivec);
                counter++;
            }
        }
        if (counter != 13) {
            cout << "counter != 13" << endl;
            continue;
        }

        // start smearing
        // double newpt = 0;
        Smearer(bJetFromTop);
        Smearer(bbarJetFromTop);
        Smearer(qFromW);
        Smearer(qbarFromW);
        Smearer(bJetFromH);
        Smearer(bbarJetFromH);

        SmearedFiller(1, bJetFromTop, 5);
        SmearedFiller(2, bbarJetFromTop, -5);
        SmearedFiller(3, qFromW, qFromWPID);
        SmearedFiller(4, qbarFromW, qbarFromWPID);
        SmearedFiller(5, bJetFromH, 5);
        SmearedFiller(6, bbarJetFromH, -5);

        double metPx =
            -(bJetFromTop.Px() + bbarJetFromTop.Px() + bJetFromH.Px() +
              bbarJetFromH.Px() + qFromW.Px() + qbarFromW.Px() + leptonPx);
        double metPy =
            -(bJetFromTop.Py() + bbarJetFromTop.Py() + bJetFromH.Py() +
              bbarJetFromH.Py() + qFromW.Py() + qbarFromW.Py() + leptonPy);

        P_X.at(20) = metPx;
        P_Y.at(20) = metPy;
        P_Z.at(20) = 0;
        PID.at(20) = -9999;
        E.at(20) = 0;
        M.at(20) = 0;
    }
    outtree->Fill();
    outfile->Write();
    outfile->Close();
    return;
}

/*
      XYZTLorentzVector lepton, antiLepton;
      XYZTLorentzVector neutrino, antiNeutrino;
      XYZTLorentzVector bottomQuark, antiBottomQuark;
      XYZTLorentzVector lightParton1, lightParton2;
      XYZTLorentzVector topQuark, antiTopQuark;
      XYZTLorentzVector recoTopQuark(0,0,0,0), recoAntiTopQuark(0,0,0,0);
      XYZTLorentzVector WPlus, WMinus, higgs, qFromW, qbarFromW, met;
      XYZTLorentzVector recoWPlus(0,0,0,0), recoWMinus(0,0,0,0);
      XYZTLorentzVector smearedLightQuark1, smearedLightQuark2,
smearedBottomQuark, smearedAntiBottomQuark;
      XYZTLorentzVector smearedLightParton1, smearedLightParton2, smearedLepton,
smearedAntiLepton;
      XYZTLorentzVector smearedOtherLightParton;
      std::vector<XYZTLorentzVector> smearedOtherLightPartons;
      bool leptonFlag = 0; //1 if lepton and antineutrino (i.e. tbar branch is
the leptonic one), 0 otherwise.

      double METx(0.), METy(0.);
      int iJet = 0;
//      for(int iParticle = 0 ; iParticle < n_particles; iParticle++)
//	{
      smearedLightQuark1.SetPxPyPzE(P_X->at(3),P_Y->at(3),P_Z->at(3),E->at(3));
      smearedLightQuark2.SetPxPyPzE(P_X->at(4),P_Y->at(4),P_Z->at(4),E->at(4));
      smearedBottomQuark.SetPxPyPzE(P_X->at(1),P_Y->at(1),P_Z->at(1),E->at(1));
      smearedAntiBottomQuark.SetPxPyPzE(P_X->at(2),P_Y->at(2),P_Z->at(2),E->at(2));
      smearedLightParton1.SetPxPyPzE(P_X->at(5),P_Y->at(5),P_Z->at(5),E->at(5));
      smearedLightParton2.SetPxPyPzE(P_X->at(6),P_Y->at(6),P_Z->at(6),E->at(6));
      //smearedLepton.SetPxPyPzE(P_X->at(0),P_Y->at(0),P_Z->at(0),E->at(0));

      topQuark.SetPxPyPzE(P_X->at(14),P_Y->at(14),P_Z->at(14),E->at(14));
      antiTopQuark.SetPxPyPzE(P_X->at(15),P_Y->at(15),P_Z->at(15),E->at(15));
      WPlus.SetPxPyPzE(P_X->at(16),P_Y->at(16),P_Z->at(16),E->at(16));
      WMinus.SetPxPyPzE(P_X->at(17),P_Y->at(17),P_Z->at(17),E->at(17));
          //lepton.SetPxPyPzE(P_X->at(7),P_Y->at(7),P_Z->at(7),E->at(7));
      //neutrino.SetPxPyPzE(P_X->at(18),P_Y->at(18),P_Z->at(18),E->at(18));
      bottomQuark.SetPxPyPzE(P_X->at(8),P_Y->at(8),P_Z->at(8),E->at(8));
      antiBottomQuark.SetPxPyPzE(P_X->at(9),P_Y->at(9),P_Z->at(9),E->at(9));
      qFromW.SetPxPyPzE(P_X->at(10),P_Y->at(10),P_Z->at(10),E->at(10));
      qbarFromW.SetPxPyPzE(P_X->at(11),P_Y->at(11),P_Z->at(11),E->at(11));
      lightParton1.SetPxPyPzE(P_X->at(12),P_Y->at(12),P_Z->at(12),E->at(12));
      lightParton2.SetPxPyPzE(P_X->at(13),P_Y->at(13),P_Z->at(13),E->at(13));
      higgs.SetPxPyPzE(P_X->at(19),P_Y->at(19),P_Z->at(19),E->at(19));
          met.SetPxPyPzE(P_X->at(20),P_Y->at(20),P_Z->at(20),E->at(20));
          std::cout<<"after set met"<<std::endl;

        XYZTLorentzVector allExternal;
        allExternal.SetPxPyPzE(0,0,0,0);

          for ( int ii=0; ii < ( (int)PID->size() - 21 ); ii++){
              smearedOtherLightParton.SetPxPyPzE(P_X->at(ii+21), P_Y->at(ii+21),
P_Z->at(ii+21), E->at(ii+21) );
              smearedOtherLightPartons.push_back(smearedOtherLightParton);

              allExternal += smearedOtherLightParton;
          }
          std::cout<<"SIZEEEEEEEEEEEEEEEEEEE = " <<
smearedOtherLightPartons.size() << std::endl;

        XYZTLorentzVector leptonToAdd;

          if (PID->at(0) == 13) {
              smearedLepton.SetPxPyPzE(P_X->at(0),P_Y->at(0),P_Z->at(0),E->at(0));
                  lepton.SetPxPyPzE(P_X->at(7),P_Y->at(7),P_Z->at(7),E->at(7));
              antiNeutrino.SetPxPyPzE(P_X->at(18),P_Y->at(18),P_Z->at(18),E->at(18));
                  leptonFlag = 1;
                  leptonToAdd = smearedLepton;
          } else if (PID->at(0) == -13){
              smearedAntiLepton.SetPxPyPzE(P_X->at(0),P_Y->at(0),P_Z->at(0),E->at(0));
                  antiLepton.SetPxPyPzE(P_X->at(7),P_Y->at(7),P_Z->at(7),E->at(7));
              neutrino.SetPxPyPzE(P_X->at(18),P_Y->at(18),P_Z->at(18),E->at(18));
                  leptonFlag = 0;
                  leptonToAdd = smearedAntiLepton;
          }

    XYZTLorentzVector allVisible = smearedLightQuark1 + smearedLightQuark2 +
smearedBottomQuark + smearedAntiBottomQuark + smearedLightParton1 +
smearedLightParton2 + leptonToAdd + allExternal;

    std::cout<< "comparison "<<std::endl;
    std::cout<< allVisible.Px() <<std::endl;
    std::cout<< met.Px() << std::endl;
    //std::cout<<met.Py()<<std::endl;
    //std::cout<<met.Pz()<<std::endl;
    //std::cout<<met.E()<<std::endl;

    //hack-test to form an extra particle to absorb the discrepancy between met
and sum of all visible particles
    XYZTLorentzVector extraParticle;
    extraParticle.SetPxPyPzE(-allVisible.Px()-met.Px(),
-allVisible.Py()-met.Py(), 0, pow( (allVisible.Px()+met.Px()),2) + pow(
(allVisible.Py()+met.Py()),2) );

    std::cout<< (allVisible.Px()+extraParticle.Px()) <<std::endl;
    std::cout<< met.Px() << std::endl;



          //	} //end loop over particles

      //cout << "true W+             mass is " << WPlus.M()        << endl;
      //cout << "true top quark      mass is " << topQuark.M()     << endl;
      //cout << "true W-             mass is " << WMinus.M()       << endl;
      //cout << "true anti top quark mass is " << antiTopQuark.M() << endl;


      //cout << "true neutrino momentum:" << endl;
      //printVector(neutrino);
      //cout << "true anti-neutrino momentum:" << endl;
      //printVector(antiNeutrino);


      //RESET EVERYTHING
      recoWMinus.SetPxPyPzE(0,0,0,0);
      recoAntiTopQuark.SetPxPyPzE(0,0,0,0);
      recoWPlus.SetPxPyPzE(0,0,0,0);
      recoTopQuark.SetPxPyPzE(0,0,0,0);
      METx = 0.;
      METy = 0.;
      nonTopObjects.clear();
      nonTopObjectPtWidths.clear();
      nonTopObjectEtaWidths.clear();
      nonTopObjectPhiWidths.clear();


      //for fully leptonic decays
      //METx -= lepton.Px() + antiLepton.Px();
      //METy -= lepton.Py() + antiLepton.Py();
      //simulate hadronically decaying \bar{t} quark, so smear the lepton and
subtract it from the MET separately
      METx -= antiLepton.Px();
      METy -= antiLepton.Py();
*/

/*      XYZTLorentzVector smearedBottomQuark;
      XYZTLorentzVector smearedAntiBottomQuark;
      XYZTLorentzVector smearedLightQuark1;
      XYZTLorentzVector smearedLightQuark2;
      XYZTLorentzVector smearedLightParton1;
      XYZTLorentzVector smearedLightParton2;
      XYZTLorentzVector tempP4;
*/
/*      double bJet1PtSmear = log(1.+0.1)*rand.Gaus();
      double bJet2PtSmear = log(1.+0.1)*rand.Gaus();
      double lightJet1PtSmear = log(1.+0.1)*rand.Gaus();
      double lightJet2PtSmear = log(1.+0.1)*rand.Gaus();
      double bJet1PhiSmear = 0.01*rand.Gaus();
      double bJet2PhiSmear = 0.01*rand.Gaus();
      double lightJet1PhiSmear = 0.01*rand.Gaus();
      double lightJet2PhiSmear = 0.01*rand.Gaus();
      double bJet1EtaSmear = 0.01*rand.Gaus();
      double bJet2EtaSmear = 0.01*rand.Gaus();
      //smear the lepton and anti-neutrino as if they were jets, to simulate the
   \bar{t} quark decaying hadronically
      double leptonPtSmear = log(1.+0.1)*rand.Gaus();
      double leptonPhiSmear = 0.01*rand.Gaus();
      double leptonEtaSmear = 0.01*rand.Gaus();
      double antiNeutrinoPtSmear = log(1.+0.1)*rand.Gaus();
      double antiNeutrinoPhiSmear = 0.01*rand.Gaus();
      double antiNeutrinoEtaSmear = 0.01*rand.Gaus();
*/

/*
      double pt(bottomQuark.Pt()*exp(bJet1PtSmear));
      double phi(bottomQuark.Phi()+bJet1PhiSmear);
      double eta(bottomQuark.Eta()+bJet1EtaSmear);
      double
   energy(sqrt((bottomQuark.M2()+bottomQuark.Pt()*bottomQuark.Pt()*cosh(eta)*cosh(eta)))*exp(bJet1PtSmear));

      smearedBottomQuark.SetPx(pt*cos(phi));
      smearedBottomQuark.SetPy(pt*sin(phi));
      smearedBottomQuark.SetPz(pt*sinh(eta));
      smearedBottomQuark.SetE(energy);
*/
//      recoTopQuark+=smearedBottomQuark;
/*
      METx -= smearedBottomQuark.Px();
      METy -= smearedBottomQuark.Py();

      pt = antiBottomQuark.Pt()*exp(bJet2PtSmear);
      phi = antiBottomQuark.Phi()+bJet2PhiSmear;
      eta = antiBottomQuark.Eta()+bJet2EtaSmear;
      energy =
   sqrt((antiBottomQuark.M2()+antiBottomQuark.Pt()*antiBottomQuark.Pt()*cosh(eta)*cosh(eta)))*exp(bJet2PtSmear);

      smearedAntiBottomQuark.SetPx(pt*cos(phi));
      smearedAntiBottomQuark.SetPy(pt*sin(phi));
      smearedAntiBottomQuark.SetPz(pt*sinh(eta));
      smearedAntiBottomQuark.SetE(energy);
*/
//      recoAntiTopQuark+=smearedAntiBottomQuark;
/*
      METx -= smearedAntiBottomQuark.Px();
      METy -= smearedAntiBottomQuark.Py();

      //smear the lepton to simulate a hadronic \bar{t} decay
      pt = lepton.Pt()*exp(leptonPtSmear);
      phi = lepton.Phi()+leptonPhiSmear;
      eta = lepton.Eta()+leptonEtaSmear;
      energy =
   sqrt((lepton.M2()+lepton.Pt()*lepton.Pt()*cosh(eta)*cosh(eta)))*exp(leptonPtSmear);

      smearedLightQuark1.SetPx(pt*cos(phi));
      smearedLightQuark1.SetPy(pt*sin(phi));
      smearedLightQuark1.SetPz(pt*sinh(eta));
      smearedLightQuark1.SetE(energy);

      METx -= smearedLightQuark1.Px();
      METy -= smearedLightQuark1.Py();

      //also smear the anti-neutrino
      pt = antiNeutrino.Pt()*exp(antiNeutrinoPtSmear);
      phi = antiNeutrino.Phi()+antiNeutrinoPhiSmear;
      eta = antiNeutrino.Eta()+antiNeutrinoEtaSmear;
      energy =
   sqrt((antiNeutrino.M2()+antiNeutrino.Pt()*antiNeutrino.Pt()*cosh(eta)*cosh(eta)))*exp(antiNeutrinoPtSmear);

      smearedLightQuark2.SetPx(pt*cos(phi));
      smearedLightQuark2.SetPy(pt*sin(phi));
      smearedLightQuark2.SetPz(pt*sinh(eta));
      smearedLightQuark2.SetE(energy);

      METx -= smearedLightQuark2.Px();
      METy -= smearedLightQuark2.Py();


      pt = lightParton1.Pt()*exp(lightJet1PtSmear);
      phi = lightParton1.Phi()+lightJet1PhiSmear;
      eta = lightParton1.Eta();
      energy = lightParton1.E()*exp(lightJet1PtSmear);

      smearedLightParton1.SetPx(pt*cos(phi));
      smearedLightParton1.SetPy(pt*sin(phi));
      smearedLightParton1.SetPz(pt*sinh(eta));
      smearedLightParton1.SetE(energy);
      METx -= smearedLightParton1.Px();
      METy -= smearedLightParton1.Py();


      pt = lightParton2.Pt()*exp(lightJet2PtSmear);
      phi = lightParton2.Phi()+lightJet2PhiSmear;
      eta = lightParton2.Eta();
      energy = lightParton2.E()*exp(lightJet2PtSmear);

      smearedLightParton2.SetPx(pt*cos(phi));
      smearedLightParton2.SetPy(pt*sin(phi));
      smearedLightParton2.SetPz(pt*sinh(eta));
      smearedLightParton2.SetE(energy);
      METx -= smearedLightParton2.Px();
      METy -= smearedLightParton2.Py();

*/

//      if(smearedLightParton2.Pt()>smearedLightParton1.Pt())
//	{
//	  //tempP4 = smearedLightParton1;
//	  //smearedLightParton1 = smearedLightParton2;
//	  //smearedLightParton2 = tempP4;
//	  smearingSwitchedLightJetOrdering = true;
//	}
//      else
//	{
//	  smearingSwitchedLightJetOrdering = false;
//	}

/*

      nonTopObjects.push_back(smearedLightParton1);
      nonTopObjectPtWidths.push_back(sqrt(smearedLightParton1.pt()));
      nonTopObjectEtaWidths.push_back(0.01);
      nonTopObjectPhiWidths.push_back(0.01);

      nonTopObjects.push_back(smearedLightParton2);
      nonTopObjectPtWidths.push_back(sqrt(smearedLightParton2.pt()));
      nonTopObjectEtaWidths.push_back(0.01);
      nonTopObjectPhiWidths.push_back(0.01);

      for (int ii=0; ii<(int)smearedOtherLightPartons.size(); ii++){
        nonTopObjects.push_back(smearedOtherLightPartons.at(ii));
        nonTopObjectPtWidths.push_back( sqrt( (smearedOtherLightPartons.at(ii)).pt() ) );
//        nonTopObjectPtWidths.push_back(0.00001);
        nonTopObjectEtaWidths.push_back(0.01);
        nonTopObjectPhiWidths.push_back(0.01);
      }

      //hack test of extra particle
      nonTopObjects.push_back(extraParticle);
      nonTopObjectPtWidths.push_back( sqrt(extraParticle.pt()) );
      nonTopObjectEtaWidths.push_back(0.01);
      nonTopObjectPhiWidths.push_back(0.01);

      std::cout<<"SIZEEEEEEEEEEEEEEEEEEEEEEE = "<<nonTopObjects.size() << std::endl;

          


      if (leptonFlag == 0){
      recoWPlus += smearedAntiLepton;
      //recoWPlus += smearedNeutrino;
      recoWMinus += smearedLightQuark1;
      recoWMinus += smearedLightQuark2;
      } else if (leptonFlag == 1){
      recoWMinus += smearedLepton;
      //recoWMinus += smearedAntiNeutrino;
      recoWPlus += smearedLightQuark1;
      recoWPlus += smearedLightQuark2;
      }

      recoTopQuark+=recoWPlus;
      recoAntiTopQuark+=recoWMinus;



      //cout << "METx is " << METx << endl;				     
      //cout << "real METx is " << neutrino.Px() + antiNeutrino.Px() << endl;
      //cout << "METy is " << METy << endl;				     
      //cout << "real METy is " << neutrino.Py() + antiNeutrino.Py() << endl;
      //cout << "neutrino       mass is      " << neutrino.M() << endl;
      //printVector(neutrino);
      //cout << "anti-lepton    mass is      " << antiLepton.M() << endl;
      //printVector(antiLepton);
      //cout << "reco W+        mass is      " << recoWPlus.M() << endl;
      //printVector(recoWPlus);
      //cout << "b quark        mass is      " << smearedBottomQuark.M() << endl;
      //printVector(smearedBottomQuark);
      //cout << "reco top quark mass is      " << recoTopQuark.M() << endl;
      //printVector(recoTopQuark);
      //cout << "anti-neutrino  mass is      " << antiNeutrino.M() << endl;
      //printVector(antiNeutrino);
      //cout << "lepton         mass is      " << lepton.M() << endl;
      //printVector(lepton);
      //cout << "smeared   1st light quark mass is " << smearedLightQuark1.M() << endl;
      //printVector(smearedLightQuark1);
      //cout << "smeared   2nd light quark mass is " << smearedLightQuark2.M() << endl;
      //printVector(smearedLightQuark2);
      //cout << "reco W-        mass is      " << recoWMinus.M() << endl;
      //printVector(recoWMinus);
      //cout << "b-bar quark    mass is      " << smearedAntiBottomQuark.M() << endl;
      //printVector(smearedAntiBottomQuark);
      //cout << "reco anti top quark mass is " << recoAntiTopQuark.M() << endl;
      //printVector(recoAntiTopQuark);
      //cout << "first light parton  mass is " << smearedLightParton1.M() << endl;
      //printVector(smearedLightParton1);
      //cout << "second light parton mass is " << smearedLightParton2.M() << endl;
      //printVector(smearedLightParton2);


	  

      //if(jentry!=1&&jentry!=4&&jentry!=7) continue;
      //if(jentry!=8) continue;
      


      double bJet2Px = smearedAntiBottomQuark.Px();
      double bJet2Py = smearedAntiBottomQuark.Py();
      double bJet2Pz = smearedAntiBottomQuark.Pz();
      double bJet2E  = smearedAntiBottomQuark.E() ;
      double bJet2Pt = smearedAntiBottomQuark.Pt();

      double bJet1Px = smearedBottomQuark.Px();
      double bJet1Py = smearedBottomQuark.Py();
      double bJet1Pz = smearedBottomQuark.Pz();
      double bJet1E  = smearedBottomQuark.E() ;
      double bJet1Pt = smearedBottomQuark.Pt();

      double lightJet1Px = smearedLightQuark1.Px();
      double lightJet1Py = smearedLightQuark1.Py();
      double lightJet1Pz = smearedLightQuark1.Pz();
      double lightJet1E  = smearedLightQuark1.E() ;
      double lightJet1Pt = smearedLightQuark1.Pt();

      double lightJet2Px = smearedLightQuark2.Px();
      double lightJet2Py = smearedLightQuark2.Py();
      double lightJet2Pz = smearedLightQuark2.Pz();
      double lightJet2E  = smearedLightQuark2.E() ;
      double lightJet2Pt = smearedLightQuark2.Pt();

      double antiLeptonPx, antiLeptonPy, antiLeptonPz, antiLeptonE, antiLeptonPt;
      if (leptonFlag == 0){
       antiLeptonPx = smearedAntiLepton.Px();
       antiLeptonPy = smearedAntiLepton.Py();
       antiLeptonPz = smearedAntiLepton.Pz();
       antiLeptonE  = smearedAntiLepton.E() ;
       antiLeptonPt = smearedAntiLepton.Pt();
      }

      double leptonPx, leptonPy, leptonPz, leptonE, leptonPt;
      if (leptonFlag == 1){
       leptonPx = smearedLepton.Px();
       leptonPy = smearedLepton.Py();
       leptonPz = smearedLepton.Pz();
       leptonE  = smearedLepton.E() ;
       leptonPt = smearedLepton.Pt();
      }

      double mTop=173., mW=80.4;
      double sigmaMTop=2.0, sigmaMW=2.085;
//      double sigmaPtJet=0.1, sigmaPhiJet=0., sigmaEtaJet=0.;
      double sigmaPtLep=0. , sigmaPhiLep=0.  , sigmaEtaLep=0.  ;
      double sigmaPtJet=0.1, sigmaPhiJet=0.01, sigmaEtaJet=0.01;
//      double sigmaPtLep=0.1 , sigmaPhiLep=0.01  , sigmaEtaLep=0.01  ;

      topEventMinimizer ev(nonTopObjects,
			   nonTopObjectPtWidths,
			   nonTopObjectEtaWidths,
			   nonTopObjectPhiWidths,
			   mTop,sigmaMTop,
			   mW,sigmaMW
			   );

      std::cout<<"before add tops"<<std::endl;
      if (leptonFlag == 0){
      ev.addLeptonicTop(bJet1Px,bJet1Py,bJet1Pz,bJet1E,
                        //sigmaPtJet,sigmaEtaJet,sigmaPhiJet,
                        sqrt(bJet1Pt),sigmaEtaJet,sigmaPhiJet,
                        antiLeptonPx,antiLeptonPy,antiLeptonPz,antiLeptonE,
                        sigmaPtLep,sigmaEtaLep,sigmaPhiLep,
                        mTop,sigmaMTop,
                        mW,sigmaMW);
      
      ev.addHadronicTop(bJet2Px,bJet2Py,bJet2Pz,bJet2E,
			sqrt(bJet2Pt),sigmaEtaJet,sigmaPhiJet,
			//sigmaPtJet,sigmaEtaJet,sigmaPhiJet,
			//0.000001, sigmaEtaJet,sigmaPhiJet,
                        lightJet1Px,lightJet1Py,lightJet1Pz,lightJet1E,
			//lightJet2Px,lightJet2Py,lightJet2Pz,lightJet2E,
			sqrt(lightJet1Pt),sigmaEtaJet,sigmaPhiJet,
			//sigmaPtJet,sigmaEtaJet,sigmaPhiJet,
			lightJet2Px,lightJet2Py,lightJet2Pz,lightJet2E,
			//lightJet1Px,lightJet1Py,lightJet1Pz,lightJet1E,
			//sigmaPtJet,sigmaEtaJet,sigmaPhiJet,
		        sqrt(lightJet2Pt), 0.01, 0.01,
                        mTop,sigmaMTop,
			mW,sigmaMW);
      }


      if (leptonFlag == 1){
      ev.addHadronicTop(bJet1Px,bJet1Py,bJet1Pz,bJet1E,
			sqrt(bJet1Pt),sigmaEtaJet,sigmaPhiJet,
			//sigmaPtJet,sigmaEtaJet,sigmaPhiJet,
			//0.000001, sigmaEtaJet,sigmaPhiJet,
			lightJet1Px,lightJet1Py,lightJet1Pz,lightJet1E,
			//lightJet2Px,lightJet2Py,lightJet2Pz,lightJet2E,
			sqrt(lightJet1Pt),sigmaEtaJet,sigmaPhiJet,
			//sigmaPtJet,sigmaEtaJet,sigmaPhiJet,
			lightJet2Px,lightJet2Py,lightJet2Pz,lightJet2E,
			//lightJet1Px,lightJet1Py,lightJet1Pz,lightJet1E,
			//sigmaPtJet,sigmaEtaJet,sigmaPhiJet,
			sqrt(lightJet2Pt), 0.01, 0.01,
                        mTop,sigmaMTop,
			mW,sigmaMW);

      ev.addLeptonicTop(bJet2Px,bJet2Py,bJet2Pz,bJet2E,
                        //sigmaPtJet,sigmaEtaJet,sigmaPhiJet,
                        sqrt(bJet2Pt),sigmaEtaJet,sigmaPhiJet,
                        leptonPx,leptonPy,leptonPz,leptonE,
                        sigmaPtLep,sigmaEtaLep,sigmaPhiLep,
                        mTop,sigmaMTop,
                        mW,sigmaMW);
      }

      std::cout<<"after add tops"<<std::endl;
      //ev.printTopConstituents();
      ev.initializeDeltas();
      //ev.calcTopMassRanges();
      std::cout<<"after initialise deltas"<<std::endl;

      //ev.findStartingValues(50);
      //TString thisPlotName="ellipses_";
      //thisPlotName+=jentry;
      //ev.plotEllipses(thisPlotName);

      //ev.minimizeNonTopChiSquare();

      ev.minimizeTotalChiSquare();
      std::cout<<"after minimise"<<std::endl;
      
      //h_chi2.Fill(ev.getChiSquare());      


      //Fill branches

      eventNumber=jentry;

      //Minimizer stuff
       nonTopChi2=ev.getNonTopChiSquare();
       hadronicChi2=ev.getHadronicChiSquare();
       topMassChi2=ev.getTopMassChiSquare();
       topSystemChi2=ev.getTopChiSquare();
       totalChi2=ev.getChiSquare();

      //cout << "Total chi^2 is " << totalChi2_ << endl;
      //cout << "Hadronic chi^2 is " << hadronicChi2_ << endl;
      //cout << "Top mass chi^2 is " << topMassChi2_ << endl;
      //cout << "Non-top chi^2 is " << nonTopChi2_ << endl;
      //cout << "Top system chi^2 is " << topSystemChi2_ << endl;

       innerMinStatus=ev.getInnerMinimizerStatus();
       outerMinStatus=ev.getOuterMinimizerStatus();
       outerMinEdm=ev.getOuterMinimizerEdm();

      //True values
       top1TruePt        =topQuark.Pt()    ; 
       top1TrueEta       =topQuark.Eta()   ; 
       top1TruePhi       =topQuark.Phi()   ;
       top1TrueE         =topQuark.E()     ;
       top1TrueMass      =sqrt(topQuark.M2()); 
       W1TruePt          =WPlus.Pt()       ;
       W1TrueEta         =WPlus.Eta()      ;
       W1TruePhi         =WPlus.Phi()      ;
       W1TrueE           =WPlus.E()        ;
       W1TrueMass        =sqrt(WPlus.M2()) ;
       bJet1TruePt       =bottomQuark.Pt() ;
       bJet1TrueEta      =bottomQuark.Eta();
       bJet1TruePhi      =bottomQuark.Phi();
       bJet1TrueE        =bottomQuark.E()  ;
      if (leptonFlag == 0){
       W1Daughter1TruePt =antiLepton.Pt()  ;
       W1Daughter1TrueEta=antiLepton.Eta() ;
       W1Daughter1TruePhi=antiLepton.Phi() ;
       W1Daughter1TrueE  =antiLepton.E()   ;
       W1Daughter2TruePt =neutrino.Pt()    ;
       W1Daughter2TrueEta=neutrino.Eta()   ;
       W1Daughter2TruePhi=neutrino.Phi()   ;
       W1Daughter2TrueE  =neutrino.E()     ;
      }
      else if (leptonFlag == 1){
       W1Daughter1TruePt =qFromW.Pt()  ;
       W1Daughter1TrueEta=qFromW.Eta() ;
       W1Daughter1TruePhi=qFromW.Phi() ;
       W1Daughter1TrueE  =qFromW.E()   ;
       W1Daughter2TruePt =qbarFromW.Pt()    ;
       W1Daughter2TrueEta=qbarFromW.Eta()   ;
       W1Daughter2TruePhi=qbarFromW.Phi()   ;
       W1Daughter2TrueE  =qbarFromW.E()     ;
      }
       lightJet1TruePx   =lightParton1.Px(); 
       lightJet1TruePy   =lightParton1.Py();
       lightJet1TruePz   =lightParton1.Pz();
       lightJet1TrueE    =lightParton1.E() ;
       lightJet1TruePt   =lightParton1.Pt();


       top2TruePt        =antiTopQuark.Pt()    ;
       top2TrueEta       =antiTopQuark.Eta()   ;
       top2TruePhi       =antiTopQuark.Phi()   ;
       top2TrueE         =antiTopQuark.E()     ;
       top2TrueMass      =sqrt(antiTopQuark.M2());
       W2TruePt          =WMinus.Pt()          ;
       W2TrueEta         =WMinus.Eta()         ;
       W2TruePhi         =WMinus.Phi()         ;
       W2TrueE           =WMinus.E()           ;
       W2TrueMass        =sqrt(WMinus.M2())    ;
       bJet2TruePt       =antiBottomQuark.Pt() ;
       bJet2TrueEta      =antiBottomQuark.Eta();
       bJet2TruePhi      =antiBottomQuark.Phi();
       bJet2TrueE        =antiBottomQuark.E()  ;
      if (leptonFlag == 1){
       W2Daughter2TruePt =antiNeutrino.Pt()    ;
       W2Daughter2TrueEta=antiNeutrino.Eta()   ;
       W2Daughter2TruePhi=antiNeutrino.Phi()   ;
       W2Daughter2TrueE  =antiNeutrino.E()     ;
       W2Daughter1TruePt =lepton.Pt()          ;
       W2Daughter1TrueEta=lepton.Eta()         ;
       W2Daughter1TruePhi=lepton.Phi()         ;
       W2Daughter1TrueE  =lepton.E()           ;
      } else if (leptonFlag == 0){
       W2Daughter1TruePt =qFromW.Pt()    ;
       W2Daughter1TrueEta=qFromW.Eta()   ;
       W2Daughter1TruePhi=qFromW.Phi()   ;
       W2Daughter1TrueE  =qFromW.E()     ;
       W2Daughter2TruePt =qbarFromW.Pt()          ;
       W2Daughter2TrueEta=qbarFromW.Eta()         ;
       W2Daughter2TruePhi=qbarFromW.Phi()         ;
       W2Daughter2TrueE  =qbarFromW.E()           ;
      }
       lightJet2TruePx   =lightParton2.Px()    ;
       lightJet2TruePy   =lightParton2.Py()    ;
       lightJet2TruePz   =lightParton2.Pz()    ;
       lightJet2TrueE    =lightParton2.E()     ;		
       lightJet2TruePt   =lightParton2.Pt()    ;



      //Smeared values
       top1SmearedPt        =recoTopQuark.Pt()       ;
       top1SmearedEta       =recoTopQuark.Eta()      ;
       top1SmearedPhi       =recoTopQuark.Phi()      ;
       top1SmearedE         =recoTopQuark.E()        ;
       top1SmearedMass      =sqrt(recoTopQuark.M2()) ;
       W1SmearedPt          =recoWPlus.Pt()          ;
       W1SmearedEta         =recoWPlus.Eta()         ;
       W1SmearedPhi         =recoWPlus.Phi()         ;
       W1SmearedE           =recoWPlus.E()           ;
       W1SmearedMass        =sqrt(recoWPlus.M2())    ;
       bJet1SmearedPt       =smearedBottomQuark.Pt() ;
       bJet1SmearedEta      =smearedBottomQuark.Eta();
       bJet1SmearedPhi      =smearedBottomQuark.Phi();
       bJet1SmearedE        =smearedBottomQuark.E()  ;
      if (leptonFlag == 0){
       W1Daughter1SmearedPt =smearedAntiLepton.Pt()         ;
       W1Daughter1SmearedEta=smearedAntiLepton.Eta()        ;
       W1Daughter1SmearedPhi=smearedAntiLepton.Phi()        ;
       W1Daughter1SmearedE  =smearedAntiLepton.E()          ;
       W1Daughter2SmearedPt =-1.                     ;
       W1Daughter2SmearedEta=-1.                     ;
       W1Daughter2SmearedPhi=-1.                     ;
       W1Daughter2SmearedE  =-1.                     ;
      } else if (leptonFlag == 1){
       W1Daughter1SmearedPt =smearedLightQuark1.Pt()         ;
       W1Daughter1SmearedEta=smearedLightQuark1.Eta()        ;
       W1Daughter1SmearedPhi=smearedLightQuark1.Phi()        ;
       W1Daughter1SmearedE  =smearedLightQuark1.E()          ;
       W1Daughter2SmearedPt =smearedLightQuark2.Pt()                     ;
       W1Daughter2SmearedEta=smearedLightQuark2.Eta()                     ;
       W1Daughter2SmearedPhi=smearedLightQuark2.Phi()                     ;
       W1Daughter2SmearedE  =smearedLightQuark2.E()                    ;
      }
       lightJet1SmearedPx   =smearedLightParton1.Px();
       lightJet1SmearedPy   =smearedLightParton1.Py();
       lightJet1SmearedPz   =smearedLightParton1.Pz();
       lightJet1SmearedE    =smearedLightParton1.E() ;
       lightJet1SmearedPt   =smearedLightParton1.Pt();

       top2SmearedPt        =recoAntiTopQuark.Pt()       ;
       top2SmearedEta       =recoAntiTopQuark.Eta()      ;
       top2SmearedPhi       =recoAntiTopQuark.Phi()      ;
       top2SmearedE         =recoAntiTopQuark.E()        ;
       top2SmearedMass      =sqrt(recoAntiTopQuark.M2()) ;
       W2SmearedPt          =recoWMinus.Pt()             ;
       W2SmearedEta         =recoWMinus.Eta()            ;
       W2SmearedPhi         =recoWMinus.Phi()            ;
       W2SmearedE           =recoWMinus.E()              ;
       W2SmearedMass        =sqrt(recoWMinus.M2())       ;
       bJet2SmearedPt       =smearedAntiBottomQuark.Pt() ;
       bJet2SmearedEta      =smearedAntiBottomQuark.Eta();
       bJet2SmearedPhi      =smearedAntiBottomQuark.Phi();
       bJet2SmearedE        =smearedAntiBottomQuark.E()  ;
      if (leptonFlag == 0){
       W2Daughter1SmearedPt =smearedLightQuark1.Pt()     ;
       W2Daughter1SmearedEta=smearedLightQuark1.Eta()    ;
       W2Daughter1SmearedPhi=smearedLightQuark1.Phi()    ;
       W2Daughter1SmearedE  =smearedLightQuark1.E()      ;
       W2Daughter2SmearedPt =smearedLightQuark2.Pt()     ;
       W2Daughter2SmearedEta=smearedLightQuark2.Eta()    ;
       W2Daughter2SmearedPhi=smearedLightQuark2.Phi()    ;
       W2Daughter2SmearedE  =smearedLightQuark2.E()      ;
      } else if (leptonFlag == 1){
       W2Daughter1SmearedPt =smearedLepton.Pt()         ;
       W2Daughter1SmearedEta=smearedLepton.Eta()        ;
       W2Daughter1SmearedPhi=smearedLepton.Phi()        ;
       W2Daughter1SmearedE  =smearedLepton.E()          ;
       W2Daughter2SmearedPt =-1.                     ;
       W2Daughter2SmearedEta=-1.                     ;
       W2Daughter2SmearedPhi=-1.                     ;
       W2Daughter2SmearedE  =-1.                     ;
      }
       lightJet2SmearedPx   =smearedLightParton2.Px()    ;
       lightJet2SmearedPy   =smearedLightParton2.Py()    ;
       lightJet2SmearedPz   =smearedLightParton2.Pz()    ;
       lightJet2SmearedE    =smearedLightParton2.E()     ;
       lightJet2SmearedPt   =smearedLightParton2.Pt()    ;


      std::cout<<"before get best values"<<std::endl;
      //Best values


      //Minimizer deltas

      vector<double> bJetPtDeltas; 
      vector<double> bJetPhiDeltas;
      vector<double> bJetEtaDeltas;
      vector<double> firstWDaughterPtDeltas;
      vector<double> firstWDaughterPhiDeltas;
      vector<double> firstWDaughterEtaDeltas;
      vector<double> secondWDaughterPtDeltas;
      vector<double> secondWDaughterPhiDeltas;
      vector<double> secondWDaughterEtaDeltas;
      vector<double> WMassDeltas;
      vector<double> topMassDeltas;
      vector<double> lightJetPxDeltas;
      vector<double> lightJetPyDeltas;

      ev.getBestDeltas(bJetPtDeltas, bJetPhiDeltas, bJetEtaDeltas,
		       firstWDaughterPtDeltas, firstWDaughterPhiDeltas, firstWDaughterEtaDeltas,
		       secondWDaughterPtDeltas, secondWDaughterPhiDeltas, secondWDaughterEtaDeltas,
		       WMassDeltas, topMassDeltas,
		       lightJetPxDeltas, lightJetPyDeltas);


       bJet1PtDelta =bJetPtDeltas .at(0); 
       bJet1PhiDelta=bJetPhiDeltas.at(0); 
       bJet1EtaDelta=bJetEtaDeltas.at(0);

       bJet2PtDelta =bJetPtDeltas .at(1);
       bJet2PhiDelta=bJetPhiDeltas.at(1);
       bJet2EtaDelta=bJetEtaDeltas.at(1);

       W1Daughter1PtDelta =firstWDaughterPtDeltas .at(0); 
       W1Daughter1PhiDelta=firstWDaughterPhiDeltas.at(0);
       W1Daughter1EtaDelta=firstWDaughterEtaDeltas.at(0);

       W2Daughter1PtDelta =firstWDaughterPtDeltas .at(1);
       W2Daughter1PhiDelta=firstWDaughterPhiDeltas.at(1);
       W2Daughter1EtaDelta=firstWDaughterEtaDeltas.at(1);

       W1Daughter2PtDelta =secondWDaughterPtDeltas .at(0);
       W1Daughter2PhiDelta=secondWDaughterPhiDeltas.at(0);
       W1Daughter2EtaDelta=secondWDaughterEtaDeltas.at(0);

       W2Daughter2PtDelta =secondWDaughterPtDeltas .at(1);
       W2Daughter2PhiDelta=secondWDaughterPhiDeltas.at(1);
       W2Daughter2EtaDelta=secondWDaughterEtaDeltas.at(1);

       W1MassDelta=WMassDeltas.at(0);
       W2MassDelta=WMassDeltas.at(1);

       top1MassDelta=topMassDeltas.at(0);
       top2MassDelta=topMassDeltas.at(1);

       lightJet1PxDelta=lightJetPxDeltas.at(0);
       lightJet2PxDelta=lightJetPxDeltas.at(1);

       lightJet1PyDelta=lightJetPyDeltas.at(0);
       lightJet2PyDelta=lightJetPyDeltas.at(1);


      //Object momenta

      ev.getTop         (0,top1BestPx       ,top1BestPy       ,top1BestPz       ,top1BestE       );
      ev.getBJet        (0,bJet1BestPx      ,bJet1BestPy      ,bJet1BestPz      ,bJet1BestE      );
      ev.getW           (0,W1BestPx         ,W1BestPy         ,W1BestPz         ,W1BestE         );
      ev.getWDaughter1  (0,W1Daughter1BestPx,W1Daughter1BestPy,W1Daughter1BestPz,W1Daughter1BestE);
      ev.getWDaughter2  (0,W1Daughter2BestPx,W1Daughter2BestPy,W1Daughter2BestPz,W1Daughter2BestE);
      ev.getNonTopObject(0,lightJet1BestPx,lightJet1BestPy);

      XYZTLorentzVector top1Best;
      top1Best.SetPxPyPzE(top1BestPx,top1BestPy,top1BestPz,top1BestE);
      XYZTLorentzVector bJet1Best;
      bJet1Best.SetPxPyPzE(bJet1BestPx,bJet1BestPy,bJet1BestPz,bJet1BestE);
      XYZTLorentzVector W1Best;
      W1Best.SetPxPyPzE(W1BestPx,W1BestPy,W1BestPz,W1BestE);
      XYZTLorentzVector W1Daughter1Best;
      W1Daughter1Best.SetPxPyPzE(W1Daughter1BestPx,W1Daughter1BestPy,W1Daughter1BestPz,W1Daughter1BestE);
      XYZTLorentzVector W1Daughter2Best;
      W1Daughter2Best.SetPxPyPzE(W1Daughter2BestPx,W1Daughter2BestPy,W1Daughter2BestPz,W1Daughter2BestE);
      XYZTLorentzVector lightJet1Best;
      lightJet1BestPz=smearedLightParton1.Pz();
      lightJet1BestE =sqrt(smearedLightParton1.M2()+pow(lightJet1BestPx,2)+pow(lightJet1BestPy,2)+pow(smearedLightParton1.Pz(),2));
      lightJet1Best.SetPxPyPzE(lightJet1BestPx,lightJet1BestPy,lightJet1BestPz,lightJet1BestE);

       top1BestPt        =top1Best.Pt() ;
       top1BestEta       =top1Best.Eta();
       top1BestPhi       =top1Best.Phi();
       top1BestMass      =max(0.,sqrt(top1Best.M2()));
       bJet1BestPt       =bJet1Best.Pt() ;
       bJet1BestEta      =bJet1Best.Eta();
       bJet1BestPhi      =bJet1Best.Phi();
       W1Daughter1BestPt =W1Daughter1Best.Pt() ;
       W1Daughter1BestEta=W1Daughter1Best.Eta();
       W1Daughter1BestPhi=W1Daughter1Best.Phi();
       W1Daughter2BestPt =W1Daughter2Best.Pt() ;
       W1Daughter2BestEta=W1Daughter2Best.Eta();
       W1Daughter2BestPhi=W1Daughter2Best.Phi();
       W1BestPt          =W1Best.Pt() ;
       W1BestEta         =W1Best.Eta();
       W1BestPhi         =W1Best.Phi();
       W1BestMass        =max(0.,sqrt(W1Best.M2()));
       lightJet1BestPt   =lightJet1Best.Pt();


      ev.getTop         (1,top2BestPx       ,top2BestPy       ,top2BestPz       ,top2BestE       );
      ev.getBJet        (1,bJet2BestPx      ,bJet2BestPy      ,bJet2BestPz      ,bJet2BestE      );
      ev.getWDaughter1  (1,W2Daughter1BestPx,W2Daughter1BestPy,W2Daughter1BestPz,W2Daughter1BestE);
      ev.getWDaughter2  (1,W2Daughter2BestPx,W2Daughter2BestPy,W2Daughter2BestPz,W2Daughter2BestE);
      ev.getW           (1,W2BestPx         ,W2BestPy         ,W2BestPz         ,W2BestE         );
      ev.getNonTopObject(1,lightJet2BestPx,lightJet2BestPy);


      XYZTLorentzVector top2Best;
      top2Best.SetPxPyPzE(top2BestPx,top2BestPy,top2BestPz,top2BestE);
      XYZTLorentzVector bJet2Best;
      bJet2Best.SetPxPyPzE(bJet2BestPx,bJet2BestPy,bJet2BestPz,bJet2BestE);
      XYZTLorentzVector W2Best;
      W2Best.SetPxPyPzE(W2BestPx,W2BestPy,W2BestPz,W2BestE);
      XYZTLorentzVector W2Daughter1Best;
      W2Daughter1Best.SetPxPyPzE(W2Daughter1BestPx,W2Daughter1BestPy,W2Daughter1BestPz,W2Daughter1BestE);
      XYZTLorentzVector W2Daughter2Best;
      W2Daughter2Best.SetPxPyPzE(W2Daughter2BestPx,W2Daughter2BestPy,W2Daughter2BestPz,W2Daughter2BestE);
      XYZTLorentzVector lightJet2Best;
      lightJet2BestPz=smearedLightParton2.Pz();
      lightJet2BestE =sqrt(smearedLightParton2.M2()+pow(lightJet2BestPx,2)+pow(lightJet2BestPy,2)+pow(smearedLightParton2.Pz(),2));
      lightJet2Best.SetPxPyPzE(lightJet2BestPx,lightJet2BestPy,lightJet2BestPz,lightJet2BestE);

       top2BestPt        =top2Best.Pt() ;
       top2BestEta       =top2Best.Eta();
       top2BestPhi       =top2Best.Phi();
       top2BestMass      =max(0.,sqrt(top2Best.M2()));
       bJet2BestPt       =bJet2Best.Pt() ;
       bJet2BestEta      =bJet2Best.Eta();
       bJet2BestPhi      =bJet2Best.Phi();
       W2Daughter1BestPt =W2Daughter1Best.Pt() ;
       W2Daughter1BestEta=W2Daughter1Best.Eta();
       W2Daughter1BestPhi=W2Daughter1Best.Phi();
       W2Daughter2BestPt =W2Daughter2Best.Pt() ;
       W2Daughter2BestEta=W2Daughter2Best.Eta();
       W2Daughter2BestPhi=W2Daughter2Best.Phi();
       W2BestPt          =W2Best.Pt() ;
       W2BestEta         =W2Best.Eta();
       W2BestPhi         =W2Best.Phi();
       W2BestMass        =max(0.,sqrt(W2Best.M2()));
      lightJet2BestPt   =lightJet2Best.Pt();

      std::cout<<"before compare"<<std::endl;

      //Compare true and smeared
       top1DeltaPtTrueSmeared       =top1TruePt-top1SmearedPt; 
       top1DeltaRTrueSmeared        =deltaR(topQuark,recoTopQuark); 
       top1DeltaMTrueSmeared        =top1TrueMass-top1SmearedMass;
       W1DeltaPtTrueSmeared         =W1TruePt-W1SmearedPt; 
       W1DeltaRTrueSmeared          =deltaR(WPlus,recoWPlus); 
       W1DeltaMTrueSmeared          =W1TrueMass-W1SmearedMass;
       bJet1DeltaPtTrueSmeared      =bJet1TruePt-bJet1SmearedPt; 
       bJet1DeltaRTrueSmeared       =deltaR(bottomQuark,smearedBottomQuark);
      if (leptonFlag == 0){
       W1Daughter1DeltaPtTrueSmeared=W1Daughter1TruePt-W1Daughter1SmearedPt;
       W1Daughter1DeltaRTrueSmeared =deltaR(antiLepton,smearedAntiLepton);
       W1Daughter2DeltaPtTrueSmeared=W1Daughter2TruePt-met.pt();
       W1Daughter2DeltaRTrueSmeared =deltaR(neutrino,met);
      } else if (leptonFlag == 1){
       W1Daughter1DeltaPtTrueSmeared=W1Daughter1TruePt-W1Daughter1SmearedPt;
       W1Daughter1DeltaRTrueSmeared =deltaR(qFromW,smearedLightQuark1);
       W1Daughter2DeltaPtTrueSmeared=W1Daughter2TruePt-W1Daughter2SmearedPt;
       W1Daughter2DeltaRTrueSmeared =deltaR(qbarFromW,smearedLightQuark2);
      }
       lightJet1DeltaPxTrueSmeared  =lightJet1TruePx-lightJet1SmearedPx;
       lightJet1DeltaPyTrueSmeared  =lightJet1TruePy-lightJet1SmearedPy;
       lightJet1DeltaRTrueSmeared   =deltaR(lightParton1,smearedLightParton1);

       top2DeltaPtTrueSmeared       =top2TruePt-top2SmearedPt;
       top2DeltaRTrueSmeared        =deltaR(antiTopQuark,recoAntiTopQuark);
       top2DeltaMTrueSmeared        =top2TrueMass-top2SmearedMass;
       W2DeltaPtTrueSmeared         =W2TruePt-W2SmearedPt;
       W2DeltaRTrueSmeared          =deltaR(WMinus,recoWMinus);
       W2DeltaMTrueSmeared          =W2TrueMass-W2SmearedMass;
       bJet2DeltaPtTrueSmeared      =bJet2TruePt-bJet2SmearedPt;
       bJet2DeltaRTrueSmeared       =deltaR(antiBottomQuark,smearedAntiBottomQuark);
      if (leptonFlag == 0){
       W2Daughter1DeltaPtTrueSmeared=W2Daughter1TruePt-W2Daughter1SmearedPt;
       W2Daughter1DeltaRTrueSmeared =deltaR(qFromW,smearedLightQuark1);
       W2Daughter2DeltaPtTrueSmeared=W2Daughter2TruePt-W2Daughter2SmearedPt;
       W2Daughter2DeltaRTrueSmeared =deltaR(qbarFromW,smearedLightQuark2);
      } else if (leptonFlag == 1){
       W2Daughter1DeltaPtTrueSmeared=W2Daughter1TruePt-W2Daughter1SmearedPt;
       W2Daughter1DeltaRTrueSmeared =deltaR(lepton,smearedLepton);
       W2Daughter2DeltaPtTrueSmeared=W2Daughter2TruePt-met.pt();
       W2Daughter2DeltaRTrueSmeared =deltaR(antiNeutrino,met);
      }
       lightJet2DeltaPxTrueSmeared  =lightJet2TruePx-lightJet2SmearedPx;
       lightJet2DeltaPyTrueSmeared  =lightJet2TruePy-lightJet2SmearedPy;
       lightJet2DeltaRTrueSmeared   =deltaR(lightParton2,smearedLightParton2);



      //Compare true and corrected
       top1DeltaPtTrueBest       =top1TruePt-top1BestPt;
       top1DeltaRTrueBest        =deltaR(topQuark,top1Best);
       top1DeltaMTrueBest        =top1TrueMass-top1BestMass;
       W1DeltaPtTrueBest         =W1TruePt-W1BestPt;
       W1DeltaRTrueBest          =deltaR(WPlus,W1Best);
       W1DeltaMTrueBest          =W1TrueMass-W1BestMass;
       bJet1DeltaPtTrueBest      =bJet1TruePt-bJet1BestPt;
       bJet1DeltaRTrueBest       =deltaR(bottomQuark,bJet1Best);
      if (leptonFlag == 0){
       W1Daughter1DeltaPtTrueBest=W1Daughter1TruePt-W1Daughter1BestPt;
       W1Daughter1DeltaRTrueBest =deltaR(antiLepton,W1Daughter1Best);
       W1Daughter2DeltaPtTrueBest=W1Daughter2TruePt-W1Daughter2BestPt;
       W1Daughter2DeltaRTrueBest =deltaR(neutrino,W1Daughter2Best);
      } else if (leptonFlag == 1){
       W1Daughter1DeltaPtTrueBest=W1Daughter1TruePt-W1Daughter1BestPt;
       W1Daughter1DeltaRTrueBest =deltaR(qFromW,W1Daughter1Best);
       W1Daughter2DeltaPtTrueBest=W1Daughter2TruePt-W1Daughter2BestPt;
       W1Daughter2DeltaRTrueBest =deltaR(qbarFromW,W1Daughter2Best);
      }
       lightJet1DeltaPxTrueBest  =lightJet1TruePx-lightJet1BestPx;
       lightJet1DeltaPyTrueBest  =lightJet1TruePy-lightJet1BestPy;
       lightJet1DeltaRTrueBest   =deltaR(lightParton1,lightJet1Best);

       top2DeltaPtTrueBest       =top2TruePt-top2BestPt;
       top2DeltaRTrueBest        =deltaR(antiTopQuark,top2Best);
       top2DeltaMTrueBest        =top2TrueMass-top2BestMass;
       W2DeltaPtTrueBest         =W2TruePt-W2BestPt;
       W2DeltaRTrueBest          =deltaR(WMinus,W2Best);
       W2DeltaMTrueBest          =W2TrueMass-W2BestMass;
       bJet2DeltaPtTrueBest      =bJet2TruePt-bJet2BestPt;
       bJet2DeltaRTrueBest       =deltaR(antiBottomQuark,bJet2Best);
      if (leptonFlag == 0){
       W2Daughter1DeltaPtTrueBest=W2Daughter1TruePt-W2Daughter1BestPt;
       W2Daughter1DeltaRTrueBest =deltaR(qFromW,W2Daughter1Best);
       W2Daughter2DeltaPtTrueBest=W2Daughter2TruePt-W2Daughter2BestPt;
       W2Daughter2DeltaRTrueBest =deltaR(qbarFromW,W2Daughter2Best);
      } else if (leptonFlag == 1){
       W2Daughter1DeltaPtTrueBest=W2Daughter1TruePt-W2Daughter1BestPt;
       W2Daughter1DeltaRTrueBest =deltaR(antiNeutrino,W2Daughter1Best);
       W2Daughter2DeltaPtTrueBest=W2Daughter2TruePt-W2Daughter2BestPt;
       W2Daughter2DeltaRTrueBest =deltaR(lepton,W2Daughter2Best);
      }
       lightJet2DeltaPxTrueBest  =lightJet2TruePx-lightJet2BestPx;
       lightJet2DeltaPyTrueBest  =lightJet2TruePy-lightJet2BestPy;
       lightJet2DeltaRTrueBest   =deltaR(lightParton2,lightJet2Best);


      //Compare smeared and corrected
       top1DeltaPtSmearedBest       =top1SmearedPt-top1BestPt;
       top1DeltaRSmearedBest        =deltaR(recoTopQuark,top1Best);
       top1DeltaMSmearedBest        =top1SmearedMass-top1BestMass;
       W1DeltaPtSmearedBest         =W1SmearedPt-W1BestPt;
       W1DeltaRSmearedBest          =deltaR(recoWPlus, W1Best);
       W1DeltaMSmearedBest          =W1SmearedMass-W1BestMass;
       bJet1DeltaPtSmearedBest      =bJet1SmearedPt-bJet1BestPt;
       bJet1DeltaRSmearedBest       =deltaR(smearedBottomQuark,bJet1Best);
      if (leptonFlag == 0){
       W1Daughter1DeltaPtSmearedBest=W1Daughter1SmearedPt-W1Daughter1BestPt;
       W1Daughter1DeltaRSmearedBest =deltaR(smearedAntiLepton,W1Daughter1Best);
       W1Daughter2DeltaPtSmearedBest=met.pt()-W1Daughter2BestPt;
       W1Daughter2DeltaRSmearedBest =deltaR(met,W1Daughter2Best);
      } else if (leptonFlag == 1){
       W1Daughter1DeltaPtSmearedBest=W1Daughter1SmearedPt-W1Daughter1BestPt;
       W1Daughter1DeltaRSmearedBest =deltaR(smearedLightQuark1,W1Daughter1Best);
       W1Daughter2DeltaPtSmearedBest=W1Daughter2SmearedPt-W1Daughter2BestPt;
       W1Daughter2DeltaRSmearedBest =deltaR(smearedLightQuark2,W1Daughter2Best);
      }
       lightJet1DeltaPxSmearedBest  =lightJet1SmearedPx-lightJet1BestPx;
       lightJet1DeltaPySmearedBest  =lightJet1SmearedPy-lightJet1BestPy;
       lightJet1DeltaRSmearedBest   =deltaR(smearedLightParton1,lightJet1Best);

       top2DeltaPtSmearedBest       =top2SmearedPt-top2BestPt;
       top2DeltaRSmearedBest        =deltaR(recoAntiTopQuark,top2Best);
       top2DeltaMSmearedBest        =top2SmearedMass-top2BestMass;
       W2DeltaPtSmearedBest         =W2SmearedPt-W2BestPt;
       W2DeltaRSmearedBest          =deltaR(recoWMinus,W2Best);
       W2DeltaMSmearedBest          =W2SmearedMass-W2BestMass;
       bJet2DeltaPtSmearedBest      =bJet2SmearedPt-bJet2BestPt;
       bJet2DeltaRSmearedBest       =deltaR(smearedAntiBottomQuark,bJet2Best);
      if (leptonFlag == 0){
       W2Daughter1DeltaPtSmearedBest=W2Daughter1SmearedPt-W2Daughter1BestPt;
       W2Daughter1DeltaRSmearedBest =deltaR(smearedLightQuark1,W2Daughter1Best);
       W2Daughter2DeltaPtSmearedBest=W2Daughter2SmearedPt-W2Daughter2BestPt;
       W2Daughter2DeltaRSmearedBest =deltaR(smearedLightQuark2,W2Daughter2Best);
      } else if (leptonFlag == 1){
       W2Daughter1DeltaPtSmearedBest=W2Daughter1SmearedPt-W2Daughter1BestPt;
       W2Daughter1DeltaRSmearedBest =deltaR(smearedLepton,W2Daughter1Best);
       W2Daughter2DeltaPtSmearedBest=met.pt()-W2Daughter2BestPt;
       W2Daughter2DeltaRSmearedBest =deltaR(met,W2Daughter2Best);
      }
       lightJet2DeltaPxSmearedBest  =lightJet2SmearedPx-lightJet2BestPx;
       lightJet2DeltaPySmearedBest  =lightJet2SmearedPy-lightJet2BestPy;
       lightJet2DeltaRSmearedBest   =deltaR(smearedLightParton2,lightJet2Best);

      
      outTree->Fill();
      std::cout<<"hello!"<<std::endl;

      if (leptonFlag == 0 and innerMinStatus == 0 and outerMinStatus == 0){
      //fill TS and TC hists
      leptonicBottomPtTS->Fill(bJet1DeltaPtTrueSmeared);
      leptonicBottomPtTC->Fill(bJet1DeltaPtTrueBest);
      leptonPtTS->Fill(W1Daughter1DeltaPtTrueSmeared);
      leptonPtTC->Fill(W1Daughter1DeltaPtTrueBest);
      neutrinoPtTS->Fill(W1Daughter2DeltaPtTrueSmeared);
      neutrinoPtTC->Fill(W1Daughter2DeltaPtTrueBest);
      leptonicTopPtTS->Fill(top1DeltaPtTrueSmeared);
      leptonicTopPtTC->Fill(top1DeltaPtTrueBest);
      leptonicWPtTS->Fill(W1DeltaPtTrueSmeared);
      leptonicWPtTC->Fill(W1DeltaPtTrueBest);

      hadronicBottomPtTS->Fill(bJet2DeltaPtTrueSmeared);
      hadronicBottomPtTC->Fill(bJet2DeltaPtTrueBest);
      hadronicQuarkPtTS->Fill(W2Daughter1DeltaPtTrueSmeared);
      hadronicQuarkPtTC->Fill(W2Daughter1DeltaPtTrueBest);
      hadronicAntiQuarkPtTS->Fill(W2Daughter2DeltaPtTrueSmeared);
      hadronicAntiQuarkPtTC->Fill(W2Daughter2DeltaPtTrueBest);
      hadronicTopPtTS->Fill(top2DeltaPtTrueSmeared);
      hadronicTopPtTC->Fill(top2DeltaPtTrueBest);
      hadronicWPtTS->Fill(W2DeltaPtTrueSmeared);
      hadronicWPtTC->Fill(W2DeltaPtTrueBest);

      lightBottomPxTS->Fill(lightJet1DeltaPxTrueSmeared);
      lightBottomPxTC->Fill(lightJet1DeltaPxTrueBest);
      lightAntiBottomPxTS->Fill(lightJet2DeltaPxTrueSmeared);
      lightAntiBottomPxTC->Fill(lightJet2DeltaPxTrueBest);

      std::cout<<"after fill TC TS" <<std::endl;
      //fill SC hist
      leptonicBottomPtSC->Fill(bJet1DeltaPtSmearedBest);
      std::cout<<"after fill first SC" <<std::endl;
      leptonPtSC->Fill(W1Daughter1DeltaPtSmearedBest);
      neutrinoPtSC->Fill(W1Daughter2DeltaPtSmearedBest);
      leptonicTopPtSC->Fill(top1DeltaPtSmearedBest);
      leptonicWPtSC->Fill(W1DeltaPtSmearedBest);

      hadronicBottomPtSC->Fill(bJet2DeltaPtSmearedBest);
      hadronicQuarkPtSC->Fill(W2Daughter1DeltaPtSmearedBest);
      hadronicAntiQuarkPtSC->Fill(W2Daughter2DeltaPtSmearedBest);
      hadronicTopPtSC->Fill(top2DeltaPtSmearedBest);
      hadronicWPtSC->Fill(W2DeltaPtSmearedBest);

      lightBottomPxSC->Fill(lightJet1DeltaPxSmearedBest);
      lightAntiBottomPxSC->Fill(lightJet2DeltaPxSmearedBest);

      std::cout<<"after fill last SC"<<std::endl;
      }

      if (leptonFlag == 1 and innerMinStatus == 0 and outerMinStatus == 0){
      //fill TS and TC hists
      leptonicBottomPtTS->Fill(bJet2DeltaPtTrueSmeared);
      leptonicBottomPtTC->Fill(bJet2DeltaPtTrueBest);
      leptonPtTS->Fill(W2Daughter1DeltaPtTrueSmeared);
      leptonPtTC->Fill(W2Daughter1DeltaPtTrueBest);
      neutrinoPtTS->Fill(W2Daughter2DeltaPtTrueSmeared);
      neutrinoPtTC->Fill(W2Daughter2DeltaPtTrueBest);
      leptonicTopPtTS->Fill(top2DeltaPtTrueSmeared);
      leptonicTopPtTC->Fill(top2DeltaPtTrueBest);
      leptonicWPtTS->Fill(W2DeltaPtTrueSmeared);
      leptonicWPtTC->Fill(W2DeltaPtTrueBest);

      hadronicBottomPtTS->Fill(bJet1DeltaPtTrueSmeared);
      hadronicBottomPtTC->Fill(bJet1DeltaPtTrueBest);
      hadronicQuarkPtTS->Fill(W1Daughter1DeltaPtTrueSmeared);
      hadronicQuarkPtTC->Fill(W1Daughter1DeltaPtTrueBest);
      hadronicAntiQuarkPtTS->Fill(W1Daughter2DeltaPtTrueSmeared);
      hadronicAntiQuarkPtTC->Fill(W1Daughter2DeltaPtTrueBest);
      hadronicTopPtTS->Fill(top1DeltaPtTrueSmeared);
      hadronicTopPtTC->Fill(top1DeltaPtTrueBest);
      hadronicWPtTS->Fill(W1DeltaPtTrueSmeared);
      hadronicWPtTC->Fill(W1DeltaPtTrueBest);

      lightBottomPxTS->Fill(lightJet1DeltaPxTrueSmeared);
      lightBottomPxTC->Fill(lightJet1DeltaPxTrueBest);
      lightAntiBottomPxTS->Fill(lightJet2DeltaPxTrueSmeared);
      lightAntiBottomPxTC->Fill(lightJet2DeltaPxTrueBest);

      std::cout<<"after fill TC TS" <<std::endl;
      //fill SC hist
      leptonicBottomPtSC->Fill(bJet2DeltaPtSmearedBest);
      std::cout<<"after fill first SC" <<std::endl;
      leptonPtSC->Fill(W2Daughter1DeltaPtSmearedBest);
      neutrinoPtSC->Fill(W2Daughter2DeltaPtSmearedBest);
      leptonicTopPtSC->Fill(top2DeltaPtSmearedBest);
      leptonicWPtSC->Fill(W2DeltaPtSmearedBest);

      hadronicBottomPtSC->Fill(bJet1DeltaPtSmearedBest);
      hadronicQuarkPtSC->Fill(W1Daughter1DeltaPtSmearedBest);
      hadronicAntiQuarkPtSC->Fill(W1Daughter2DeltaPtSmearedBest);
      hadronicTopPtSC->Fill(top1DeltaPtSmearedBest);
      hadronicWPtSC->Fill(W1DeltaPtSmearedBest);

      lightBottomPxSC->Fill(lightJet1DeltaPxSmearedBest);
      lightAntiBottomPxSC->Fill(lightJet2DeltaPxSmearedBest);

      std::cout<<"after fill last SC"<<std::endl;
      }

      //std::cout<<"bJet1SmearedE = "<< bJet1SmearedE<<std::endl;
      //std::cout<<"bJet1TrueE = "<<bJet1TrueE<<std::endl;
   } //end loop over events

   int mc1 = 5;
   int mc2 = 4;
   int mc3 = 3;
   int mc0 = 0;

   //compare TS and TC
   std::cout<<"drawing TS TC canvases"<<std::endl;
   leptonicBottomPtTS->SetFillColor(mc1);
   leptonicBottomPtTS->SetLineColor(mc1);
   leptonicBottomPtTC->SetLineColor(mc2);
   cleptonicBottomPt->cd();
   leptonicBottomPtTS->GetXaxis()->SetTitle("Pt res (GeV)");
   leptonicBottomPtTS->GetYaxis()->SetTitle("Events");
   leptonicBottomPtTS->SetTitle("Leptonic Bottom");
   leptonicBottomPtTS->SetMaximum( max( leptonicBottomPtTS->GetMaximum(), leptonicBottomPtTC->GetMaximum())+1 );
   leptonicBottomPtTS->Draw("HIST");
   leptonicBottomPtTC->Draw("SAMES");
   //cout<<"before func call"<<endl;
   //gPad->Update();
   moveStatsBox(leptonicBottomPtTC);
   //cout<<"after func call"<<endl;
   cleptonicBottomPt->Write();
   cleptonicBottomPt->ls();
   //compare SC
   std::cout<<"drawing SC canvases"<<std::endl;
   leptonicBottomPtSC->SetFillColor(mc3);
   leptonicBottomPtSC->SetLineColor(mc0);
   cleptonicBottomPtSC->cd();
   leptonicBottomPtSC->GetXaxis()->SetTitle("Pt_reco - Pt_corr (GeV)");
   leptonicBottomPtSC->GetYaxis()->SetTitle("Events");
   leptonicBottomPtSC->SetTitle("Leptonic Bottom");
   leptonicBottomPtSC->DrawCopy("HIST");
   cleptonicBottomPtSC->Write();
   cleptonicBottomPtSC->ls();

   //compare TS and TC
   std::cout<<"drawing TS TC canvases"<<std::endl;
   leptonicTopPtTS->SetFillColor(mc1);
   leptonicTopPtTS->SetLineColor(mc1);
   leptonicTopPtTC->SetLineColor(mc2);
   cleptonicTopPt->cd();
   leptonicTopPtTS->GetXaxis()->SetTitle("Pt res (GeV)");
   leptonicTopPtTS->GetYaxis()->SetTitle("Events");
   leptonicTopPtTS->SetTitle("Leptonic Top");
   leptonicTopPtTS->SetMaximum( max( leptonicTopPtTS->GetMaximum(), leptonicTopPtTC->GetMaximum())+1 );
   leptonicTopPtTS->Draw("HIST");
   leptonicTopPtTC->Draw("SAMES");
   moveStatsBox(leptonicTopPtTC);
   cleptonicTopPt->Write();
   cleptonicTopPt->ls();
   //compare SC
   std::cout<<"drawing SC canvases"<<std::endl;
   leptonicTopPtSC->SetFillColor(mc3);
   leptonicTopPtSC->SetLineColor(mc0);
   cleptonicTopPtSC->cd();
   leptonicTopPtSC->GetXaxis()->SetTitle("Pt_reco - Pt_corr (GeV)");
   leptonicTopPtSC->GetYaxis()->SetTitle("Events");
   leptonicTopPtSC->SetTitle("Leptonic Top");
   leptonicTopPtSC->DrawCopy("HIST");
   cleptonicTopPtSC->Write();
   cleptonicTopPtSC->ls();

   //compare TS and TC
   std::cout<<"drawing TS TC canvases"<<std::endl;
   leptonicWPtTS->SetFillColor(mc1);
   leptonicWPtTS->SetLineColor(mc1);
   leptonicWPtTC->SetLineColor(mc2);
   cleptonicWPt->cd();
   leptonicWPtTS->GetXaxis()->SetTitle("Pt res (GeV)");
   leptonicWPtTS->GetYaxis()->SetTitle("Events");
   leptonicWPtTS->SetTitle("Leptonic W");
   leptonicWPtTS->SetMaximum( max( leptonicWPtTS->GetMaximum(), leptonicWPtTC->GetMaximum())+1 );
   leptonicWPtTS->Draw("HIST");
   leptonicWPtTC->Draw("SAMES");   
   moveStatsBox(leptonicWPtTC);
   cleptonicWPt->Write();
   cleptonicWPt->ls();
   //compare SC
   std::cout<<"drawing SC canvases"<<std::endl;
   leptonicWPtSC->SetFillColor(mc3);
   leptonicWPtSC->SetLineColor(mc0);
   cleptonicWPtSC->cd();
   leptonicWPtSC->GetXaxis()->SetTitle("Pt_reco - Pt_corr (GeV)");
   leptonicWPtSC->GetYaxis()->SetTitle("Events");
   leptonicWPtSC->SetTitle("Leptonic W");
   leptonicWPtSC->DrawCopy("HIST");
   cleptonicWPtSC->Write();
   cleptonicWPtSC->ls();

   //compare TS and TC
   std::cout<<"drawing TS TC canvases"<<std::endl;
   leptonPtTS->SetFillColor(mc1);
   leptonPtTS->SetLineColor(mc1);
   leptonPtTC->SetLineColor(mc2);
   cleptonPt->cd();
   leptonPtTS->GetXaxis()->SetTitle("Pt res (GeV)");
   leptonPtTS->GetYaxis()->SetTitle("Events");
   leptonPtTS->SetTitle("Lepton");
   leptonPtTS->SetMaximum( max( leptonPtTS->GetMaximum(), leptonPtTC->GetMaximum())+1 );
   leptonPtTS->Draw("HIST");
   leptonPtTC->Draw("SAMES");   
   moveStatsBox(leptonPtTC);
   cleptonPt->Write();
   cleptonPt->ls();
   //compare SC
   std::cout<<"drawing SC canvases"<<std::endl;
   leptonPtSC->SetFillColor(mc3);
   leptonPtSC->SetLineColor(mc0);
   cleptonPtSC->cd();
   leptonPtSC->GetXaxis()->SetTitle("Pt_reco - Pt_corr (GeV)");
   leptonPtSC->GetYaxis()->SetTitle("Events");
   leptonPtSC->SetTitle("Lepton");
   leptonPtSC->DrawCopy("HIST");
   cleptonPtSC->Write();
   cleptonPtSC->ls();

   //compare TS and TC
   std::cout<<"drawing TS TC canvases"<<std::endl;
   neutrinoPtTS->SetFillColor(mc1);
   neutrinoPtTS->SetLineColor(mc1);
   neutrinoPtTC->SetLineColor(mc2);
   cneutrinoPt->cd();
   neutrinoPtTS->GetXaxis()->SetTitle("Pt res (GeV)");
   neutrinoPtTS->GetYaxis()->SetTitle("Events");
   neutrinoPtTS->SetTitle("Neutrino");
   neutrinoPtTS->SetMaximum( max( neutrinoPtTS->GetMaximum(), neutrinoPtTC->GetMaximum())+1 );
   neutrinoPtTS->Draw("HIST");
   neutrinoPtTC->Draw("SAMES");   
   moveStatsBox(neutrinoPtTC);
   cneutrinoPt->Write();
   cneutrinoPt->ls();
   //compare SC
   std::cout<<"drawing SC canvases"<<std::endl;
   neutrinoPtSC->SetFillColor(mc3);
   neutrinoPtSC->SetLineColor(mc0);
   cneutrinoPtSC->cd();
   neutrinoPtSC->GetXaxis()->SetTitle("Pt_met - Pt_corr (GeV)");
   neutrinoPtSC->GetYaxis()->SetTitle("Events");
   neutrinoPtSC->SetTitle("Neutrino");
   neutrinoPtSC->DrawCopy("HIST");
   cneutrinoPtSC->Write();
   cneutrinoPtSC->ls();

   //compare TS and TC
   std::cout<<"drawing TS TC canvases"<<std::endl;
   hadronicBottomPtTS->SetFillColor(mc1);
   hadronicBottomPtTS->SetLineColor(mc1);
   hadronicBottomPtTC->SetLineColor(mc2);
   chadronicBottomPt->cd();
   hadronicBottomPtTS->GetXaxis()->SetTitle("Pt res (GeV)");
   hadronicBottomPtTS->GetYaxis()->SetTitle("Events");
   hadronicBottomPtTS->SetTitle("Hadronic Bottom");
   hadronicBottomPtTS->SetMaximum( max( hadronicBottomPtTS->GetMaximum(), hadronicBottomPtTC->GetMaximum())+1 );
   hadronicBottomPtTS->Draw("HIST");
   hadronicBottomPtTC->Draw("SAMES");   
   moveStatsBox(hadronicBottomPtTC);
   chadronicBottomPt->Write();
   chadronicBottomPt->ls();
   //compare SC
   std::cout<<"drawing SC canvases"<<std::endl;
   hadronicBottomPtSC->SetFillColor(mc3);
   hadronicBottomPtSC->SetLineColor(mc0);
   chadronicBottomPtSC->cd();
   hadronicBottomPtSC->GetXaxis()->SetTitle("Pt_reco - Pt_corr (GeV)");
   hadronicBottomPtSC->GetYaxis()->SetTitle("Events");
   hadronicBottomPtSC->SetTitle("Hadronic Bottom");
   hadronicBottomPtSC->DrawCopy("HIST");
   chadronicBottomPtSC->Write();
   chadronicBottomPtSC->ls();

   //compare TS and TC
   std::cout<<"drawing TS TC canvases"<<std::endl;
   hadronicTopPtTS->SetFillColor(mc1);
   hadronicTopPtTS->SetLineColor(mc1);
   hadronicTopPtTC->SetLineColor(mc2);
   chadronicTopPt->cd();
   hadronicTopPtTS->GetXaxis()->SetTitle("Pt res (GeV)");
   hadronicTopPtTS->GetYaxis()->SetTitle("Events");
   hadronicTopPtTS->SetTitle("Hadronic Top");
   hadronicTopPtTS->SetMaximum( max( hadronicTopPtTS->GetMaximum(), hadronicTopPtTC->GetMaximum())+1 );
   hadronicTopPtTS->Draw("HIST");
   hadronicTopPtTC->Draw("SAMES");   
   moveStatsBox(hadronicTopPtTC);
   chadronicTopPt->Write();
   chadronicTopPt->ls();
   //compare SC
   std::cout<<"drawing SC canvases"<<std::endl;
   hadronicTopPtSC->SetFillColor(mc3);
   hadronicTopPtSC->SetLineColor(mc0);
   chadronicTopPtSC->cd();
   hadronicTopPtSC->GetXaxis()->SetTitle("Pt_reco - Pt_corr (GeV)");
   hadronicTopPtSC->GetYaxis()->SetTitle("Events");
   hadronicTopPtSC->SetTitle("Hadronic Top");
   hadronicTopPtSC->DrawCopy("HIST");
   chadronicTopPtSC->Write();
   chadronicTopPtSC->ls();

   //compare TS and TC
   std::cout<<"drawing TS TC canvases"<<std::endl;
   hadronicWPtTS->SetFillColor(mc1);
   hadronicWPtTS->SetLineColor(mc1);
   hadronicWPtTC->SetLineColor(mc2);
   chadronicWPt->cd();
   hadronicWPtTS->GetXaxis()->SetTitle("Pt res (GeV)");
   hadronicWPtTS->GetYaxis()->SetTitle("Events");
   hadronicWPtTS->SetTitle("Hadronic W");
   hadronicWPtTS->SetMaximum( max( hadronicWPtTS->GetMaximum(), hadronicWPtTC->GetMaximum())+1 );
   hadronicWPtTS->Draw("HIST");
   hadronicWPtTC->Draw("SAMES");   
   moveStatsBox(hadronicWPtTC);
   chadronicWPt->Write();
   chadronicWPt->ls();
   //compare SC
   std::cout<<"drawing SC canvases"<<std::endl;
   hadronicWPtSC->SetFillColor(mc3);
   hadronicWPtSC->SetLineColor(mc0);
   chadronicWPtSC->cd();
   hadronicWPtSC->GetXaxis()->SetTitle("Pt_reco - Pt_corr (GeV)");
   hadronicWPtSC->GetYaxis()->SetTitle("Events");
   hadronicWPtSC->SetTitle("Hadronic W");
   hadronicWPtSC->DrawCopy("HIST");
   chadronicWPtSC->Write();
   chadronicWPtSC->ls();

   //compare TS and TC
   std::cout<<"drawing TS TC canvases"<<std::endl;
   hadronicQuarkPtTS->SetFillColor(mc1);
   hadronicQuarkPtTS->SetLineColor(mc1);
   hadronicQuarkPtTC->SetLineColor(mc2);
   chadronicQuarkPt->cd();
   hadronicQuarkPtTS->GetXaxis()->SetTitle("Pt res (GeV)");
   hadronicQuarkPtTS->GetYaxis()->SetTitle("Events");
   hadronicQuarkPtTS->SetTitle("Quark from W");
   hadronicQuarkPtTS->SetMaximum( max( hadronicQuarkPtTS->GetMaximum(), hadronicQuarkPtTC->GetMaximum())+1 );
   hadronicQuarkPtTS->Draw("HIST");
   hadronicQuarkPtTC->Draw("SAMES");   
   moveStatsBox(hadronicQuarkPtTC);
   chadronicQuarkPt->Write();
   chadronicQuarkPt->ls();
   //compare SC
   std::cout<<"drawing SC canvases"<<std::endl;
   hadronicQuarkPtSC->SetFillColor(mc3);
   hadronicQuarkPtSC->SetLineColor(mc0);
   chadronicQuarkPtSC->cd();
   hadronicQuarkPtSC->GetXaxis()->SetTitle("Pt_reco - Pt_corr (GeV)");
   hadronicQuarkPtSC->GetYaxis()->SetTitle("Events");
   hadronicQuarkPtSC->SetTitle("Quark from W");
   hadronicQuarkPtSC->DrawCopy("HIST");
   chadronicQuarkPtSC->Write();
   chadronicQuarkPtSC->ls();

   //compare TS and TC
   std::cout<<"drawing TS TC canvases"<<std::endl;
   hadronicAntiQuarkPtTS->SetFillColor(mc1);
   hadronicAntiQuarkPtTS->SetLineColor(mc1);
   hadronicAntiQuarkPtTC->SetLineColor(mc2);
   chadronicAntiQuarkPt->cd();
   hadronicAntiQuarkPtTS->GetXaxis()->SetTitle("Pt res (GeV)");
   hadronicAntiQuarkPtTS->GetYaxis()->SetTitle("Events");
   hadronicAntiQuarkPtTS->SetTitle("Antiquark from W");
   hadronicAntiQuarkPtTS->SetMaximum( max( hadronicAntiQuarkPtTS->GetMaximum(), hadronicAntiQuarkPtTC->GetMaximum())+1 );
   hadronicAntiQuarkPtTS->Draw("HIST");
   hadronicAntiQuarkPtTC->Draw("SAMES");   
   moveStatsBox(hadronicAntiQuarkPtTC);
   chadronicAntiQuarkPt->Write();
   chadronicAntiQuarkPt->ls();
   //compare SC
   std::cout<<"drawing SC canvases"<<std::endl;
   hadronicAntiQuarkPtSC->SetFillColor(mc3);
   hadronicAntiQuarkPtSC->SetLineColor(mc0);
   chadronicAntiQuarkPtSC->cd();
   hadronicAntiQuarkPtSC->GetXaxis()->SetTitle("Pt_reco - Pt_corr (GeV)");
   hadronicAntiQuarkPtSC->GetYaxis()->SetTitle("Events");
   hadronicAntiQuarkPtSC->SetTitle("Antiquark from W");
   hadronicAntiQuarkPtSC->DrawCopy("HIST");
   chadronicAntiQuarkPtSC->Write();
   chadronicAntiQuarkPtSC->ls();

   //compare TS and TC
   std::cout<<"drawing TS TC canvases"<<std::endl;
   lightBottomPxTS->SetFillColor(mc1);
   lightBottomPxTS->SetLineColor(mc1);
   lightBottomPxTC->SetLineColor(mc2);
   clightBottomPx->cd();
   lightBottomPxTS->GetXaxis()->SetTitle("Px res (GeV)");
   lightBottomPxTS->GetYaxis()->SetTitle("Events");
   lightBottomPxTS->SetTitle("Bottom from H");
   lightBottomPxTS->SetMaximum( max( lightBottomPxTS->GetMaximum(), lightBottomPxTC->GetMaximum())+1 );
   lightBottomPxTS->Draw("HIST");
   lightBottomPxTC->Draw("SAMES");   
   moveStatsBox(lightBottomPxTC);
   clightBottomPx->Write();
   clightBottomPx->ls();
   //compare SC
   std::cout<<"drawing SC canvases"<<std::endl;
   lightBottomPxSC->SetFillColor(mc3);
   lightBottomPxSC->SetLineColor(mc0);
   clightBottomPxSC->cd();
   lightBottomPxSC->GetXaxis()->SetTitle("Px_reco - Px_corr (GeV)");
   lightBottomPxSC->GetYaxis()->SetTitle("Events");
   lightBottomPxSC->SetTitle("Bottom from H");
   lightBottomPxSC->DrawCopy("HIST");
   clightBottomPxSC->Write();
   clightBottomPxSC->ls();

   //compare TS and TC
   std::cout<<"drawing TS TC canvases"<<std::endl;
   lightAntiBottomPxTS->SetFillColor(mc1);
   lightAntiBottomPxTS->SetLineColor(mc1);
   lightAntiBottomPxTC->SetLineColor(mc2);
   clightAntiBottomPx->cd();
   lightAntiBottomPxTS->GetXaxis()->SetTitle("Px res (GeV)");
   lightAntiBottomPxTS->GetYaxis()->SetTitle("Events");
   lightAntiBottomPxTS->SetTitle("Anti Bottom from H");
   lightAntiBottomPxTS->SetMaximum( max( lightAntiBottomPxTS->GetMaximum(), lightAntiBottomPxTC->GetMaximum())+1 );
   lightAntiBottomPxTS->Draw("HIST");
   lightAntiBottomPxTC->Draw("SAMES");   
   moveStatsBox(lightAntiBottomPxTC);
   clightAntiBottomPx->Write();
   clightAntiBottomPx->ls();
   //compare SC
   std::cout<<"drawing SC canvases"<<std::endl;
   lightAntiBottomPxSC->SetFillColor(mc3);
   lightAntiBottomPxSC->SetLineColor(mc0);
   clightAntiBottomPxSC->cd();
   lightAntiBottomPxSC->GetXaxis()->SetTitle("Px_reco - Px_corr (GeV)");
   lightAntiBottomPxSC->GetYaxis()->SetTitle("Events");
   lightAntiBottomPxSC->SetTitle("Anti Bottom from H");
   lightAntiBottomPxSC->DrawCopy("HIST");
   clightAntiBottomPxSC->Write();
   clightAntiBottomPxSC->ls();

   std::cout<<"writing file"<<std::endl;



   outFile->cd();
   outFile->Write();
   outFile->Close(); 

   //c1.cd();
   //h_chi2.GetXaxis()->SetTitle("total #chi^{2} minimum");
   //h_chi2.Draw();
   //c1.SaveAs("minTotalChi2.pdf");
   
    //return 0;
   
}

void topReconstructionFromLHE::moveStatsBox(TH1F *hist){
    cout<<"in stats box func"<<endl;
    gPad->Update();
    cout<<"after gpad update"<<endl;
    TPaveStats *s = (TPaveStats*)hist->FindObject("stats");
    cout<<"1"<< endl;
*/ /*(    if (s == NULL) {
        //return;
        cout <<"null pointer"<<endl;
        return;
    }*/
// float x1 = s->GetX1NDC();
// float x2 = s->GetX2NDC();
/*    float y1 = s->GetY1NDC();
    float y2 = s->GetY2NDC();
cout<<"2"<<endl;
    s->SetY1NDC(y1 - (y2-y1) );
    s->SetY2NDC(y2 - (y2-y1) );
cout<<"3"<<endl;
}*/
