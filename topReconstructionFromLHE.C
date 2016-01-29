#define topReconstructionFromLHE_cxx
#include "topReconstructionFromLHE.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TString.h>
#include <TH1F.h>
#include <TCanvas.h>
//#include "WDaughterEllipseCalculator.h"
#include "topEventMinimizer.h"
#include <TPaveStats.h>
#include "Math/GenVector/LorentzVector.h"

using namespace std;
using namespace ROOT::Math;


void topReconstructionFromLHE::printVector(XYZTLorentzVector& v)
{
  cout << "px = " << v.Px() << endl;
  cout << "py = " << v.Py() << endl;
  cout << "pz = " << v.Pz() << endl;
  cout << "E  = " << v.E () << endl;
}

int main(){

    topReconstructionFromLHE t;
    t.Loop("output_files",0,1);

    return 0;
}

vector <string> topReconstructionFromLHE::particleNames;

void topReconstructionFromLHE::DeclareOutBranches(){
    outTree->Branch( "test", &top1SmearedPt);
    outTree->Branch( "testvec", &testvec);
    outTree->Branch( "stringtest", &(*evh_outside.smearedParticles["top"]) );
}

void topReconstructionFromLHE::DeclareMaps(){
    //typedef void (*FnPtr)(int, double, double, double, double);
    //getBestMap["top"] = topEventMinimizer::getTop;
}

void topReconstructionFromLHE::DeclareHists(){
    double lbound = -1;
    double rbound = 1;

    for ( vector<string>::iterator name = names.begin(); name != names.end(); name++){
        for ( vector<string>::iterator vartype = varTypes.begin(); vartype != varTypes.end(); vartype++){
            if (*vartype == "Pt") {lbound = -150; rbound = 150;}
            if (*vartype == "Eta") {lbound = -5; rbound = 5;}
            if (*vartype == "Phi") {lbound = -5; rbound = 5;}

            for ( vector<string>::iterator diftype = difTypes.begin(); diftype != difTypes.end(); diftype++){
                histdif[*diftype][*vartype][*name] = new TH1D( (*name + "_" + *vartype + "_" + *diftype).c_str(),
                        (*name + "_" + *vartype + "_" + *diftype).c_str(), 100, lbound, rbound);
            }
        }
    }
       
}

void topReconstructionFromLHE::DeclareCanvases(){
    for ( vector<string>::iterator name = names.begin(); name != names.end(); name++){
        for ( vector<string>::iterator vartype = varTypes.begin(); vartype != varTypes.end(); vartype++){
            canvasdif[*vartype][*name] = new TCanvas( ("c_" + *name + "_" + *vartype).c_str(), ("c_" + *name + "_" + *vartype).c_str(), 700, 700);
        }
    }
}

void topReconstructionFromLHE::FillHists(handleEvent evh, bool leptonFlag){
    for( hmap3::iterator h3 = histdif.begin(); h3 != histdif.end(); h3++){
        hmap2 histdif2 = h3->second;
        string diftype = h3->first;
            for( hmap2::iterator h2 = histdif2.begin(); h2 != histdif2.end(); h2++){
                hmap1 histdif1 = h2->second;
                string vartype = h2->first;
                for( hmap1::iterator h1 = histdif1.begin(); h1 != histdif1.end(); h1++){
                    //TH1D* hist = h1->second;
                    string name = h1->first;
                    string leptonFlagStr = static_cast<ostringstream*>( &(ostringstream() << leptonFlag) )->str();
                    string pname = "noname";

                    cout<< name <<" "<< vartype<<endl;

                    for ( vector< vector< string>>::iterator nameMapIt = nameMap.begin(); nameMapIt != nameMap.end(); nameMapIt++){
                        if (leptonFlagStr == nameMapIt->at(0) and name == nameMapIt->at(1)){
                            pname = nameMapIt->at(2);
                            cout<< pname <<endl;
                        }
                    }

                    if( diftype == "smearedTrue" ){
                        if( vartype == "Pt" ){
                            histdif[diftype][vartype][name]->Fill( evh.smearedParticles[pname]->Pt() - evh.trueParticles[pname]->Pt() );
                        } else if( vartype == "Eta" ){
                            histdif[diftype][vartype][name]->Fill( evh.smearedParticles[pname]->Eta() - evh.trueParticles[pname]->Eta() );
                        } else if( vartype == "Phi" ){
                            histdif[diftype][vartype][name]->Fill( evh.smearedParticles[pname]->Phi() - evh.trueParticles[pname]->Phi() );
                        }

                    } else if( diftype == "bestTrue" ){
                        if( vartype == "Pt" ){
                            histdif[diftype][vartype][name]->Fill( evh.bestParticles[pname]->Pt() - evh.trueParticles[pname]->Pt() );
                        } else if( vartype == "Eta" ){
                            histdif[diftype][vartype][name]->Fill( evh.bestParticles[pname]->Eta() - evh.trueParticles[pname]->Eta() );
                        } else if( vartype == "Phi" ){
                            histdif[diftype][vartype][name]->Fill( evh.bestParticles[pname]->Phi() - evh.trueParticles[pname]->Phi() );
                        }
                    }
                }
            }
    }

}

void topReconstructionFromLHE::PlotHists(){
    std::cout<<"Drawing histograms..."<<std::endl;
    int mc1 = 5;
    int mc2 = 4;
    int mc3 = 3;
    int mc0 = 0;
    string unit;

    for ( vector<string>::iterator name = names.begin(); name != names.end(); name++){
        for ( vector<string>::iterator vartype = varTypes.begin(); vartype != varTypes.end(); vartype++){
            if ( *vartype == "Pt" ) {unit = "(GeV)";}
            else {unit = "";}
            histdif["smearedTrue"][*vartype][*name]->SetFillColor(mc1);
            histdif["smearedTrue"][*vartype][*name]->SetLineColor(mc1);
            histdif["bestTrue"][*vartype][*name]->SetLineColor(mc2);

            histdif["smearedTrue"][*vartype][*name]->GetXaxis()->SetTitle( (*vartype + "Resolution " + unit).c_str() );
            histdif["smearedTrue"][*vartype][*name]->GetYaxis()->SetTitle( "Events" );
            histdif["smearedTrue"][*vartype][*name]->SetTitle( (*name + "_" + *vartype).c_str());
            histdif["smearedTrue"][*vartype][*name]->SetMaximum( max( histdif["smearedTrue"][*vartype][*name]->GetMaximum(), histdif["bestTrue"][*vartype][*name]->GetMaximum() ) +1 );

            canvasdif[*vartype][*name]->cd();

            histdif["smearedTrue"][*vartype][*name]->Draw("HIST");
            histdif["bestTrue"][*vartype][*name]->Draw("SAMES");

            moveStatsBox( histdif["smearedTrue"][*vartype][*name] );

            canvasdif[*vartype][*name]->Write();
            canvasdif[*vartype][*name]->ls();


        }
    }
  
}

void topReconstructionFromLHE::Loop(TString dir, int whichLoop, int maxLoops)
{
    //int whichLoop = 1; int maxLoops = 1;
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
//by  b_branchname->GetEntry(ientry); //read only this branch
    std::cout<<"before fchain"<<std::endl;
   if (fChain == 0) return;
   TRandom3 rand;
   rand.SetSeed(10);


   initOutput(dir,whichLoop);

   Long64_t nentries = fChain->GetEntriesFast();
   //Long64_t nentries = intree->GetEntries();
   
   Long64_t nbytes = 0, nb = 0;

   vector<XYZTLorentzVector> nonTopObjects;
   vector<double> nonTopObjectPtWidths;
   vector<double> nonTopObjectEtaWidths;
   vector<double> nonTopObjectPhiWidths;

   //gStyle->SetOptStat(101110);
   //TCanvas c1("c1","c1",500,500);
   //TH1D h_chi2("h_chi2","",100,0,100);
   std::cout<<"hi!"<<std::endl;

   //int jStart = whichLoop*(nentries/maxLoops) + ((whichLoop>(maxLoops-nentries%maxLoops))?(whichLoop+nentries%maxLoops-maxLoops):0);
int jStart = 0;
   //int jFinish = jStart + (nentries+whichLoop)/maxLoops;
   std::cout<<"number of entries = "<<nentries<<std::endl;
   int jFinish = 1000;

   std::cout<<jStart<<std::endl;
   std::cout<<"nentries = "<< nentries<<std::endl;
   std::cout<<jFinish<<std::endl;

    handleEvent evtemp;
    evh_outside = evtemp;

   DeclareHists();
   DeclareOutBranches();


   for (Long64_t jentry=jStart; jentry<jFinish; jentry++) {
     cout << "BEGINNING BRANCH NUMBER " << jentry << endl;
      Long64_t ientry = LoadTree(jentry);
      //intree->GetEntry(jentry);
      if (ientry < 0) break;
      //if (jentry > 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      XYZTLorentzVector smearedOtherLightParton;
      std::vector<XYZTLorentzVector> smearedOtherLightPartons;
      bool leptonFlag = 0; //1 if lepton and antineutrino (i.e. tbar branch is the leptonic one), 0 otherwise.

      //double METx(0.), METy(0.);
      int iJet = 0;

        //Declare structure (used as a container) to store all the particles involved in this event
        handleEvent evh;
cout<<"bottompx before set = "<<evh.smearedParticles["bottom"]->Px();


	  evh.smearedParticles["qFromW"]->SetPxPyPzE(P_X->at(3),P_Y->at(3),P_Z->at(3),E->at(3));
	  evh.smearedParticles["qbarFromW"]->SetPxPyPzE(P_X->at(4),P_Y->at(4),P_Z->at(4),E->at(4));
	  evh.smearedParticles["bottom"]->SetPxPyPzE(P_X->at(1),P_Y->at(1),P_Z->at(1),E->at(1));
	  evh.smearedParticles["antiBottom"]->SetPxPyPzE(P_X->at(2),P_Y->at(2),P_Z->at(2),E->at(2));
	  evh.smearedParticles["bFromH"]->SetPxPyPzE(P_X->at(5),P_Y->at(5),P_Z->at(5),E->at(5));
	  evh.smearedParticles["bbarFromH"]->SetPxPyPzE(P_X->at(6),P_Y->at(6),P_Z->at(6),E->at(6));
	  //smearedLepton.SetPxPyPzE(P_X->at(0),P_Y->at(0),P_Z->at(0),E->at(0));

	  evh.trueParticles["top"]->SetPxPyPzE(P_X->at(14),P_Y->at(14),P_Z->at(14),E->at(14));
	  evh.trueParticles["antiTop"]->SetPxPyPzE(P_X->at(15),P_Y->at(15),P_Z->at(15),E->at(15));
	  evh.trueParticles["Wplus"]->SetPxPyPzE(P_X->at(16),P_Y->at(16),P_Z->at(16),E->at(16));
	  evh.trueParticles["Wminus"]->SetPxPyPzE(P_X->at(17),P_Y->at(17),P_Z->at(17),E->at(17));
          //lepton.SetPxPyPzE(P_X->at(7),P_Y->at(7),P_Z->at(7),E->at(7));
	  //neutrino.SetPxPyPzE(P_X->at(18),P_Y->at(18),P_Z->at(18),E->at(18));
	  evh.trueParticles["bottom"]->SetPxPyPzE(P_X->at(8),P_Y->at(8),P_Z->at(8),E->at(8));
	  evh.trueParticles["antiBottom"]->SetPxPyPzE(P_X->at(9),P_Y->at(9),P_Z->at(9),E->at(9));
	  evh.trueParticles["qFromW"]->SetPxPyPzE(P_X->at(10),P_Y->at(10),P_Z->at(10),E->at(10));
	  evh.trueParticles["qbarFromW"]->SetPxPyPzE(P_X->at(11),P_Y->at(11),P_Z->at(11),E->at(11));
	  evh.trueParticles["bFromH"]->SetPxPyPzE(P_X->at(12),P_Y->at(12),P_Z->at(12),E->at(12));
	  evh.trueParticles["bbarFromH"]->SetPxPyPzE(P_X->at(13),P_Y->at(13),P_Z->at(13),E->at(13));
	  evh.trueParticles["higgs"]->SetPxPyPzE(P_X->at(19),P_Y->at(19),P_Z->at(19),E->at(19));


          std::cout<<"after set met"<<std::endl;
cout<<"bottompx after set = "<<evh.smearedParticles["bottom"]->Px();

        cout<<"Setting light partons"<<endl;
        for ( int ii=0; ii < ( (int)PID->size() - 21 ); ii++){
              smearedOtherLightParton.SetPxPyPzE(P_X->at(ii+21), P_Y->at(ii+21), P_Z->at(ii+21), E->at(ii+21) );
              smearedOtherLightPartons.push_back(smearedOtherLightParton);

          }
          std::cout<<"SIZEEEEEEEEEEEEEEEEEEE = " << smearedOtherLightPartons.size() << std::endl;

        XYZTLorentzVector leptonToAdd;
        XYZTLorentzVector leptonToAddGen;

        cout<<"Setting leptons and neutrinos"<<endl;
          if (PID->at(0) == 13) {
        	  evh.smearedParticles["lepton"]->SetPxPyPzE(P_X->at(0),P_Y->at(0),P_Z->at(0),E->at(0));
                  //define smeared anti neutrino as the MET; Pz set to zero; E set based on zero mass and Px and Py 
                  evh.smearedParticles["antiNeutrino"]->SetPxPyPzE(P_X->at(20),P_Y->at(20),0, sqrt( pow(P_X->at(20),2) + pow(P_Y->at(20),2) ) );
                  evh.trueParticles["lepton"]->SetPxPyPzE(P_X->at(7),P_Y->at(7),P_Z->at(7),E->at(7));
        	  evh.trueParticles["antiNeutrino"]->SetPxPyPzE(P_X->at(18),P_Y->at(18),P_Z->at(18),E->at(18));
                *(evh.smearedParticles["Wminus"]) = *(evh.smearedParticles["lepton"]) + *(evh.smearedParticles["antiNeutrino"]);
                *(evh.smearedParticles["Wplus"]) = *(evh.smearedParticles["qFromW"]) + *(evh.smearedParticles["qbarFromW"]);
                *(evh.smearedParticles["top"]) = *(evh.smearedParticles["Wplus"]) + *(evh.smearedParticles["bottom"]);
                *(evh.smearedParticles["antiTop"]) = *(evh.smearedParticles["Wminus"]) + *(evh.smearedParticles["antiBottom"]);

                  leptonFlag = 1;
                  leptonToAdd = *evh.smearedParticles["lepton"];
                  leptonToAddGen = *evh.trueParticles["lepton"];
          } else if (PID->at(0) == -13){
        	  evh.smearedParticles["antiLepton"]->SetPxPyPzE(P_X->at(0),P_Y->at(0),P_Z->at(0),E->at(0));
                  //define smeared neutrino as the MET; Pz set to zero; E set based on zero mass and Px and Py 
                  evh.smearedParticles["neutrino"]->SetPxPyPzE(P_X->at(20),P_Y->at(20),0, sqrt( pow(P_X->at(20),2) + pow(P_Y->at(20),2) ) );
                  evh.trueParticles["antiLepton"]->SetPxPyPzE(P_X->at(7),P_Y->at(7),P_Z->at(7),E->at(7));
        	  evh.trueParticles["neutrino"]->SetPxPyPzE(P_X->at(18),P_Y->at(18),P_Z->at(18),E->at(18));
                *(evh.smearedParticles["Wplus"]) = *(evh.smearedParticles["antiLepton"]) + *(evh.smearedParticles["neutrino"]);
                *(evh.smearedParticles["Wminus"]) = *(evh.smearedParticles["qFromW"]) + *(evh.smearedParticles["qbarFromW"]);
                *(evh.smearedParticles["top"]) = *(evh.smearedParticles["Wplus"]) + *(evh.smearedParticles["bottom"]);
                *(evh.smearedParticles["antiTop"]) = *(evh.smearedParticles["Wminus"]) + *(evh.smearedParticles["antiBottom"]);

                  leptonFlag = 0;
                  leptonToAdd = *evh.smearedParticles["antiLepton"];
                  leptonToAddGen = *evh.trueParticles["antiLepton"];
          }

          //Set smeared higgs by adding bFromH and bbarFromH
          *(evh.smearedParticles["higgs"]) = *(evh.smearedParticles["bFromH"]) + *(evh.smearedParticles["bbarFromH"]);

          std::cout<<"Setting non-top objects"<<endl;
      nonTopObjects.clear();
      nonTopObjectPtWidths.clear();
      nonTopObjectEtaWidths.clear();
      nonTopObjectPhiWidths.clear();
      
      //nonTopObjects.push_back( dynamic_cast<XYZTLorentzVector> (*(evh.smearedParticles["bFromH"])) );
      nonTopObjects.push_back( *(evh.smearedParticles["bFromH"]));
      nonTopObjectPtWidths.push_back(sqrt(evh.smearedParticles["bFromH"]->Pt()));
      nonTopObjectEtaWidths.push_back(0.01);
      nonTopObjectPhiWidths.push_back(0.01);

//      nonTopObjects.push_back( dynamic_cast<XYZTLorentzVector> (*(evh.smearedParticles["bbarFromH"])) );
      nonTopObjects.push_back( *(evh.smearedParticles["bbarFromH"]));
      nonTopObjectPtWidths.push_back(sqrt(evh.smearedParticles["bbarFromH"]->Pt()));
      nonTopObjectEtaWidths.push_back(0.01);
      nonTopObjectPhiWidths.push_back(0.01);

      for (int ii=0; ii<(int)smearedOtherLightPartons.size(); ii++){
        nonTopObjects.push_back(smearedOtherLightPartons.at(ii));
        nonTopObjectPtWidths.push_back( sqrt( (smearedOtherLightPartons.at(ii)).pt() ) );
//        nonTopObjectPtWidths.push_back(0.00001);
        nonTopObjectEtaWidths.push_back(0.01);
        nonTopObjectPhiWidths.push_back(0.01);
      }

      std::cout<<"SIZEEEEEEEEEEEEEEEEEEEEEEE = "<<nonTopObjects.size() << std::endl;

      std::cout<<"Setup top event minimizer"<<endl;
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
            ev.addLeptonicTop(evh.smearedParticles["bottom"]->Px(), evh.smearedParticles["bottom"]->Py(), evh.smearedParticles["bottom"]->Pz(), evh.smearedParticles["bottom"]->E(),
                    sqrt(evh.smearedParticles["bottom"]->Pt()), sigmaEtaJet, sigmaPhiJet,
                    evh.smearedParticles["antiLepton"]->Px(), evh.smearedParticles["antiLepton"]->Py(), evh.smearedParticles["antiLepton"]->Pz(), evh.smearedParticles["antiLepton"]->E(),
                    sigmaPtLep, sigmaEtaLep, sigmaPhiLep,
                    mTop, sigmaMTop, mW, sigmaMW);

            ev.addHadronicTop(evh.smearedParticles["antiBottom"]->Px(), evh.smearedParticles["antiBottom"]->Py(), evh.smearedParticles["antiBottom"]->Pz(), evh.smearedParticles["antiBottom"]->E(),
                    sqrt(evh.smearedParticles["antiBottom"]->Pt()), sigmaEtaJet, sigmaPhiJet,
                    evh.smearedParticles["qFromW"]->Px(), evh.smearedParticles["qFromW"]->Py(), evh.smearedParticles["qFromW"]->Pz(), evh.smearedParticles["qFromW"]->E(),
                    sqrt(evh.smearedParticles["qFromW"]->Pt()), sigmaEtaJet, sigmaPhiJet,
                    evh.smearedParticles["qbarFromW"]->Px(), evh.smearedParticles["qbarFromW"]->Py(), evh.smearedParticles["qbarFromW"]->Pz(), evh.smearedParticles["qbarFromW"]->E(),
                    sqrt(evh.smearedParticles["qbarFromW"]->Pt()), sigmaEtaJet, sigmaPhiJet,
                    mTop, sigmaMTop, mW, sigmaMW);


      }


      if (leptonFlag == 1){
            ev.addHadronicTop(evh.smearedParticles["bottom"]->Px(), evh.smearedParticles["bottom"]->Py(), evh.smearedParticles["bottom"]->Pz(), evh.smearedParticles["bottom"]->E(),
                    sqrt(evh.smearedParticles["bottom"]->Pt()), sigmaEtaJet, sigmaPhiJet,
                    evh.smearedParticles["qFromW"]->Px(), evh.smearedParticles["qFromW"]->Py(), evh.smearedParticles["qFromW"]->Pz(), evh.smearedParticles["qFromW"]->E(),
                    sqrt(evh.smearedParticles["qFromW"]->Pt()), sigmaEtaJet, sigmaPhiJet,
                    evh.smearedParticles["qbarFromW"]->Px(), evh.smearedParticles["qbarFromW"]->Py(), evh.smearedParticles["qbarFromW"]->Pz(), evh.smearedParticles["qbarFromW"]->E(),
                    sqrt(evh.smearedParticles["qbarFromW"]->Pt()), sigmaEtaJet, sigmaPhiJet,
                    mTop, sigmaMTop, mW, sigmaMW);

            ev.addLeptonicTop(evh.smearedParticles["antiBottom"]->Px(), evh.smearedParticles["antiBottom"]->Py(), evh.smearedParticles["antiBottom"]->Pz(), evh.smearedParticles["antiBottom"]->E(),
                    sqrt(evh.smearedParticles["antiBottom"]->Pt()), sigmaEtaJet, sigmaPhiJet,
                    evh.smearedParticles["lepton"]->Px(), evh.smearedParticles["lepton"]->Py(), evh.smearedParticles["lepton"]->Pz(), evh.smearedParticles["lepton"]->E(),
                    sigmaPtLep, sigmaEtaLep, sigmaPhiLep,
                    mTop, sigmaMTop, mW, sigmaMW);


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

    //NEW Get best values
       std::cout<<"Getting best values"<<endl;

        *evh.bestParticles["top"] = ev.getConverter("getTop", 0);
        *evh.bestParticles["bottom"] = ev.getConverter("getBJet", 0);
        *evh.bestParticles["Wplus"] = ev.getConverter("getW", 0);
        if (leptonFlag == 0){
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
        if (leptonFlag == 0){
            *evh.bestParticles["qFromW"] = ev.getConverter("getWDaughter1", 1);
            *evh.bestParticles["qbarFromW"] = ev.getConverter("getWDaughter2", 1);
        } else {
            *evh.bestParticles["lepton"] = ev.getConverter("getWDaughter1", 1);
            *evh.bestParticles["antiNeutrino"] = ev.getConverter("getWDaughter2", 1);
        }
        *evh.bestParticles["bbarFromH"] = ev.getConverter("getNonTopObject4", 1);

        //evh.bestParticles["higgs"]->SetPxPyPzE( evh.bestParticles["bFromH"]->Px() + evh.bestParticles["bbarFromH"]->Px(), evh.bestParticles["bFromH"]->Py() + evh.bestParticles["bbarFromH"]->Py(), evh.bestParticles["bFromH"]->Pz() + evh.bestParticles["bbarFromH"]->Pz(), evh.bestParticles["bFromH"]->E() + evh.bestParticles["bbarFromH"]->E() );
            
        //Set best higgs by adding bFromH and bbarFromH
          *(evh.bestParticles["higgs"]) = *(evh.bestParticles["bFromH"]) + *(evh.bestParticles["bbarFromH"]);
 


        //Fill Hists
        std::cout<<"Filling hists"<<endl;
        if( innerMinStatus == 0 and outerMinStatus == 0 ){
            FillHists(evh, leptonFlag);
        }

    cout <<"Printing P's"<<endl;
    cout << evh.smearedParticles["top"]->Px()<<endl;
    //cout << evh.smearedParticles["top"]->Pt()<<endl;

    evh_outside = evh;
    cout << "outside "<< evh_outside.smearedParticles["top"]->Px()<<endl;
    testvec.SetPxPyPzE(11.4, 2,3,4);

    outTree->Fill();

   } //end loop over events

   cout<<"Opening outFile"<<endl;
    outFile->cd();
    DeclareCanvases();
    PlotHists();
    outFile->Write();
    outFile->Close();

    cout<<"Done!"<<endl;
   
}

//void topReconstructionFromLHE::getBestObjects(){}

void topReconstructionFromLHE::moveStatsBox(TH1D *hist){
    cout<<"in stats box func"<<endl;
    gPad->Update();
    cout<<"after gpad update"<<endl;
    TPaveStats *s = (TPaveStats*)hist->FindObject("stats");
    cout<<"1"<< endl;
/*(    if (s == NULL) {
        //return;
        cout <<"null pointer"<<endl;
        return;
    }*/
    //float x1 = s->GetX1NDC();
    //float x2 = s->GetX2NDC();
    float y1 = s->GetY1NDC();
    float y2 = s->GetY2NDC();
cout<<"2"<<endl;
    s->SetY1NDC(y1 - (y2-y1) );
    s->SetY2NDC(y2 - (y2-y1) );
cout<<"3"<<endl;
}
