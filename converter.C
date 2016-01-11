#define converter_cxx
#include "converter.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TString.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TPaveStats.h>
#include <TLorentzVector.h>

using namespace std;
using namespace ROOT::Math;



int main(){

    converter t;
    //t.Loop("output_files",0,1);
    t.Loop();

    return 0;
}

void converter::Filler(int i, int j){
    PID.at(i) = Pid.at(j);
    P_X.at(i) = Px.at(j);
    P_Y.at(i) = Py.at(j);
    P_Z.at(i) = Pz.at(j);
    E.at(i) = En.at(j);
    M.at(i) = Ma.at(j);

    return;

}

void converter::Smearer( TLorentzVector &jet ){
    //cout<<"in smearer"<<endl;
    double pt = jet.Pt();
    double phi = jet.Phi();
    double eta = jet.Eta();
    double m = jet.M();
    double newpt = pt + sqrt(pt)*rand.Gaus();
    double newphi = phi + 0.01*rand.Gaus();
    jet.SetPtEtaPhiM( newpt, eta, newphi, m );
    cout << "Mass " << m << " " << jet.M() <<endl;
    return;
}

void converter::SmearedFiller( int i, TLorentzVector &jet, int pid ){
    //cout<<"in smearedfiller"<<endl;
    PID.at(i) = pid;
    P_X.at(i) = jet.Px();
    P_Y.at(i) = jet.Py();
    P_Z.at(i) = jet.Pz();
    E.at(i) = jet.Energy();
    M.at(i) = jet.M();
    return;
}



void converter::Loop()
{
//   In a ROOT session, you can do:
//      root> .L converter.C
//      root> converter t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
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

   TFile *outfile = new TFile("./fixedsinglemuntuple.root","RECREATE","tree");
   TTree *outtree = new TTree("physics","physics");

   //vector<double> P_X, P_Y, P_Z, E, M;
   //vector<int> PID;

   outtree->Branch("P_X", &P_X);
   outtree->Branch("P_Y", &P_Y);
   outtree->Branch("P_Z", &P_Z);
   outtree->Branch("E", &E);
   outtree->Branch("M", &M);
   outtree->Branch("PID", &PID);
   outtree->Branch("run", &run);
   outtree->Branch("lumi", &lumi);
   outtree->Branch("event", &event);
   outtree->Branch("status", &status);
   outtree->Branch("n_particles", &n_particles);

   for (int i=0; i<21; i++){
       P_X.push_back(0);
       P_Y.push_back(0);
       P_Z.push_back(0);
       E.push_back(0);
       M.push_back(0);
       PID.push_back(-9999);
   }


   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

    Pid.clear();
    Px.clear();
    Py.clear();
    Pz.clear();
    En.clear();
    Ma.clear();
    mother1.clear();
    mother2.clear();

    if (Particle_size > 15){
        cout<<"particle_size = " <<Particle_size<<endl;
    }

    for (int ipart = 0; ipart < Particle_size; ipart++){
        Pid.push_back(Particle_PID[ipart]);
        Px.push_back(Particle_Px[ipart]);
        Py.push_back(Particle_Py[ipart]);
        Pz.push_back(Particle_Pz[ipart]);
        En.push_back(Particle_E[ipart]);
        Ma.push_back(Particle_M[ipart]);
        mother1.push_back(Particle_Mother1[ipart]);
        mother2.push_back(Particle_Mother2[ipart]);
        //cout<<"mother1 = "<<Particle_Mother1[ipart]<<endl;
        //cout<<"mother2 = "<<Particle_Mother2[ipart]<<endl;
    }

    //cout << "Pid = "<<Particle_PID[0]<<endl;

    numpar = Particle_size;
    n_particles = 7;
    status = 1;
    event = Event_;
    lumi = 0;
    run = 1;



    TLorentzVector bJetFromTop, bbarJetFromTop, bJetFromH, bbarJetFromH, qFromW, qbarFromW, muon;
    double leptonPx, leptonPy;
    int qFromWPID, qbarFromWPID;
//cout<<"haha"<<endl;
      int counter = 0;

      for (int ivec=0; ivec < numpar; ivec++){
        if (Pid.at(ivec) == 5){
            if (Pid.at( mother1.at(ivec) ) == 6 or Pid.at( mother2.at(ivec) ) == 6 ){
                Filler( 8, ivec );
                counter++;
                //cout<<"pid1"<<endl;
                bJetFromTop.SetPxPyPzE( Px.at(ivec), Py.at(ivec), Pz.at(ivec), En.at(ivec) );
            }
            if (Pid.at( mother1.at(ivec) ) == 25 or Pid.at( mother2.at(ivec) ) == 25 ){
                Filler( 12, ivec );
                counter++;
                 //cout<<"pid2"<<endl;
                bJetFromH.SetPxPyPzE( Px.at(ivec), Py.at(ivec), Pz.at(ivec), En.at(ivec) );
             }
        } else if (Pid.at(ivec) == -5){
            if (Pid.at( mother1.at(ivec) ) == -6 or Pid.at( mother2.at(ivec) ) == -6 ){
                Filler( 9, ivec );
                counter++;
                //cout<<"pid3"<<endl;
                 bbarJetFromTop.SetPxPyPzE( Px.at(ivec), Py.at(ivec), Pz.at(ivec), En.at(ivec) );
             }
            if (Pid.at( mother1.at(ivec) ) == 25 or Pid.at( mother2.at(ivec) ) == 25 ){
                Filler( 13, ivec );
                counter++;
                //cout<<"pid4"<<endl;
                bbarJetFromH.SetPxPyPzE( Px.at(ivec), Py.at(ivec), Pz.at(ivec), En.at(ivec) );
             }
        } else if (Pid.at(ivec) == 13){
            Filler( 7, ivec );
            counter++;
                //cout<<"pid5"<<endl;
            muon.SetPxPyPzE( Px.at(ivec), Py.at(ivec), Pz.at(ivec), En.at(ivec) );
            SmearedFiller( 0, muon, 13 );
            leptonPx = Px.at(ivec);
            leptonPy = Py.at(ivec);
        } else if (Pid.at(ivec) == 1 or Pid.at(ivec) == 2 or Pid.at(ivec) == 3 or Pid.at(ivec) == 4 ){
            if (mother1.at(ivec)>=0 and mother1.at(ivec)<numpar and mother2.at(ivec)>=0 and mother2.at(ivec)<numpar){
            if (Pid.at( mother1.at(ivec) ) == 24 or Pid.at( mother2.at(ivec) ) == 24 or Pid.at(mother1.at(ivec) )==-24 or Pid.at(mother2.at(ivec))==-24 ){
                Filler ( 10, ivec );
                counter++;
                //cout<<"pid6"<<endl;
                qFromW.SetPxPyPzE( Px.at(ivec), Py.at(ivec), Pz.at(ivec), En.at(ivec) );
                qFromWPID = Pid.at(ivec);
            }
            }
         } else if (Pid.at(ivec) == -1 or Pid.at(ivec) == -2 or Pid.at(ivec) == -3 or Pid.at(ivec) == -4 ){
            if (mother1.at(ivec)>=0 and mother1.at(ivec)<numpar and mother2.at(ivec)>=0 and mother2.at(ivec)<numpar){
            if (Pid.at( mother1.at(ivec) ) == 24 or Pid.at( mother2.at(ivec) ) == 24 or Pid.at(mother1.at(ivec) )==-24 or Pid.at(mother2.at(ivec))==-24 ){
                Filler (11, ivec);
                counter++;
                //cout<<"pid7"<<endl;
                qbarFromW.SetPxPyPzE( Px.at(ivec), Py.at(ivec), Pz.at(ivec), En.at(ivec) );
                qbarFromWPID = Pid.at(ivec);
            }
            }
         } else if (Pid.at(ivec) == -14){
            Filler (18, ivec);
            counter++;
                //cout<<"pid8"<<endl;
        } else if (Pid.at(ivec) == 6){
            Filler (14, ivec);
            counter++;
                //cout<<"pid9"<<endl;
        } else if (Pid.at(ivec) == -6){
            Filler (15, ivec);
            counter++;
                 //cout<<"pid10"<<endl;
       } else if (Pid.at(ivec) == 24){
            Filler (16, ivec);
            counter++;
                 //cout<<"pid11"<<endl;
       } else if (Pid.at(ivec) == -24){
            Filler (17, ivec);
            counter++;
                 //cout<<"pid12"<<endl;
       } else if (Pid.at(ivec) == 25){
            Filler (19, ivec);
            counter++;
                 //cout<<"pid13"<<endl;
       }
      }
      //cout<<"hehe"<<endl;
      if (counter != 13){
          cout <<"counter != 13" <<endl;
          cout<<counter<<endl;
          cout<<"jentry = "<<jentry<<endl;
          continue;
      }


      //cout<<"hoho"<<endl;
    //start smearing
    //double newpt = 0;
    Smearer( bJetFromTop );
    Smearer( bbarJetFromTop );
    Smearer( qFromW );
    Smearer( qbarFromW );
    Smearer( bJetFromH );
    Smearer( bbarJetFromH );

    SmearedFiller( 1, bJetFromTop, 5 );
    SmearedFiller( 2, bbarJetFromTop, -5 );
    SmearedFiller( 3, qFromW, qFromWPID );
    SmearedFiller( 4, qbarFromW, qbarFromWPID );
    SmearedFiller( 5, bJetFromH, 5 );
    SmearedFiller( 6, bbarJetFromH, -5 );

    double metPx = -( bJetFromTop.Px() + bbarJetFromTop.Px() + bJetFromH.Px() + bbarJetFromH.Px() + qFromW.Px() + qbarFromW.Px() + leptonPx );
    double metPy = -( bJetFromTop.Py() + bbarJetFromTop.Py() + bJetFromH.Py() + bbarJetFromH.Py() + qFromW.Py() + qbarFromW.Py() + leptonPy );

    P_X.at(20) = metPx;
    P_Y.at(20) = metPy;
    P_Z.at(20) = 0;
    PID.at(20) = 9999;
    E.at(20) = 0;
    M.at(20) = 0;

    //cout<<"endloop"<<endl;
  
   outtree->Fill();

   }

   outfile->Write();
    outfile->Close();
    return;

   
}
