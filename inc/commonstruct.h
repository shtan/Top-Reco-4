#ifndef commonstruct_h
#define commonstruct_h

#include <vector>
#include "Math/GenVector/LorentzVector.h"
#include "TLorentzVector.h"

using namespace std;

namespace commonstruct
{

//struct eventInfo {
    

//};

struct top_system {
/*
    //INPUT -----------------------------------------------------------
    //b jet
    const double b_pt_input, b_eta_input, b_phi_input, b_e_input;
    const double b_pt_width, b_eta_width, b_phi_width, b_e_width;

    //W Daughter 1
    const double Wd1_pt_input, Wd1_eta_input, Wd1_phi_input, Wd1_e_input;
    const double Wd1_pt_width, Wd1_eta_width, Wd1_phi_width, Wd1_e_width;

    //W Daughter 2
    const double Wd2_pt_input, Wd2_eta_input, Wd2_phi_input, Wd2_e_input;
    const double Wd2_pt_width, Wd2_eta_width, Wd2_phi_width, Wd2_e_width;

    //Central value and width for Breit Wigner function for Top and W masses
    const double mTop_central, mW_central, mTop_width, mW_width;
*/

    //VARIABLES THAT ARE CHANGED DURING MINIMISATION ----------------------

    public:

        struct inputs {
            //b jet
            const double b_pt, b_eta, b_phi, b_e;
            const double b_pt_width, b_eta_width, b_phi_width;

            //W Daughter 1
            const double Wd1_pt, Wd1_phi, Wd1_eta, Wd1_e;
            const double Wd1_pt_width, Wd1_eta_width, Wd1_phi_width;

            //W Daughter 2
            const double Wd2_pt, Wd2_phi, Wd2_eta, Wd2_e;
            const double Wd2_pt_width, Wd2_eta_width, Wd2_phi_width;

            //Central value and width for Breit Wigner function for Top and W masses
            const double mTop_central, mW_central, mTop_width, mW_width;

            //Leptonic flag
            const bool leptonic;
           
            
/*            const double b_px;
            const double b_py = 7;
            const double b_pz = 8;
            const double b_e = 3;
            double btest = 4;*/
        /*    const double b_pt_width, b_eta_width, b_phi_width, b_e_width;

            //W Daughter 1
            const double Wd1_px, Wd1_py, Wd1_pz, Wd1_e;
            const double Wd1_pt_width, Wd1_eta_width, Wd1_phi_width, Wd1_e_width;

            //W Daughter 2
            const double Wd2_px, Wd2_py, Wd2_pz, Wd2_e;
            const double Wd2_pt_width, Wd2_eta_width, Wd2_phi_width, Wd2_e_width;

            //Central value and width for Breit Wigner function for Top and W masses
            const double mTop_central, mW_central, mTop_width, mW_width;
        */
           // inputs( const double &inputconst2 ) : b_px(inputconst2) {}

            inputs( const double &bPt, const double &bEta, const double &bPhi, const double &bE,
                    const double &bPtWidth, const double &bEtaWidth, const double &bPhiWidth,
                    const double &Wd1Pt, const double &Wd1Eta, const double &Wd1Phi, const double &Wd1E,
                    const double &Wd1PtWidth, const double &Wd1EtaWidth, const double &Wd1PhiWidth,
                    const double &Wd2Pt, const double &Wd2Eta, const double &Wd2Phi, const double &Wd2E,
                    const double &Wd2PtWidth, const double &Wd2EtaWidth, const double &Wd2PhiWidth,
                    const double &mTopCentral, const double &mTopWidth, const double &mWCentral, const double &mWWidth)
                :   b_pt(bPt), b_eta(bEta), b_phi(bPhi), b_e(bE),
                    b_pt_width(bPtWidth), b_eta_width(bEtaWidth), b_phi_width(bPhiWidth),
                    Wd1_pt(Wd1Pt), Wd1_eta(Wd1Eta), Wd1_phi(Wd1Phi), Wd1_e(Wd1E),
                    Wd1_ptWidth(Wd1PtWidth), Wd1_etaWidth(Wd1EtaWidth), Wd1_phiWidth(Wd1PhiWidth)                    
                    Wd2_pt(Wd2Pt), Wd2_eta(Wd2Eta), Wd2_phi(Wd2Phi), Wd2_e(Wd2E),
                    Wd2_ptWidth(Wd2PtWidth), Wd2_etaWidth(Wd2EtaWidth), Wd2_phiWidth(Wd2PhiWidth)
                    mTop_central(mTopCentral), mTop_width(mTopWidth), mW_central(mWCentral), mW_width(mWWidth) {}

        };

        inputs input;


        struct variables {
            inputs &input_;

            public:
                variables(inputs & bip) : input_(bip) {}
                double b_delta_pt;

                double b_pt(){ return b_delta_pt + input_.b_px; }
        };

        variables vars;

/*        struct bestfit{
            inputs &input_;
            variables &vars_;

            public:
                bestfit(inputs & bip, variables & var) : input_(bip), vars_(var) {}
                double lars(){ return input_.b_px + vars_.b_pt(); }

                double bestchi;
        };

        bestfit best;*/

        //top_system( const double &inputconst ) : input( inputconst ), vars(input), best(input, vars) {}

        top_system( const double &bPt, const double &bEta, const double &bPhi, const double &bE,
                    const double &bPtWidth, const double &bEtaWidth, const double &bPhiWidth,
                    const double &Wd1Pt, const double &Wd1Eta, const double &Wd1Phi, const double &Wd1E,
                    const double &Wd1PtWidth, const double &Wd1EtaWidth, const double &Wd1PhiWidth,
                    const double &Wd2Pt, const double &Wd2Eta, const double &Wd2Phi, const double &Wd2E,
                    const double &Wd2PtWidth, const double &Wd2EtaWidth, const double &Wd2PhiWidth,
                    const double &mTopCentral, const double &mTopWidth, const double &mWCentral, const double &mWWidth)
                :   input( bPt, bEta, bPhi, bE,
                    bPtWidth, bEtaWidth, bPhiWidth,
                    Wd1Pt, Wd1Eta, Wd1Phi, Wd1E,
                    Wd1PtWidth, Wd1EtaWidth, Wd1PhiWidth                    
                    Wd2Pt, Wd2Eta, Wd2Phi, Wd2E,
                    Wd2PtWidth, Wd2EtaWidth, Wd2PhiWidth
                    mTopCentral, mTopWidth, mWCentral, mWWidth ),
                    vars(input) {}

}; //END top_system definition

struct hadronic_top_system : public top_system
{
    public:
        double bb2;
        const double blah;
        const double optional1_;

        hadronic_top_system( const double &inputconst, const double th, const double optional1 = -99 ) : top_system( inputconst ), blah(th), optional1_(optional1) {}

};

//double adder(top_system &);
//double tester();

//inline double adder(top_system & cha){ return cha.input.b_px + cha.vars.b_pt() + cha.best.lars(); }
inline double adder(top_system & cha){
    TLorentzVector bob;
    bob.SetPxPyPzE(cha.input.b_px, cha.input.b_py, cha.input.b_pz, cha.input.b_e);
    return bob.Pt();
}
    





}

#endif
