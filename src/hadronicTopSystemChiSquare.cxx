#include "hadronicTopSystemChiSquare.h"

using namespace commonstruct;

hadronicTopSystemChiSquare::hadronicTopSystemChiSquare( top_system &topsystem )
    : topSystemChiSquare( topsystem )
{
}


hadronicTopSystemChiSquare::~hadronicTopSystemChiSquare()
{
}

void hadronicTopSystemChiSquare::printTopConstituents()
{
    cout << "Hadronic top decay products:" << endl;
    //  cout << "In the derived class:" << endl;
    cout
        << "b-jet: "
        << "\npx = " << topsys.calc.b_px() << "\npy = " << topsys.calc.b_py() << "\npz = " << topsys.calc.b_pz()
        << "\ne  = " << topsys.calc.b_e()
        //       << "\nm  = " <<
        //       sqrt(max(0.,topsys.calc.b_e()*topsys.calc.b_e()-topsys.calc.b_px()*topsys.calc.b_px()-topsys.calc.b_py()*topsys.calc.b_py()-topsys.calc.b_pz()*topsys.calc.b_pz()))
        << endl;

    cout
        << "first light quark: "
        << "\npx = " << topsys.calc.Wd1_px() << "\npy = " << topsys.calc.Wd1_py()
        << "\npz = " << topsys.calc.Wd1_pz() << "\ne  = " << topsys.calc.Wd1_e()
        //       << "\nm  = " <<
        //       sqrt(max(0.,topsys.calc.Wd1_e()*topsys.calc.Wd1_e()-topsys.calc.Wd1_px()*topsys.calc.Wd1_px()-topsys.calc.Wd1_py()*topsys.calc.Wd1_py()-topsys.calc.Wd1_pz()*topsys.calc.Wd1_pz()))
        << endl;

    cout
        << "second light quark: "
        << "\npx = " << topsys.vars.Wd2_px << "\npy = " << topsys.vars.Wd2_py
        << "\npz = " << topsys.vars.Wd2_pz << "\ne  = " << topsys.calc.Wd2_e()
        //       << "\nm  = " <<
        //       sqrt(max(0.,WDaughter2E_*WDaughter2E_-topsys.calc.Wd2_px()*topsys.calc.Wd2_px()-WDaughter2Py_*WDaughter2Py_-WDaughter2Pz_*WDaughter2Pz_))
        //       << "\nreco m = " << sqrt(reconstructed_WDaughter2Mass2_)
        << endl;

    cout << "top mass: " << topsys.calc.mTop() << endl;

    cout << "W mass: " << topsys.calc.mW() << endl;

}

void hadronicTopSystemChiSquare::setEllipseAngle()
{
    //theta_ = theta;
    resetWDaughter2();
    //calcWDaughter2Deltas();
    //calcTopMomentum();
}

void hadronicTopSystemChiSquare::printWDaughter2()
{
    cout << "Current second light quark: "
         << "\npt  = " << topsys.calc.Wd2_pt() << "\nphi = " << topsys.calc.Wd2_phi()
         << "\neta = " << topsys.calc.Wd2_eta() << endl;

    cout << "theta is " << topsys.vars.theta << endl;

    cout << "dif pt  = " << topsys.calc.Wd2_dif_pt()
         << "\ndif phi = " << topsys.calc.Wd2_dif_phi()
         << "\ndif eta = " << topsys.calc.Wd2_dif_eta() << endl;

}


