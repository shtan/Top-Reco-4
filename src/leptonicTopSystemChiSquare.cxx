#include "leptonicTopSystemChiSquare.h"

using namespace commonstruct;

leptonicTopSystemChiSquare::leptonicTopSystemChiSquare(
        top_system & topsystem)
    : topSystemChiSquare(topsystem) {}

leptonicTopSystemChiSquare::~leptonicTopSystemChiSquare()
{
}

void leptonicTopSystemChiSquare::printTopConstituents()
{
    cout << "Leptonic top decay products:" << endl;
    cout
        << "b-jet: "
        << "\npx = " << topsys.calc.b_px() << "\npy = " << topsys.calc.b_py() << "\npz = " << topsys.calc.b_pz()
        << "\ne  = " << topsys.calc.b_e()
        //       << "\nm  = " <<
        //       sqrt(max(0.,bJetE_*bJetE_-bJetPx_*bJetPx_-bJetPy_*bJetPy_-bJetPz_*bJetPz_))
        << endl;

    cout
        << "lepton: "
        << "\npx = " << topsys.calc.Wd1_px() << "\npy = " << topsys.calc.Wd1_px()
        << "\npz = " << topsys.calc.Wd1_pz() << "\ne  = " << topsys.calc.Wd1_e()
        //       << "\nm  = " <<
        //       sqrt(max(0.,WDaughter1E_*WDaughter1E_-WDaughter1Px_*WDaughter1Px_-WDaughter1Py_*WDaughter1Py_-WDaughter1Pz_*WDaughter1Pz_))
        << endl;

}

void leptonicTopSystemChiSquare::setEllipseAngle()
{
    //theta_ = theta;
    resetWDaughter2();
    //calcTopMomentum();
}
