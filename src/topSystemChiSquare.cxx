#include "topSystemChiSquare.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TH1D.h"
#include "Math/QuantFuncMathCore.h"

using namespace commonstruct;

topSystemChiSquare::topSystemChiSquare( top_system & topsystem )
    : topsys(topsystem), WDaughter2Calc_( topsys )
{
    //setupWDaughter2Ellipse();
}

/*topSystemChiSquare::topSystemChiSquare(
    const double &bJetPx, const double &bJetPy, const double &bJetPz,
    const double &bJetE, const double &bJetPtWidth, const double &bJetEtaWidth,
    const double &bJetPhiWidth, const double &WDaughter1Px,
    const double &WDaughter1Py, const double &WDaughter1Pz,
    const double &WDaughter1E, const double &WDaughter1PtWidth,
    const double &WDaughter1EtaWidth, const double &WDaughter1PhiWidth,
    const double &mTop, const double &sigmaMTop, const double &mW,
    const double &sigmaMW, const double &mWDaughter2, top_system &topsys)
    : bJetPx_(bJetPx), bJetPy_(bJetPy), bJetPz_(bJetPz), bJetE_(bJetE),
      WDaughter1Px_(WDaughter1Px), WDaughter1Py_(WDaughter1Py),
      WDaughter1Pz_(WDaughter1Pz), WDaughter1E_(WDaughter1E), mTop_(mTop),
      mW_(mW), mWDaughter2_(mWDaughter2), bJetPtWidth_(bJetPtWidth),
      bJetEtaWidth_(bJetEtaWidth), bJetPhiWidth_(bJetPhiWidth),
      WDaughter1PtWidth_(WDaughter1PtWidth),
      WDaughter1EtaWidth_(WDaughter1EtaWidth),
      WDaughter1PhiWidth_(WDaughter1PhiWidth), sigmaMTop_(sigmaMTop),
      sigmaMW_(sigmaMW), rangeFlag_(false), hasHighEdge_(false),
      mTopEdgeLow_(-1.), mTopEdgeHigh_(-1.), deltaMTopRangeLow_(-1.),
      deltaMTopRangeHigh_(-1.),
      //WDaughter2Calc_(bJetPx_, bJetPy_, bJetPz_, bJetE_, WDaughter1Px_,
                      //WDaughter1Py_, WDaughter1Pz_, WDaughter1E_, mTop_, mW_,
                      //mWDaughter2_),
      WDaughter2Calc_(topsys),
      WDaughter2Px_(0.), WDaughter2Py_(0.), WDaughter2Pz_(0.), WDaughter2E_(0.),
      WDaughter2Pt_(0.), WDaughter2Phi_(0.), WDaughter2Eta_(0.), topPx_(0.),
      topPy_(0.), topPz_(0.), topE_(0.)
{
    // cout << "Setting up a top system chi square object" << endl;

    setBJet(bJetPx, bJetPy, bJetPz, bJetE);
    setObservedBJet(bJetPx, bJetPy, bJetPz, bJetE);
    setBJetDeltas(0., 0., 0.);

    setWDaughter1(WDaughter1Px, WDaughter1Py, WDaughter1Pz, WDaughter1E);
    setObservedWDaughter1(WDaughter1Px, WDaughter1Py, WDaughter1Pz,
                          WDaughter1E);
    setWDaughter1Deltas(0., 0., 0.);

    setTopMassDelta(0.);
    setWMassDelta(0.);

    setupWDaughter2Ellipse();
}

topSystemChiSquare::topSystemChiSquare(const topSystemChiSquare &other)
    : bJetPx_(other.bJetPx_), bJetPy_(other.bJetPy_), bJetPz_(other.bJetPz_),
      bJetE_(other.bJetE_), WDaughter1Px_(other.WDaughter1Px_),
      WDaughter1Py_(other.WDaughter1Py_), WDaughter1Pz_(other.WDaughter1Pz_),
      WDaughter1E_(other.WDaughter1E_), mTop_(other.mTop_), mW_(other.mW_),
      mWDaughter2_(other.mWDaughter2_), bJetPtWidth_(other.bJetPtWidth_),
      bJetEtaWidth_(other.bJetEtaWidth_), bJetPhiWidth_(other.bJetPhiWidth_),
      WDaughter1PtWidth_(other.WDaughter1PtWidth_),
      WDaughter1EtaWidth_(other.WDaughter1EtaWidth_),
      WDaughter1PhiWidth_(other.WDaughter1PhiWidth_),
      sigmaMTop_(other.sigmaMTop_), sigmaMW_(other.sigmaMW_),
      rangeFlag_(other.rangeFlag_), hasHighEdge_(other.hasHighEdge_),
      mTopEdgeLow_(other.mTopEdgeLow_), mTopEdgeHigh_(other.mTopEdgeHigh_),
      deltaMTopRangeLow_(other.deltaMTopRangeLow_),
      deltaMTopRangeHigh_(other.deltaMTopRangeHigh_),
      topsys(other.topsys),
      //WDaughter2Calc_(other.bJetPx_, other.bJetPy_, other.bJetPz_, other.bJetE_,
                      //other.WDaughter1Px_, other.WDaughter1Py_,
                      //other.WDaughter1Pz_, other.WDaughter1E_, other.mTop_,
                      //other.mW_, other.mWDaughter2_),
      WDaughter2Calc_(topsys),
      WDaughter2Px_(other.WDaughter2Px_), WDaughter2Py_(other.WDaughter2Py_),
      WDaughter2Pz_(other.WDaughter2Pz_), WDaughter2E_(other.WDaughter2E_),
      WDaughter2Pt_(other.WDaughter2Pt_), WDaughter2Phi_(other.WDaughter2Phi_),
      WDaughter2Eta_(other.WDaughter2Eta_), topPx_(other.topPx_),
      topPy_(other.topPy_), topPz_(other.topPz_), topE_(other.topE_)

{
    setBJetDeltas(0., 0., 0.);
    setWDaughter1Deltas(0., 0., 0.);
    setTopMassDelta(0.);
    setWMassDelta(0.);
}
*/
topSystemChiSquare::~topSystemChiSquare() {}

/*void topSystemChiSquare::setBJet(double px, double py, double pz, double e)
{
    bJetPx_ = px;
    bJetPy_ = py;
    bJetPz_ = pz;
    bJetE_ = e;
}*/

/*void topSystemChiSquare::getBJet(double &px, double &py, double &pz, double &e)
{
    //  double phiNew = reconstructed_bJetPhi_+bJetPhiDelta_*bJetPhiWidth_;
    //  double etaNew = reconstructed_bJetEta_+bJetEtaDelta_*bJetEtaWidth_;
    //  double ptNew  = reconstructed_bJetPt_*exp(bJetPtDelta_*bJetPtWidth_) ;
    //  double pxNew = ptNew*cos(phiNew);
    //  double pyNew = ptNew*sin(phiNew);
    //  double pzNew = ptNew*sinh(etaNew);
    //  double pNew  = ptNew*cosh(etaNew);
    //  double eNew  = sqrt(reconstructed_bJetMass2_ + pNew*pNew );
    //
    //  cout << "b-jet momentum calculated using deltas"
    //       << "\npx = " << pxNew
    //       << "\npy = " << pyNew
    //       << "\npz = " << pzNew
    //       << "\ne  = " << eNew  << endl;
    //
    //  cout <<   "bJetPx_ = " << bJetPx_
    //       << "\nbJetPy_ = " << bJetPy_
    //       << "\nbJetPz_ = " << bJetPz_
    //       << "\nbJetE_  = " << bJetE_  << endl;

    px = bJetPx_;
    py = bJetPy_;
    pz = bJetPz_;
    e = bJetE_;
}*/

/*void topSystemChiSquare::setObservedBJet(double px, double py, double pz,
                                         double e)
{
    reconstructed_bJetLorentzVector_ = XYZTLorentzVector(px, py, pz, e);
    reconstructed_bJetPt_ = reconstructed_bJetLorentzVector_.Pt();
    reconstructed_bJetPhi_ = reconstructed_bJetLorentzVector_.Phi();
    reconstructed_bJetEta_ = reconstructed_bJetLorentzVector_.Eta();
    reconstructed_bJetMass2_ = reconstructed_bJetLorentzVector_.M2();
}

void topSystemChiSquare::setBJetWidths(double sigmaPt, double sigmaPhi,
                                       double sigmaEta)
{
    bJetPtWidth_ = sigmaPt;
    bJetPhiWidth_ = sigmaPhi;
    bJetEtaWidth_ = sigmaEta;
}

void topSystemChiSquare::setWDaughter1(double px, double py, double pz,
                                       double e)
{
    WDaughter1Px_ = px;
    WDaughter1Py_ = py;
    WDaughter1Pz_ = pz;
    WDaughter1E_ = e;
}

void topSystemChiSquare::getWDaughter1(double &px, double &py, double &pz,
                                       double &e)
{
    //  double phiNew =
    //  reconstructed_WDaughter1Phi_+WDaughter1PhiDelta_*WDaughter1PhiWidth_;
    //  double etaNew =
    //  reconstructed_WDaughter1Eta_+WDaughter1EtaDelta_*WDaughter1EtaWidth_;
    //  double ptNew  =
    //  reconstructed_WDaughter1Pt_*exp(WDaughter1PtDelta_*WDaughter1PtWidth_) ;
    //  double pxNew = ptNew*cos(phiNew);
    //  double pyNew = ptNew*sin(phiNew);
    //  double pzNew = ptNew*sinh(etaNew);
    //  double pNew  = ptNew*cosh(etaNew);
    //  double eNew  = sqrt(reconstructed_WDaughter1Mass2_ + pNew*pNew );
    //
    //  cout << "first W daughter momentum calculated using deltas"
    //       << "\npx = " << pxNew
    //       << "\npy = " << pyNew
    //       << "\npz = " << pzNew
    //       << "\ne  = " << eNew  << endl;
    //
    //  cout <<   "WDaughter1Px_ = " << WDaughter1Px_
    //       << "\nWDaughter1Py_ = " << WDaughter1Py_
    //       << "\nWDaughter1Pz_ = " << WDaughter1Pz_
    //       << "\nWDaughter1E_  = " << WDaughter1E_  << endl;

    px = WDaughter1Px_;
    py = WDaughter1Py_;
    pz = WDaughter1Pz_;
    e = WDaughter1E_;
}

void topSystemChiSquare::setObservedWDaughter1(double px, double py, double pz,
                                               double e)
{
    reconstructed_WDaughter1LorentzVector_ = XYZTLorentzVector(px, py, pz, e);

    reconstructed_WDaughter1Pt_ = reconstructed_WDaughter1LorentzVector_.Pt();
    reconstructed_WDaughter1Phi_ = reconstructed_WDaughter1LorentzVector_.Phi();
    reconstructed_WDaughter1Eta_ = reconstructed_WDaughter1LorentzVector_.Eta();
    reconstructed_WDaughter1Mass2_ =
        reconstructed_WDaughter1LorentzVector_.M2();
}

void topSystemChiSquare::setWDaughter1Widths(double sigmaPt, double sigmaPhi,
                                             double sigmaEta)
{
    WDaughter1PtWidth_ = sigmaPt;
    WDaughter1PhiWidth_ = sigmaPhi;
    WDaughter1EtaWidth_ = sigmaEta;
}*/

// void topSystemChiSquare::resetAll()
//{
//  resetBJet();
//  resetWDaughter1();
//}

/*void topSystemChiSquare::resetBJet()
{
    double phiNew = reconstructed_bJetPhi_ + bJetPhiDelta_ * bJetPhiWidth_;
    double etaNew = reconstructed_bJetEta_ + bJetEtaDelta_ * bJetEtaWidth_;
    // double ptNew = reconstructed_bJetPt_ * exp(bJetPtDelta_ * bJetPtWidth_);
    double ptNew = reconstructed_bJetPt_ + bJetPtDelta_ * bJetPtWidth_;
    double pxNew = ptNew * cos(phiNew);
    double pyNew = ptNew * sin(phiNew);
    double pzNew = ptNew * sinh(etaNew);
    double pNew = ptNew * cosh(etaNew);
    double eNew = sqrt(reconstructed_bJetMass2_ + pNew * pNew);

    // reset values
    bJetPx_ = pxNew;
    bJetPy_ = pyNew;
    bJetPz_ = pzNew;
    bJetE_ = eNew;

    // cout << "bJetPx_ = " << bJetPx_ << endl ;
    // cout << "bJetPy_ = " << bJetPy_ << endl ;
    // cout << "bJetPz_ = " << bJetPz_ << endl ;
    // cout << "bJetE_  = " << bJetE_  << endl ;
}

void topSystemChiSquare::resetWDaughter1()
{
    double phiNew = reconstructed_WDaughter1Phi_ +
                    WDaughter1PhiDelta_ * WDaughter1PhiWidth_;
    double etaNew = reconstructed_WDaughter1Eta_ +
                    WDaughter1EtaDelta_ * WDaughter1EtaWidth_;
    // double ptNew = reconstructed_WDaughter1Pt_ *
    exp(WDaughter1PtDelta_ * WDaughter1PtWidth_);
    double ptNew =
        reconstructed_WDaughter1Pt_ + WDaughter1PtDelta_ * WDaughter1PtWidth_;
    double pxNew = ptNew * cos(phiNew);
    double pyNew = ptNew * sin(phiNew);
    double pzNew = ptNew * sinh(etaNew);
    double pNew = ptNew * cosh(etaNew);
    double eNew = sqrt(reconstructed_WDaughter1Mass2_ + pNew * pNew);

    // reset values
    WDaughter1Px_ = pxNew;
    WDaughter1Py_ = pyNew;
    WDaughter1Pz_ = pzNew;
    WDaughter1E_ = eNew;
}*/

/*void topSystemChiSquare::setBJetDeltas(const double ptDelta,
                                       const double phiDelta,
                                       const double etaDelta)
{
    //  cout << "Setting b-jet deltas to"
    //       << "\ndelta pt  = " << ptDelta
    //       << "\ndelta phi = " << phiDelta
    //       << "\ndelta eta = " << etaDelta
    //       << endl;
    bJetPtDelta_ = ptDelta;
    bJetPhiDelta_ = phiDelta;
    bJetEtaDelta_ = etaDelta;
    resetBJet();
}

void topSystemChiSquare::setWDaughter1Deltas(const double ptDelta,
                                             const double phiDelta,
                                             const double etaDelta)
{
    //  cout << "Setting first W daughter deltas to"
    //       << "\ndelta pt  = " << ptDelta
    //       << "\ndelta phi = " << phiDelta
    //       << "\ndelta eta = " << etaDelta
    //       << endl;
    WDaughter1PtDelta_ = ptDelta;
    WDaughter1PhiDelta_ = phiDelta;
    WDaughter1EtaDelta_ = etaDelta;
    resetWDaughter1();
}*/

/*void topSystemChiSquare::setDeltas(double bJetPtDelta, double bJetPhiDelta,
                                   double bJetEtaDelta,
                                   double WDaughter1PtDelta,
                                   double WDaughter1PhiDelta,
                                   double WDaughter1EtaDelta,
                                   // double topMassDelta,
                                   double WMassDelta)
{
    setBJetDeltas(bJetPtDelta, bJetPhiDelta, bJetEtaDelta);
    setWDaughter1Deltas(WDaughter1PtDelta, WDaughter1PhiDelta,
                        WDaughter1EtaDelta);
    // setTopMassDelta(topMassDelta);
    setWMassDelta(WMassDelta);
}*/

void topSystemChiSquare::setEllipseAngle()
{
//  theta_=theta;
   resetWDaughter2();
//  calcTopMomentum();
}

void topSystemChiSquare::resetWDaughter2()
{
    // cout << "Calculating the second W daughter momentum at angle: " << theta
    // << endl;
    WDaughter2Calc_.getWDaughterMomentum();
    //TVectorD *p = WDaughter2Calc_.getWDaughterMomentum(theta);
    //if (p->GetNrows() != 3)
    //    cout << "Error calculating the W daughter momentum: it does not have 3 "
    //            "components"
    //         << endl;
    //else {
    //    WDaughter2Px_ = (*p)[0];
    //    WDaughter2Py_ = (*p)[1];
    //    WDaughter2Pz_ = (*p)[2];
/*        WDaughter2Px_ = topsys.vars.Wd2_px;
        WDaughter2Py_ = topsys.vars.Wd2_py;
        WDaughter2Pz_ = topsys.vars.Wd2_pz;
        WDaughter2E_ =
            sqrt(mWDaughter2_ * mWDaughter2_ + WDaughter2Px_ * WDaughter2Px_ +
                 WDaughter2Py_ * WDaughter2Py_ + WDaughter2Pz_ * WDaughter2Pz_);

        WDaughter2Pt_ =
            sqrt(WDaughter2Px_ * WDaughter2Px_ + WDaughter2Py_ * WDaughter2Py_);
        WDaughter2Phi_ = atan2(WDaughter2Py_, WDaughter2Px_);
        WDaughter2Eta_ = asinh(WDaughter2Pz_ / WDaughter2Pt_);*/
    //}

    //  cout << "The second W daughter momentum is "
    //       << "\npx = " << WDaughter2Px_
    //       << "\npy = " << WDaughter2Py_
    //       << "\npz = " << WDaughter2Pz_
    //       << "\nE  = " << WDaughter2E_   << endl;
    //
    //  cout << "or "
    //       << "\npt  = " << WDaughter2Pt_
    //       << "\nphi = " << WDaughter2Phi_
    //       << "\neta = " << WDaughter2Eta_
    //       << endl;
}
/*
void topSystemChiSquare::getWDaughter2Momentum(double &px, double &py,
                                               double &pz)
{
    px = WDaughter2Px_;
    py = WDaughter2Py_;
    pz = WDaughter2Pz_;
    //  cout << "Second W daughter momentum at angle " << theta_ << " is "
    //       << "\npx = " << WDaughter2Px_
    //       << "\npy = " << WDaughter2Py_
    //       << "\npz = " << WDaughter2Pz_
    //       << "\nE  = " << WDaughter2E_
    //       << endl;
}

void topSystemChiSquare::getWDaughter2(double &px, double &py, double &pz,
                                       double &e)
{
    px = WDaughter2Px_;
    py = WDaughter2Py_;
    pz = WDaughter2Pz_;
    e = sqrt(max(0., mWDaughter2_ + WDaughter2Px_ * WDaughter2Px_ +
                         WDaughter2Py_ * WDaughter2Py_ +
                         WDaughter2Pz_ * WDaughter2Pz_));
    //  cout << "Second W daughter momentum is "
    //       << "\npx = " << px
    //       << "\npy = " << py
    //       << "\npz = " << pz
    //       << "\nE  = " << e
    //       << endl;
}
*/
void topSystemChiSquare::preSetupWDaughter2Ellipse()
{
    WDaughter2Calc_.preSetupEllipse();
}

void topSystemChiSquare::setupWDaughter2EllipsePart2()
{
    WDaughter2Calc_.setupEllipsePart2();
}

/*void topSystemChiSquare::setupWDaughter2Ellipse()
{
    // cout << "deltaMTop = " << deltaMTop_ << endl;
    // cout << "deltaMW   = " << deltaMW_   << endl;
    //WDaughter2Calc_.setupEllipse(mTop_ + sigmaMTop_ * deltaMTop_,
                                 //mW_ + sigmaMW_ * deltaMW_, mWDaughter2_);
    WDaughter2Calc_.setupEllipse();
}*/

void topSystemChiSquare::calcWDaughter2Ellipse()
{
    // cout << "Calculating topmassrange" <<endl;
    calcTopMassRange();
    // cout << "Calculating the ellipse in homogenous representation" << endl;
    WDaughter2Calc_.calcWDaughterEllipse();
    // cout << "Calculating the ellipse in extended representation" << endl;
    WDaughter2Calc_.calcExtendedWDaughterEllipse();
}

void topSystemChiSquare::calc_hperp_nperp()
{
    WDaughter2Calc_.calcWDaughterEllipse();
    // cout << "Calculating the ellipse in extended representation" << endl;
    WDaughter2Calc_.calcExtendedWDaughterEllipse();

}

//double topSystemChiSquare::getZ2(double mTop, double mW, double mWDaughter2)
double topSystemChiSquare::getZ2()
{
    //return WDaughter2Calc_.getZ2(mTop, mW, mWDaughter2);
    return WDaughter2Calc_.getZ2();
}

TMatrixD *topSystemChiSquare::getHomogeneousWDaughterEllipse()
{
    return WDaughter2Calc_.getHomogeneousWDaughterEllipse();
}

void topSystemChiSquare::calcTopMassRange()
{
    // Find the ranges of squared top mass values in which Z^2 is positive:
    // Z^2=0 is a polynomial equation of degree 2 in mTop^2

    // First check whether the calculation has been done
    //double currentZ2 = getZ2(topsys.input.mTop_central + topsys.input.mTop_width * topsys.vars.delta_mTop,
    //                         topsys.input.mW_central + topsys.input.mW_width * topsys.vars.delta_mW, topsys.input.Wd2_m);
    double currentZ2 = getZ2();

    if (rangeFlag_ && currentZ2 > 0.) {
        // cout << "Top mass range already calculated" << endl;
        return;
    }

    rangeFlag_ = false;

    // cout << "Using pole mass " << topsys.input.mTop_central << " GeV" << endl;

    double mW = topsys.input.mW_central + topsys.input.mW_width * topsys.vars.delta_mW;
    double mW2 = mW * mW;
    double WDaughter2Mass2 = topsys.input.Wd2_m * topsys.input.Wd2_m;

    double bJetE2 = topsys.calc.b_e() * topsys.calc.b_e();
    double bJetP2 = topsys.calc.b_px() * topsys.calc.b_px() + topsys.calc.b_py() * topsys.calc.b_py() + topsys.calc.b_pz() * topsys.calc.b_pz();
    double bJetP = sqrt(bJetP2);

    double WDaughter1P2 = topsys.calc.Wd1_px() * topsys.calc.Wd1_px() +
                          topsys.calc.Wd1_py() * topsys.calc.Wd1_py() +
                          topsys.calc.Wd1_pz() * topsys.calc.Wd1_pz();
    double WDaughter1P = sqrt(WDaughter1P2);
    double WDaughter1P4 = WDaughter1P2 * WDaughter1P2;
    double WDaughter1E2 = topsys.calc.Wd1_e() * topsys.calc.Wd1_e();
    double WDaughter1E3 = WDaughter1E2 * topsys.calc.Wd1_e();
    double WDaughter1E4 = WDaughter1E2 * WDaughter1E2;

    double c = (topsys.calc.Wd1_px() * topsys.calc.b_px() + topsys.calc.Wd1_py() * topsys.calc.b_py() +
                topsys.calc.Wd1_pz() * topsys.calc.b_pz());
    double s2 = 1. - c * c / (bJetP2 * WDaughter1P2);
    c /= (bJetP * WDaughter1P);

    double mTopEdge1 =
        (topsys.calc.b_e() * WDaughter1E3 -
         topsys.calc.b_e() * topsys.calc.Wd1_e() * (WDaughter2Mass2 - mW2 + WDaughter1P2) +
         WDaughter1E2 * (bJetE2 + mW2 - bJetP * (bJetP + c * WDaughter1P)) +
         WDaughter1P * ((-bJetE2 - mW2 + bJetP2) * WDaughter1P +
                        c * bJetP * (WDaughter2Mass2 - mW2 + WDaughter1P2)) -
         sqrt((WDaughter1E4 + pow(WDaughter2Mass2 - mW2, 2) +
               2 * (WDaughter2Mass2 + mW2) * WDaughter1P2 + WDaughter1P4 -
               2 * WDaughter1E2 * (WDaughter2Mass2 + mW2 + WDaughter1P2)) *
              (pow(c * topsys.calc.Wd1_e() * bJetP - topsys.calc.b_e() * WDaughter1P, 2) +
               bJetP2 * (topsys.calc.Wd1_e() - WDaughter1P) *
                   (topsys.calc.Wd1_e() + WDaughter1P) * s2))) /
        ((topsys.calc.Wd1_e() - WDaughter1P) * (topsys.calc.Wd1_e() + WDaughter1P));
    double mTopEdge2 =
        (topsys.calc.b_e() * WDaughter1E3 -
         topsys.calc.b_e() * topsys.calc.Wd1_e() * (WDaughter2Mass2 - mW2 + WDaughter1P2) +
         WDaughter1E2 * (bJetE2 + mW2 - bJetP * (bJetP + c * WDaughter1P)) +
         WDaughter1P * ((-bJetE2 - mW2 + bJetP2) * WDaughter1P +
                        c * bJetP * (WDaughter2Mass2 - mW2 + WDaughter1P2)) +
         sqrt((WDaughter1E4 + pow(WDaughter2Mass2 - mW2, 2) +
               2 * (WDaughter2Mass2 + mW2) * WDaughter1P2 + WDaughter1P4 -
               2 * WDaughter1E2 * (WDaughter2Mass2 + mW2 + WDaughter1P2)) *
              (pow(c * topsys.calc.Wd1_e() * bJetP - topsys.calc.b_e() * WDaughter1P, 2) +
               bJetP2 * (topsys.calc.Wd1_e() - WDaughter1P) *
                   (topsys.calc.Wd1_e() + WDaughter1P) * s2))) /
        ((topsys.calc.Wd1_e() - WDaughter1P) * (topsys.calc.Wd1_e() + WDaughter1P));

    // assert(mTopEdge1==mTopEdge1);
    // assert(mTopEdge2==mTopEdge2);

    double minRoot = min(mTopEdge1, mTopEdge2);
    double maxRoot = max(mTopEdge1, mTopEdge2);

    mTopEdge1 = minRoot;
    mTopEdge2 = maxRoot;

    //  if(mTopEdge1<=-DBL_MAX)
    if (mTopEdge1 <= 0) {
        cout << "infinite lower edge" << endl;
        // cout << "first W daughter energy is " << topsys.calc.Wd1_e() << endl;
        // cout << "first W daughter momentum is " << WDaughter1P << endl;
        // mTopEdge1 = 0;
        // return;
    }
    //  if(mTopEdge2>= DBL_MAX)
    if (mTopEdge2 >= 500) {
        cout << "infinite upper edge" << endl;
        // cout << "first W daughter energy is " << topsys.calc.Wd1_e() << endl;
        // cout << "first W daughter momentum is " << WDaughter1P << endl;
        // mTopEdge2 = 500;
        // return;
    }

    // double tempTheta=theta_;

    cout << "Printing out some starting values" << endl;

    // cout << "Current top mass is " << tempTopMass << endl;
    //cout << "Current ellipse angle is " << getEllipseAngle() << endl;

    double Z2, mTop;
    topsys.vars.has_high_edge = false;

    if (mTopEdge1 > 0 && mTopEdge2 > 0) {
        // both roots are positive -> 3 top mass intervals in which to check
        // sign of Z^2

        mTopEdge1 = sqrt(mTopEdge1);
        mTopEdge2 = sqrt(mTopEdge2);

        bool pos1 = false, pos2 = false, pos3 = false;

        mTop = 0.5 * mTopEdge1;
        topsys.vars.delta_mTop = mTop - topsys.input.mTop_central;
        Z2 = getZ2();
        if (Z2 > 0)
            pos1 = true;

        mTop = 0.5 * (mTopEdge1 + mTopEdge2);
        topsys.vars.delta_mTop = mTop - topsys.input.mTop_central;
        Z2 = getZ2();
        if (Z2 > 0)
            pos2 = true;

        mTop = 2. * mTopEdge2;
        topsys.vars.delta_mTop = mTop - topsys.input.mTop_central;
        Z2 = getZ2();
        if (Z2 > 0)
            pos3 = true;

        if (pos1 && pos2 && pos3) {
            cout << "Z^2 is positive on all 3 intervals" << endl;
        } else if (pos1 && pos2 && !pos3) {
            cout << "Z^2 is positive on 2 consecutive intervals below "
                 << mTopEdge2 << " and negative above" << endl;
        } else if (pos1 && !pos2 && pos3) {
            // cout << "Z^2 is positive below " << mTopEdge1 << " and above " <<
            // mTopEdge2 << endl;
            if (topsys.input.mTop_central < mTopEdge1) {
                rangeFlag_ = true;
                topsys.vars.mTop_edge_low = 0.;
                topsys.vars.mTop_edge_high = mTopEdge1;
                topsys.vars.has_high_edge = true;
            } else if (topsys.input.mTop_central > mTopEdge1 && topsys.input.mTop_central < mTopEdge2) {
                rangeFlag_ = true;
                if (topsys.input.mTop_central - mTopEdge1 < mTopEdge2 - topsys.input.mTop_central) {
                    topsys.vars.mTop_edge_low = 0.;
                    topsys.vars.mTop_edge_high = mTopEdge1;
                    topsys.vars.has_high_edge = true;
                } else {
                    topsys.vars.mTop_edge_low = mTopEdge2;
                    topsys.vars.mTop_edge_high = -1.;
                    topsys.vars.has_high_edge = false;
                }
            } else if (topsys.input.mTop_central > mTopEdge2) {
                rangeFlag_ = true;
                topsys.vars.mTop_edge_low = mTopEdge2;
                topsys.vars.mTop_edge_high = -1.;
                topsys.vars.has_high_edge = false;
            } else {
                cout << "Logic fail: check the roots" << endl;
            }
        } else if (pos1 && !pos2 && !pos3) {
            cout << "Z^2 is positive below " << mTopEdge1
                 << " and negative on 2 consecutive intervals above" << endl;
        } else if (!pos1 && pos2 && pos3) {
            cout << "Z^2 is negative below " << mTopEdge1
                 << " and positive on 2 consecutive intervals above" << endl;
        } else if (!pos1 && pos2 && !pos3) {
            // cout << "Z^2 is positive between " << mTopEdge1 << " and " <<
            // mTopEdge2 << endl;
            rangeFlag_ = true;
            topsys.vars.mTop_edge_low = mTopEdge1;
            topsys.vars.mTop_edge_high = mTopEdge2;
            topsys.vars.has_high_edge = true;
        } else if (!pos1 && !pos2 && pos3) {
            cout << "Z^2 is negative on 2 consecutive intervals below "
                 << mTopEdge2 << " and positive above" << endl;
        } else if (!pos1 && !pos2 && !pos3) {
            cout << "Z^2 is negative on all 3 intervals" << endl;
        } else {
            cout << "Logic fail" << endl;
        }

    } else if (mTopEdge2 > 0) {
        // only one positive root -> 2 intervals to check: between 0 and
        // sqrt(mTopEdge2), above sqrt(mTopEdge2)

        mTopEdge2 = sqrt(mTopEdge2);

        bool pos1 = false, pos2 = false;

        mTop = 0.5 * mTopEdge2;
        topsys.vars.delta_mTop = mTop - topsys.input.mTop_central;
        Z2 = getZ2();
        if (Z2 > 0)
            pos1 = true;

        mTop = 2. * mTopEdge2;
        topsys.vars.delta_mTop = mTop - topsys.input.mTop_central;
        Z2 = getZ2();
        if (Z2 > 0)
            pos2 = true;

        if (pos1 && pos2) {
            cout << "Z^2 is positive on both sides of the root: " << mTopEdge2
                 << endl;
        } else if (pos1) {
            rangeFlag_ = true;
            topsys.vars.mTop_edge_low = 0.;
            topsys.vars.mTop_edge_high = mTopEdge2;
            topsys.vars.has_high_edge = true;
        } else if (pos2) {
            rangeFlag_ = true;
            topsys.vars.mTop_edge_low = mTopEdge2;
            topsys.vars.mTop_edge_high = -1.;
            topsys.vars.has_high_edge = false;
        } else {
            cout << "Z^2 is negative on both sides of the root: " << mTopEdge2
                 << endl;
        }
    } else {
        // both roots are negative -> 1 interval to check: mTop>0
        mTop = 100;
        topsys.vars.delta_mTop = mTop - topsys.input.mTop_central;
        Z2 = getZ2();
        if (Z2 > 0) {
            rangeFlag_ = true;
            topsys.vars.mTop_edge_low = 0.;
            topsys.vars.mTop_edge_high = -1.;
            topsys.vars.has_high_edge = false;
        } else {
            cout << "Both roots are negative, and Z^2 is negative above 0"
                 << endl;
        }
    }

    cout << "before rangeflag check" << endl;
    if (!rangeFlag_) {
        cout << "No top mass range where Z^2 is positive was found" << endl;
        return;
    }
    cout << "aftercheck" << endl;

    if (topsys.vars.has_high_edge)
        cout << "Allowed top mass range is [ " << topsys.vars.mTop_edge_low << " , "
             << topsys.vars.mTop_edge_high << " ]" << endl;
    else
        cout << "Allowed top mass range is [ " << topsys.vars.mTop_edge_low << " , +inf ["
             << endl;

    if (topsys.vars.has_high_edge && topsys.vars.mTop_edge_high < topsys.vars.mTop_edge_low)
        cout << "Inverted interval" << endl;

    // Now calculate the top mass delta range

    if (topsys.vars.mTop_edge_low < mW) {
        topsys.vars.delta_mTop_range_low = (mW - topsys.input.mTop_central) / topsys.input.mTop_width;
        cout << "first " << topsys.vars.delta_mTop_range_low << endl;
    } else {
        topsys.vars.delta_mTop_range_low =
            (topsys.vars.mTop_edge_low - topsys.input.mTop_central) / topsys.input.mTop_width; // can be positive or negative
        cout << "mtopedgelow = " << topsys.vars.mTop_edge_low << endl;
        cout << "mtop = " << topsys.input.mTop_central << endl;
        cout << "sigmamtop = " << topsys.input.mTop_width << endl;
        cout << "mtopedge1 = " << mTopEdge1 << endl;
        cout << "mtopedge2 = " << mTopEdge2 << endl;
        cout << topsys.vars.delta_mTop_range_low << endl;
    }

    if (topsys.vars.has_high_edge) {
        topsys.vars.delta_mTop_range_high =
            (topsys.vars.mTop_edge_high - topsys.input.mTop_central) / topsys.input.mTop_width; // can be positive or negative
    } else {
        topsys.vars.delta_mTop_range_high = -1.;
    }

    cout << "Setting lower edge delta to " << topsys.vars.delta_mTop_range_low << endl;
    if (topsys.vars.has_high_edge)
        cout << "Setting upper edge delta to " << topsys.vars.delta_mTop_range_high << endl;
    //
    cout << "setting low edge to " << topsys.input.mTop_central + topsys.vars.delta_mTop_range_low * topsys.input.mTop_width
         << endl;
    if (topsys.vars.has_high_edge)
        cout << "setting high edge to "
             << topsys.input.mTop_central + topsys.vars.delta_mTop_range_high * topsys.input.mTop_width << endl;

    // reset to initial mass value

    if (currentZ2 <= 0) {
        cout << "Setting the top mass to a value that makes the current Z^2 "
                "positive"
             << endl;

        //SM's question: isn't the whole point of finding the range that Z2 should be positive within this range?
        //So if we set delta_mTop to be the low point of this range, the resulting Z2 should be positive:
        //shouldn't need to keep adding little bits to delta_mTop in order to achieve positive Z2.
        topsys.vars.delta_mTop = topsys.vars.delta_mTop_range_low;
        currentZ2 = getZ2(); //added by SM
        int whilecount = 0;
        while (currentZ2 <= 0 and whilecount <= 10000) {
            cout << "mtop is " << topsys.input.mTop_central << endl;
            topsys.vars.delta_mTop += 1.e-3 * topsys.input.mTop_central;
            cout << "New deltaMTop is " << topsys.vars.delta_mTop << endl;
            cout << "New top mass is " << topsys.input.mTop_central + topsys.input.mTop_width * topsys.vars.delta_mTop
                 << endl;
            //setupWDaughter2EllipsePart2();// already called in getZ2()
            // WDaughter2Calc_.setupEllipse(topsys.input.mTop_central+topsys.input.mTop_width*topsys.vars.delta_mTop,mW,topsys.input.Wd2_m);
            currentZ2 =
                getZ2();
            cout << "Current Z^2 is " << currentZ2 << endl;
            whilecount++;
        }
        // WDaughter2Calc_.calcWDaughterEllipse();
        // WDaughter2Calc_.calcExtendedWDaughterEllipse();
    }

    else {
        cout << "Resetting to starting values" << endl;

        topsys.vars.delta_mTop = 0.; //added by SM.  SM's note: not sure about this?
        setupWDaughter2EllipsePart2();
        // WDaughter2Calc_.setupEllipse(topsys.input.mTop_central+topsys.input.mTop_width*topsys.vars.delta_mTop,mW,topsys.input.Wd2_m);
        // calcWDaughter2Ellipse();
        // WDaughter2Calc_.calcWDaughterEllipse();
        // WDaughter2Calc_.calcExtendedWDaughterEllipse();
        cout << "Resetting second W daughter" << endl;
        // if(tempTheta!=0) resetWDaughter2(tempTheta);

        // printWDaughter2();

        // cout << "Re-printing values -- they should match the starting values
        // printed above" << endl;

        // cout << "Current top mass is " << WDaughter2Calc_.getTopMass() <<
        // endl;

        // if(tempTheta!=getEllipseAngle()) cout << "Issue resetting ellipse
        // angle to starting value" << endl;

        // cout << "Current ellipse angle is " << tempTheta << endl;
    }
}

/*void topSystemChiSquare::getTopMassRange()
{
    calcTopMassRange();
    calcWDaughter2Ellipse();
    //hasHighEdge = topsys.vars.has_high_edge;
    //mLow = topsys.vars.mTop_edge_low;
    //mHigh = topsys.vars.mTop_edge_high;
}*/

/*
void topSystemChiSquare::getTopMassDeltaRange(bool &hasHighEdge,
                                              double &deltaLow,
                                              double &deltaHigh)
{
    calcTopMassRange();
    calcWDaughter2Ellipse();
    hasHighEdge = topsys.vars.has_high_edge;
    deltaLow = deltaMTopRangeLow_;
    deltaHigh = deltaMTopRangeHigh_;
}*/
/*
void topSystemChiSquare::calcTopMomentum()
{
    // cout << "Calculating top momentum" << endl;
    // resetWDaughter2(theta_);

    //  cout << "Current second W daughter momentum: "
    //       << "\npx = " << WDaughter2Px_
    //       << "\npy = " << WDaughter2Py_
    //       << "\npz = " << WDaughter2Pz_
    //       << "\ne  = " << WDaughter2E_
    //       << endl;

    topPx_ = bJetPx_ + WDaughter1Px_ + WDaughter2Px_;
    topPy_ = bJetPy_ + WDaughter1Py_ + WDaughter2Py_;
    topPz_ = bJetPz_ + WDaughter1Pz_ + WDaughter2Pz_;
    topE_ = bJetE_ + WDaughter1E_ + WDaughter2E_;

    //  cout << "Top momentum is "
    //       << "\npx = " << topPx_
    //       << "\npy = " << topPy_
    //       << "\npz = " << topPz_
    //       << "\nE  = " << topE_   << endl;
}

void topSystemChiSquare::getTopMomentum(double &px, double &py, double &pz,
                                        double &e)
{
    calcTopMomentum();
    px = topPx_;
    py = topPy_;
    pz = topPz_;
    e = topE_;
}*/
/*
double topSystemChiSquare::breitWignerError(const double &mass, const double &width,
                                            const double &deltaMass)
{
    double scaledDeltaMass = deltaMass * width;
    double scaledDeltaMass2 = scaledDeltaMass * scaledDeltaMass;
    double Zscore = normal_quantile(
        0.31830988618379067154 *
                atan2(scaledDeltaMass2 + 2. * scaledDeltaMass * mass,
                      mass * width) +
            0.5,
        1.0);
    return Zscore * Zscore;
}

double topSystemChiSquare::getTopMassChiSquare()
{
    return breitWignerError(mTop_, sigmaMTop_, deltaMTop_);
    // return deltaMTop_*deltaMTop_;
}*/


void topSystemChiSquare::printTopConstituents()
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

    //  cout << "reconstructed second light quark: "
    //       << "\npx = " << reconstructed_topsys.calc.Wd2_px()
    //       << "\npy = " << reconstructed_WDaughter2Py_
    //       << "\npz = " << reconstructed_WDaughter2Pz_
    //       << "\ne  = " << reconstructed_WDaughter2E_
    //       << endl;
    //
    //  cout << "Second W daughter momentum:"
    //       << "\npx = " << topSystemChiSquare::topsys.calc.Wd2_px()
    //       << "\npy = " << topSystemChiSquare::WDaughter2Py_
    //       << "\npz = " << topSystemChiSquare::WDaughter2Pz_
    //       << endl;
    //  double low, high;
    //  getTopMassRange(low,high);
}

void topSystemChiSquare::printWDaughter2()
{
    cout << "Current second light quark: "
         << "\npt  = " << topsys.calc.Wd2_pt() << "\nphi = " << topsys.calc.Wd2_phi()
         << "\neta = " << topsys.calc.Wd2_eta() << endl;

    cout << "theta is " << topsys.vars.theta << endl;

    //  cout << "reconstructed second light quark: "
    //       << "\npt  = " << reconstructed_WDaughter2Pt_
    //       << "\nphi = " << reconstructed_WDaughter2Phi_
    //       << "\neta = " << reconstructed_WDaughter2Eta_
    //       << endl;
    //
    cout << "dif pt  = " << topsys.calc.Wd2_dif_pt()
         << "\ndif phi = " << topsys.calc.Wd2_dif_phi()
         << "\ndif eta = " << topsys.calc.Wd2_dif_eta() << endl;

    // cout << "second W daughter ellipse in homogeneous coordinates:" << endl;
    // TMatrixD* a=getHomogeneousWDaughterEllipse();
    // a->Print();
}


