#ifndef LEPTONICTOPSYSTEMCHISQUARE
#define LEPTONICTOPSYSTEMCHISQUARE

#include "topSystemChiSquare.h"

class leptonicTopSystemChiSquare : public topSystemChiSquare
{
  private:
    //double chi2_;

  public:
    leptonicTopSystemChiSquare( commonstruct::top_system & );

    /*leptonicTopSystemChiSquare(const double &, const double &, const double &,
                               const double &, const double &, const double &,
                               const double &, const double &, const double &,
                               const double &, const double &, const double &,
                               const double &, const double &, const double &,
                               const double &, const double &, const double &);*/

    ~leptonicTopSystemChiSquare();

    virtual void printTopConstituents();

    virtual void setEllipseAngle();

/*    virtual void setWDaughter2(double, double, double) { return; };
    virtual void calcWDaughter2Deltas() { return; };
    virtual void getWDaughter2Deltas(double &, double &, double &);
    virtual void printWDaughter2() { return; };

    void calcChiSquare();
    double getChiSquare();
    double getWDaughter1ChiSquare();
    double getWMassChiSquare();
    double getBChiSquare();

    virtual double getHadronicChiSquare() { return 0; };*/
};

#endif
