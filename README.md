This repository contains a C++ implementation of a method to improve object resolutions in events containing two top quarks.

This package is based on original code from Gala Kaufman and Luke Winstrom.  In this repo I have heavily restructured the code to make it more efficient and less error-prone.  I have also added additional functionality (results analysis, plotting).

This repo is no longer updated -- as of July 2016, using new version at https://github.com/shtan/ttH_framework.


## Instructions to run

The code can be run from the command line:

./topReconstructionFromLHE


## Class and file descriptions


### WDaughterEllipseCalculator

A C++ implementation of the algorithms presented in [Betchart et al, arXiv:1305.1878](http://arxiv.org/abs/1305.1878) (in Python). 
The framework has been extended to consider either leptonic or hadronic top decays.  In practice this does not introduce any changes to the class.

#### Methods of interest

- preSetupEllipse() and setupEllipsePart2()

Set the assumed top, W, and second W daughter masses and recalculate the different variables defined in the Betchart paper.

- calcWDaughterEllipse()

Calculate the second W daughter ellipse in the transverse plane px, py (Hperp).

- calcExtendedWDaughterEllipse()

Calculate the second W daughter ellipse in the extended representation px, py, pz (Nperp).

- getExtendedWDaughterEllipse()

Return Hperp.

- getHomogeneousWDaughterEllipse()

Return Nperp.

- getWDaughterMomentum()

Given a point on the second W daughter ellipse, defined by an angle theta, calculate the corresponding momentum (px, py, pz).


### lightJetMinimumChiSquareSolver

In the Betchart paper, the addition of a MET constraint is equivalent to requiring that the two neutrino ellipses intersect. This can be achieved by varying within their resolutions all the objects present in the event that are not associated with the top quarks. At a hadron collider these would be light jets. 
Note: could in theory also apply to leptons, but since these are so well measured at the LHC we simply focus on jets. However the code can deal with non-zero lepton resolutions, see below the topEventMinimizer method setNonTopObjectCollections which are input to the lightJetMinimumChiSquareSolver class.
In the general case we want the non-top objects to exactly cancel the momentum imbalance from the decay of the top quarks.

We can write and minimize a chi square variable representing the distance between the measured and the corrected light jets which make the ellipses intersect. Transforming from polar to cartesian coordinates yields an analytic solution for the minimum by inverting a covariance matrix of the sum of all the object resolutions. The purpose of this class, given an input collection of objects and their resolutions, is to calculate the minimum and corresponding corrections to the objects.


#### Methods of interest

- Eval_covariance

Transform from polar (pT, phi, eta) to Cartesian widths (px, py, pz).  Once the light jet object collections have been set, calculate the global covariance matrix and invert it.

- setupEquations

Pass in the vectors of non-top object four-momenta and resolutions. Calls setCartesianWidths and calcSigmas.

- calcMin

Calculate the vector of deltas (corrections to the light jets) defined in Eq. C.7, i.e., calculate the B_i matrices and multiply by the displacement vector. From the vector of deltas, calculate the chi square. This method needs to be called only when the displacement vector changes.


### topSystemChiSquare

Parent top system chi square class. The two different decay channels, leptonic and hadronic, are derived classes with decay-specific members and methods; the parent class contains members and methods common to both channels.

Construct a chi square for one top quark and its decay products, which is contructed from corrections (deltas) to the b-jet, one measured W daughter, top mass, and W mass. The chi square represents the degree of likelihood that a top quark candidate four-momentum can be reconstructed from the input decay product momenta. 

If the decay is leptonic, the neutrino must be reconstructed in order to obtain the top momentum, using the Betchart method. If the decay is hadronic, the b-jet and one of the measured W daughters can be used to construct the ellipse that the second W daughter must lie on, also with the Betchart method. The possible second W daughter momenta can then be compared to the measured second W daughter. 

Each instance has a WDaughterEllipseCalculator instance. The b-jet and first W daughter momenta and top and W masses are expressed as functions of the measured or assumed values, the different object widths, and corrections to the various objects. The WDaughterEllipseCalculator is instantiated with these values.

One limitation of the Betchart paper is that solutions do not exist in the case where the variable Z^2, defined in Section 2.4.2, is negative. We introduce a method to get around this issue: Z^2 == 0 is equivalent to a quadratic equation in the top mass squared. We can therefore solve for the range(s) of top mass values where Z^2 is positive, and restrict corrections to the top mass to stay within the corresponding range.

#### Methods of interest

- preSetupWDaughter2Ellipse and setupWDaughter2EllipsePart2

Calls the preSetupEllipse and setupEllipsePart2 methods of the WDaughterEllipseCalculator class. Allows to update the ellipse with each new step in the minimization.

- calcWDaughter2Ellipse

Calculates the ellipse matrix in both representations (Hperp, homogeneous and Nperp, extended). Since the extended representation involves inverting the Hperp matrix, care must be taken to ensure that it is invertible. This is equivalent to requiring Z^2>0, and so the method calcTopMassRange (see below) is first called.

- calcTopMassRange

Calculates the roots mTopEdge1 and mTopEdge 2 of a quadratic equation in mTop^2, then finds the range in which Z^2 is positive. Three cases must be considered:
	1. Both roots are negative: only need to check whether Z^2 is positive in the range mTop>0.
	2. One positive and one negative root: two intervals to check, [0; sqrt(mTopEdge2)] and [sqrt(mTopEdge2); +inf[.
	3. Both roots are positive: three intervals to check, [0; sqrt(mTopEdge1)], [sqrt(mTopEdge1); sqrt(mTopEdge2)] and [sqrt(mTopEdge2);+inf[. 
In case more than one top mass interval in which Z^2>0 exists, pick the one either containing or closest to the input assumed top mass. 
Also calculate the corresponding delta range in which the top mass will be allowed to vary during the minimization.
Update: this method is no longer called.  It seems to give very strange top mass ranges.


### topEventMinimizer

Functions at the event level to minimize the per-event chi square defined in Eq. 1, which is constructed from the sum of top system chi squares (one per top assumed to be present in the event, each represented by one topSystemChiSquare object) and the non-top object chi square (combining all measured objects in the event not assumed to originate from a top decay, and represented by one lightJetMinimiumChiSquareSolver object).

As can be seen in Eq. 3, the minimum chi square can be expressed as a series of nested minimizations. We implement an outer and an inner minimizer.

Although in practice the code has been tested using two assumed tops, it can in theory accommodate any number of tops, decaying in either mode. This is achieved with a vector of pairs of (a pointer to) a topSystemChiSquare object, and a boolean holding the mode of decay (true for leptonic, false for hadronic). The user then adds the desired number of tops to the event using the methods addLeptonicTop and addHadronicTop (see below).

#### Methods of interest

- findStartingValues

Find starting values for the ellipse angles by looping around ellipses. At each point the sum of top momenta is recalculated and the non-top objects are set to recoil against it; the sum of the non-top (light jet) and hadronic chi squares is then computed. The values corresponding to the minimum are then used as a starting point for the inner minimizer.

- minimizeNonTopChiSquare

The inner minimizer parameters are the ellipse angles and top mass deltas (one each per top system). For each top system the top mass delta parameter is restricted in range by the edges computed in the topSystemChiSquare::calcTopMassRange method (see above).  The chi square to minimize is the sum of top mass, hadronic second W daughter, and non-top object chi squares.

- minimizeTotalChiSquare

The outer minimizer parameters are the b-jet pT, eta and phi deltas, the first W daughter pT, eta and phi deltas, and the W mass deltas (one each per top system). At each step in the minimization procedure, the second W daughter ellipses are recalculated using the updated b-jet and first W daughter momenta and W mass. Next the inner piece is calculated, ie, the inner minimization performed. The total chi square is defined as the sum of the minimum of the inner chi square at this outer step, and the sum of top system chi squares.


### topReconstructionFromLHE

Class that tests the above code on a simulated ttbar sample.  Outputs original (after smearing) and reconstructed particle momenta, comparing each to the generator-level momenta.  Also outputs resolution plots comparing smeared and reconstructed momenta.  Also outputs chi-squared values for the various objects, as well as the total value.

Currently set up to do semi-leptonically decaying ttbar events.

