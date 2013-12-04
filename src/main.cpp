
/** 
   Project:     SOFTSUSY 
   File:        main.cpp
   Author:      Ben Allanach 
   Manual:      B.C. Allanach,hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
                B.C. Allanach, M.A. Bernhardt, arXiv:0903.1805, Comp. Phys. 
		Commun. 181 (2010) 232-245
   Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html
   Description: main calling program example: performs a scan of tan beta 
   (starting at CMSSM10.1.1) and prints out Higgs masses as a result in
   the format
       <tan beta>      <mh0>     <mA0>    <mH0>     <mH+/->
*/

#include <iostream>
#include "mycomplex.h"
#include "def.h"
#include "linalg.h"
#include "lowe.h"
#include "rge.h"
#include "softsusy.h"
#include "softpars.h"
#include "susy.h"
#include "utils.h"
#include "numerics.h"

int main() {
  /// Sets up exception handling
  signal(SIGFPE, FPE_ExceptionHandler); 

  try {
 /// Sets format of output: 6 decimal places
  outputCharacteristics(6);

  cerr << "SOFTSUSY" << SOFTSUSY_VERSION 
       << " test program, Ben Allanach 2002\n";
  cerr << "If you use SOFTSUSY, please refer to B.C. Allanach,\n";
  cerr << "Comput. Phys. Commun. 143 (2002) 305, hep-ph/0104145\n";

  /// Parameters used: CMSSM parameters
  double m12 = 500., a0 = -1000., mGutGuess = 2.0e16, tanb = 10.0, m0 = 125.;
  int sgnMu = 1;      ///< sign of mu parameter 
  int numPoints = 200; ///< number of scan points

  QedQcd oneset;      ///< See "lowe.h" for default definitions parameters

  /// most important Standard Model inputs: you may change these and recompile
  double alphasMZ = 0.1187, mtop = 173.5, mbmb = 4.18;
  oneset.setAlpha(ALPHAS, alphasMZ);
  oneset.setPoleMt(mtop);
  oneset.setMbMb(mbmb);

  oneset.toMz();      ///< Runs SM fermion masses to MZ

  /// Print out the SM data being used, as well as quark mixing assumption and
  /// the numerical accuracy of the solution
  cout << "# Low energy data in SOFTSUSY: MIXING=" << MIXING << " TOLERANCE=" 
       << TOLERANCE << endl << oneset << endl;

  /// Print out header line
  cout << "# m0     m12     tanb      mh     Del(m0)   Del(m12)    Del(A)    Del(mu)   Del(m3sq)   Del(yt)   DelE   \n";

  int i; 
  /// Set limits of tan beta scan
  // double startTanb = 3.0, endTanb = 50.0;
  /// Cycle through different points in the scan
  tanb = 3;
  double startm0 = 10, endm0 = 5e+3;
  double startM12 = 10, endM12 = 5e+3;

  double mx=0.0;
  double mz2 =0.0, mu2 =0.0, Blow =0.0, B0 = 0.0;  
  //Individual sensitivities
  DoubleVector Del(6), Del2(5);
  //max sensitivities with and without Yt
  double DelBG, DelBG_noYt; 
  //sum in quadrature
  double DelE2 =0.0 , DelE=0.0, DelJ =0.0;
  double stepm0 = (endm0 - startm0) / double(numPoints);
  double stepm12 = (endM12 - startM12) / double(numPoints);
     
     for (i = 1; i<=numPoints; i++) {
        for (int j = 1; j<=numPoints; j++) {
           // tanb = (endTanb - startTanb) / double(numPoints) * double(i) +
           //   startTanb; // set tan beta ready for the scan.
           
           m0 = startm0 - stepm0 + stepm0 * double(i);
           m12 = startM12 - stepm12 + stepm12 * double(j);

           
           
           
           /// Preparation for calculation: set up object and input parameters
           MssmSoftsusy r; 
           DoubleVector pars(3); 
           pars(1) = m0; pars(2) = m12; pars(3) = a0;
           bool uni = true; // MGUT defined by g1(MGUT)=g2(MGUT)
           
           /// Calculate the spectrum
           mx = r.lowOrg(sugraBcs, mGutGuess, pars, sgnMu, tanb, oneset, uni);
           mz2 = sqr(r.displayMz());
           mu2 = sqr(r.displaySusyMu());
           Blow = r.displayM3Squared() / r.displaySusyMu();
           if (!r.displayProblem().test())  {   
              Del = r.fineTune(sugraBcs, pars, mx, true, B0);
              DelBG = Del.max();
              for(int k=1; k<=6;k++){
                 DelE2 = DelE2 + sqr(Del(k));
              }
              DelE = sqrt(DelE2);
              // cout << "DelBG with mt:           " << DelBG << endl;
              //keep for test then we can remove and simply set Del(6) = 0 before talking max to get DelBG_noYt.
              Del2 = r.fineTune(sugraBcs, pars, mx, false, B0);
              //cout << "DelBG without mt:        " << DelBG << endl;
              DelBG_noYt = Del2.max();
              

              //abs( (-mz2/mu2/2.0) * (1-tb2)/(1+tb2)*B/B0 * 1.0 
              DelJ = fabs( -mz2 / ( 2.0 * mu2 ) * (1-sqr(tanb))/(1+sqr(tanb)) * Blow / B0);
              
              // DelJ = abs( (-mz2/mu2/2.0) * (1-sqr(tanb))/(1+sqr(tanb))*B/B0 * 1.0 );
              cout << m0 << "   " << m12 << "   "  << tanb << "   " << r.displayPhys().mh0(1) << "   "  <<  Del(1) << "   " << Del2(1) << "   " <<  Del(2) << "   " << Del2(2) << "   " << Del(3) << "   " << Del2(3) << "   " << Del(4) << "   " << Del2(4) << "   " << Del(5) << "   " << Del2(5) << "   " << Del(6) <<  "   "<< DelBG <<  "   "<< DelBG_noYt << "    " << DelE  << "    "  << DelJ << endl;
           }

    /// check the point in question is problem free: if so print the output
   //  if (!r.displayProblem().test()) 
  //     cout << tanb << " " << r.displayPhys().mh0(1) << " " 
  //          << r.displayPhys().mA0(1) << " " 
  //          << r.displayPhys().mh0(2) << " " 
  //          << r.displayPhys().mHpm << endl;
  //   else
  //     /// print out what the problem(s) is(are)
  //     cout << tanb << " " << r.displayProblem() << endl;
  // } 
        }
     }

  }
  catch(const string & a) { cout << a; return -1; }
  catch(const char * a) { cout << a; return -1; }
  catch(...) { cout << "Unknown type of exception caught.\n"; return -1; }

  return 0;
}

