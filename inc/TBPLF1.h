// Dmitri Ivanov <dmiivanov@gmail.com>

//
//  Class that describes the broken power law functions, includes arbitrary
//  number of break points.  Inhertis from TSPECFITF1 class that
//  customizes ROOT's TF1 for the purposes of fitting cosmic ray flux results.
// 
 
#ifndef _TBPLF1_h_
#define _TBPLF1_h_

#include <iostream>
#include <cstdlib>
#include "TString.h"
#include <vector>
#include "TSPECFITF1.h"

class TBPLF1: public TSPECFITF1
{
public:

  TBPLF1() { ; }

  TBPLF1(const char *name,         // name of the class
      Int_t nbreaks,               // number of break points
      const char *ftype,           // understands J, E3J, EJ, J>, E2J>
      Double_t scalefactor,        // scaling factor at the lowest energy
      Double_t log10en_min,        // lowest energy  (log10(E/eV))
      Double_t log10en_max,        // highest energy (log10(E/eV))
      const TString *parnames,     // parameter names as an array of TString type objects
      const Double_t *params,      // values (starting values) of the parameters
      const Double_t *parerrors    // errors (starting step sizes) of the parameters
      ) :
      TSPECFITF1(name, make_formula(nbreaks, ftype, scalefactor, log10en_min), log10en_min, log10en_max, parnames, params, parerrors)
  {
    set_default_title(ftype);
  }

  TBPLF1(const char *name,          //
      Int_t nbreaks,                //
      const char *ftype = "J",      //
      Double_t scalefactor = 1.0,   //
      Double_t log10en_min = 18.0,  //
      Double_t log10en_max = 21.0,  //
      const char *csparnames = 0,   // comma-separated parameter names as a simple C-string
      const Double_t *params = 0,   //
      const Double_t *parerrors = 0 //

      ) :
      TSPECFITF1(name, make_formula(nbreaks, ftype, scalefactor, log10en_min), log10en_min, log10en_max, csparnames, params, parerrors)
  {
    set_default_title(ftype);
  }
  TBPLF1(const char *name,         //
      Int_t nbreaks,               //
      const char *ftype,           //
      Double_t scalefactor,        //
      Double_t log10en_min,        //
      Double_t log10en_max,        //
      const char *csparnames,      //
      const char *csparams,        // comma - separated list of values (starting values) of the parameters as a single C string
      const char *csparerrors      // comma - separated list of errors (starting step sizes) of the parameters as a single C string
      ) :
      TSPECFITF1(name, make_formula(nbreaks, ftype, scalefactor, log10en_min), log10en_min, log10en_max, csparnames, csparams, csparerrors)
  {
    set_default_title(ftype);
  }

  // make flux formulas based on the function type and how many break points
  static TString make_formula(Int_t nbreaks, const char *ftype, Double_t scalefactor, Double_t log10en_min);

  // get the number of break points depending on how many parameters
  // there are in the function
  Int_t GetNbreaks()
  {
    return (GetNpar() - 2) / 2;
  }

  // get the power coefficient that applies after the break point
  // this takes into account all the prior break points
  static TString get_pcf(Int_t ibreak, Int_t nbreaaks, Double_t log10en_min);

private:

  // title according to the type of the function
  void set_default_title(const char *ftype);

ClassDef(TBPLF1,1);

};

#endif
