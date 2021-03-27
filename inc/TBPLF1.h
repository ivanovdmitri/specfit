// Dmitri Ivanov <dmiivanov@gmail.com>

//
//  Class that describes the broken power law functions, includes arbitrary
//  number of break points.  Inherits from TSPECFITF1 class that
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

  TBPLF1() :
      fBplScaleFactor(1.0)
  {
    ;
  }

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
      TSPECFITF1(name, make_formula(nbreaks, ftype, scalefactor, log10en_min), log10en_min, log10en_max, parnames, params, parerrors), fBplScaleFactor(scalefactor)
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
      TSPECFITF1(name, make_formula(nbreaks, ftype, scalefactor, log10en_min), log10en_min, log10en_max, csparnames, params, parerrors), fBplScaleFactor(scalefactor)
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
      TSPECFITF1(name, make_formula(nbreaks, ftype, scalefactor, log10en_min), log10en_min, log10en_max, csparnames, csparams, csparerrors), fBplScaleFactor(scalefactor)
  {
    set_default_title(ftype);
  }

  // Translate all parameters of the current instance into a new BPL instance except the new
  // instance can be of a different type (E.G. E^3 x J, EJ, etc)

  inline TBPLF1* NewTBPLF1(const char *newname,  // name of the new (constructed) function
      const char *ftype                          // understands J, E3J, EJ, J>, E2J>
      ) const
  {
    return new TBPLF1(newname, GetNbreaks(), ftype, GetBplScaleFactor(), GetXmin(), GetXmax(), &GetParNames().front(), GetParameters(), GetParErrors());
  }

  // make flux formulas based on the function type and how many break points
  static TString make_formula(Int_t nbreaks, const char *ftype, Double_t scalefactor, Double_t log10en_min);

  // get the number of break points depending on how many parameters
  // there are in the function
  Int_t GetNbreaks() const
  {
    return (GetNpar() - 2) / 2;
  }

  Int_t GetBplScaleFactor() const
  {
    return fBplScaleFactor;
  }

  // To re-scale the function
  void Scale(Double_t c)
  {
    TSPECFITF1::Scale(c);
    fBplScaleFactor *= c;
  }

  // Multiply this function (of log10en) by another function f (of log10en) and numerically integrate with respect to linear energy dE
  // taking into account the breaks (scales to 1 at lowest energy, then scales back
  // the result back to what it should be in order to avoid numerical rounding issues)
  // If f is a null pointer, then no function multiplication is done, just integrates the current function.
  Double_t MultiplyAndIntegrate_dE(const TF1 *f, Double_t log10en_start, Double_t log10en_end, Double_t esprel = 9.9999999999999998E-13) const;

  // Numerically integrate this function (of log10en) with respect to linear energy dE
  // taking into account the breaks (scales to 1 at lowest energy, then scales back
  // the result back to what it should be in order to avoid numerical rounding issues)
  inline Double_t Integrate_dE(Double_t log10en_start, Double_t log10en_end, Double_t esprel = 9.9999999999999998E-13) const
  {
    return MultiplyAndIntegrate_dE(0,log10en_start,log10en_end,esprel);
  }

  // get the power coefficient that applies after the break point
  // this takes into account all the prior break points
  static TString get_pcf(Int_t ibreak, Int_t nbreaaks, Double_t log10en_min);

private:

  // title according to the type of the function
  void set_default_title(const char *ftype);

  // scaling factor of the BPL function
  Double_t fBplScaleFactor;

ClassDef(TBPLF1,1)
  ;

};

#endif
