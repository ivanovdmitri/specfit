// Dmitri Ivanov <dmiivanov@gmail.com>

//
// This class customizes the ROOT's TF1 function class so that it's easier
// to use when fitting cosmic ray flux results.

#ifndef _TSPECFITF1_h_
#define _TSPECFITF1_h_

#include <iostream>
#include <cstdlib>
#include "TF1.h"
#include "TString.h"
#include <vector>

// customized TF1 class to make the base class easier to use in the spectrum fit process
class TSPECFITF1: public TF1
{
public:

  TSPECFITF1()
  {
    ;
  }

  TSPECFITF1(const char *name,     // name of the class
      const char *frm,             // formula expression
      Double_t log10en_min,        // lowest energy  (log10(E/eV))
      Double_t log10en_max,        // highest energy (log10(E/eV))
      const TString *parnames,     // parameter names as an array of TString type objects
      const Double_t *params,      // values (starting values) of the parameters
      const Double_t *parerrors    // errors (starting step sizes) of the parameters
      ) :
      TF1(name, frm, log10en_min, log10en_max)
  {
    if(parnames)
      SetParNames(parnames);
    if(params)
      SetParameters(params);
    if(parerrors)
      SetParErrors(parerrors);
  }

  TSPECFITF1(const char *name,      //
      const char *frm,              //
      Double_t log10en_min = 18.0,  //
      Double_t log10en_max = 21.0,  //
      const char *csparnames = 0,   // comma-separated parameter names as a simple C-string
      const Double_t *params = 0,   //
      const Double_t *parerrors = 0 //
      ) :
      TF1(name, frm, log10en_min, log10en_max)
  {
    if(csparnames)
      SetParNamesCS(csparnames);
    if(params)
      SetParameters(params);
    if(parerrors)
      SetParErrors(parerrors);
  }
  TSPECFITF1(const char *name,     //
      const char *frm,             //
      Double_t log10en_min,        //
      Double_t log10en_max,        //
      const char *csparnames,      //
      const char *csparams,        // comma - separated list of values (starting values) of the parameters as a single C string
      const char *csparerrors      // comma - separated list of errors (starting step sizes) of the parameters as a single C string
      ) :
      TF1(name, frm, log10en_min, log10en_max)
  {
    if(csparnames)
      SetParNamesCS(csparnames);
    if(csparams)
      SetParametersCS(csparams);
    if(csparerrors)
      SetParErrorsCS(csparerrors);
  }

  virtual ~TSPECFITF1();

  TString GetParNamesCS()
  {
    return cs_strings(GetNpar(), GetParNames());
  }

  TString GetParametersCS()
  {
    return cs_doubles(GetNpar(), GetParameters());
  }

  TString GetParErrorsCS()
  {
    return cs_doubles(GetNpar(), GetParErrors());
  }

  static TString cs_doubles(Int_t n, const Double_t *data);
  static TString cs_strings(Int_t n, const TString *data);
  // parse comma-separated array of doubles into a vector of doubles
  static std::vector<Double_t> parse_cs_doubles(const char *cs_doubles);
  // parse comma-separated strings into a vector of strings
  static std::vector<TString> parse_cs_strings(const char *cs_strings);
  
  // in addition, add a method to set parameters from the array;
  // a method with these arguments is missing
  // in the TF1 class.
  void SetParNames(const TString *parnames);

  // set names, values, and errors from a string of comma-separated values
  void SetParNamesCS(const char *csparnames);
  void SetParametersCS(const char *csparameters);
  void SetParErrorsCS(const char *csparerrors);

  // get parameter names into array of TString objects
  void GetParNames(TString *parnames);

  // get a pointer to an array of TString objects that contains the parameter names
  const TString* GetParNames();

  // combining two functions with a starting formula frm_start and coefficients c1, c2 and
  // properly combines the parameter names, values, and errors.
  static TSPECFITF1* Add(const char *newname, const char *frm_start, const TF1 *f1, const TF1 *f2, Double_t c1 = 1.0, Double_t c2 = 1.0);

  // adding function to a starting formula
  static TSPECFITF1* Add(const char *newname, const char *frm_start, const TF1 *f, Double_t c = 1.0);

  // this essentially returns a scaled copy of a function
  static TSPECFITF1* Add(const char *newname, const TF1 *f, Double_t c = 1.0)
  {
    const char *frm_start = "0.0";
    return Add(newname, frm_start, f, c);
  }

  // this copies a given function
  static TSPECFITF1* MakeCopy(const char *newname, const TF1 *f)
  {
    return Add(newname, f, 1.0);
  }

  // so that it's easier to call this method from Python
  TString GetExpFormula(Option_t *opt = "")
  {
    return TF1::GetExpFormula(opt);
  }
  
  // Get the formula and apply an offset
  TString GetExpFormula(Int_t n_offset)
  {
    return GetExpFormula(this, n_offset);
  }

  // offset the parameters in the fomrula of the function by some integer value n_offset
  static TString GetExpFormula(const TF1 *f, Int_t n_offset);

private:

  // internal storage of the parameter names
  std::vector<TString> _parnames_data_;

ClassDef(TSPECFITF1,1)
  ;

};

#endif
