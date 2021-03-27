// Dmitri Ivanov <dmiivanov@gmail.com>

// class that performs the fit to multiple cosmic ray flux results with
// appropriate energy scale corrections for each. This classes uses
// TMinuit MINUIT minimizer to find the best fit parameters that minimize
// the overall binned Poisson log likelihood.

#ifndef _TCRFluxFit_h_
#define _TCRFluxFit_h_

#include <vector>
#include <map>
#include "TString.h"
#include "TCollection.h"
#include "TObjArray.h"
#include "TMinuit.h"
#include "TObject.h"
#include "TCRFlux.h"
#include "TF1.h"
#include "specfit_uti.h"

class TCRFluxFit: public TObject
{
public:
  TCRFluxFit() :
      log10en_min(17.0), log10en_max(21.0), nfitpar(0), nfluxpar(0), nencorrpar(0), chi2(0), ndof(0), fJ(0), fE3J(0), fJ_null(0), fE3J_null(0), mFIT(0)
  {
    ;
  }

  virtual ~TCRFluxFit();

  void SetEminEmax(Double_t min_log10en = 18.8, Double_t max_log10en = 21.0)
  {
    log10en_min = min_log10en;
    log10en_max = max_log10en;
  }

  // Set the pointer to the functions that describes the fitted flux versus log10(E/eV)
  // E^3 J function is optional, it's mainly used for plotting the results
  void SetFluxFun(TF1 *fJ_set, TF1 *fE3J_set = 0);

  // Set the pointer to the functions that describes the null hypothesis flux versus log10(E/eV)
  // E^3 J function is optional, it's mainly used for plotting the results
  void SetNullFun(TF1 *fJ_null_set, TF1 *fE3J_null_set = 0);

  // add the measured flux result plus the energy correction function that corresponds to that flux
  // this function assumes that the flux has been initialized elsewhere and the class will not attempt
  // to clean it up
  Bool_t Add(TCRFlux *flux, TF1 *fEnCorr_set = 0);

  // nbins: number of bins
  // log10en_values: energies log10(E/eV) of the bin centers
  // log10en_bsize_values: log10(E/eV) bin sizes
  // nevents_values: numbers of events in each energy bin
  // exposure_values: exposure [m^2 sr s] for each energy bin center value
  Bool_t Add(const char *name, const char *title, Int_t nbins, const Double_t *log10en_values, const Double_t *log10en_bsize_values, const Double_t *nevents_values,
      const Double_t *exposure_values, TF1 *fEnCorr_set = 0);

  // add a flux result from an ASCII file that has:
  // col1: energies log10(E/eV) of the bin centers
  // col2: log10(E/eV) bin sizes
  // col3: numbers of events in each energy bin
  // col4: exposure [m^2 sr s] for each energy bin center value
  Bool_t Add(const char *name, const char *title, const char *ascii_file, TF1 *fEnCorr_set = 0);

  // this selects the desired energy range for all fluxes
  void SelectEnergyRange(Double_t log10en_min = 17.0, Double_t log10en_max = 21.0);

  // return stored flux result
  TCRFlux* GetFlux(const char *name);

  // return number of stored flux results
  Int_t GetNfluxes() const
  {
    return (Int_t) Fluxes_ordered.size();
  }

  // return stored flux result by index
  TCRFlux* GetFlux(Int_t iflux);

  // Set the parameters for the flux functions and optionally errors on parameters
  // (errors may be useful for displaying the results)
  void SetFluxPar(const Double_t *params, const Double_t *parerrors = 0);

  // All correction functions must use the same parameters (but use them in different ways).
  // Set the parameters for the correction functions and optionally errors on parameters
  // (errors may be useful for displaying the results)
  void SetEncorrPar(const Double_t *params, const Double_t *parerrors = 0);

  // To set the parameters while fitting.
  // This depends on whether just the flux or the flux and
  // the nonlinearity correction are being fitted.  Parameters for the nonlinearity
  // correction follow those of the flux in the "par" array
  // that's being adjusted by Minuit minimizer while finding the minimum.
  void SetParameters(const Double_t *par);

  // calculate the overall log likelihood
  void CalcLogLikelihood();

  // calculate the numbers of events within a certain energy range that are expected from some null hypothesis flux function
  // and that are actually observed in the data
  std::pair<Double_t, Double_t> EvalNull();

  // calculate the overall log likelihood and return the (log likelihood, number of bins) pair
  const std::pair<Double_t, Double_t>& GetLogLikelihood()
  {
    CalcLogLikelihood();
    return log_likelihood;
  }

  // Performs the fit, returns true if successful.
  Bool_t Fit(Bool_t verbose = true);

  // To obtain the Minuit pointer for whatever reason.
  // In order for it to behave correctly, the global FCN must be
  // pointed to use this instance of the class.
  TMinuit* GetMinuit();

  Double_t log10en_min; // minimum log10(E/eV) for fitting
  Double_t log10en_max; // maximum log10(E/eV) for fitting
  Int_t nfitpar;        // total number of fit parameters
  Int_t nfluxpar;       // number of flux fit parameters
  Int_t nencorrpar;     // number of energy correction parameters
  Double_t chi2;        // normalized log likelihood, which in the case of large statistics behaves like chi2
  Double_t ndof;        // number of degrees of freedom
  std::map<TString, TCRFlux*> Fluxes;
  std::vector<TCRFlux*> Fluxes_ordered;
  std::pair<Double_t, Double_t> log_likelihood;
  std::pair<Double_t, Double_t> log_likelihood_nonzero;
  std::pair<Double_t, Double_t> log_likelihood_restricted;

  std::vector<Double_t> fit_parameters; // combined fit parameters
  std::vector<Double_t> fit_parerrors;  // uncertainties on combined fit parameters


  Double_t GetParameter(Int_t ipar) const
  {
    return ipar >=0 && ipar < (Int_t) fit_parameters.size() ? fit_parameters[ipar] : 0.0;
  }

  Double_t GetParError(Int_t ipar) const
  {
    return ipar >=0 && ipar < (Int_t) fit_parerrors.size() ? fit_parerrors[ipar] : 0.0;
  }

  Double_t GetChisquare() const
  {
    return log_likelihood.first;
  }

  Double_t GetNDF() const
  {
    return log_likelihood.second;
  }

  Int_t GetNpar() const
  {
    return (Int_t)fit_parameters.size();
  }

  // ipar is the parameter index
  // npts, par_lo, par_up are the number of points to consider and upper and lower limits.
  // Default values will lead to using Minuit's default settings
  // calc_deltas will calculate the offsets for the parameter from its best fitted value and
  // chi2 (normalized log likelihood) from its smallest value
  TGraph* scan_parameter(Int_t ipar, Int_t npts = 40, Double_t par_lo = 0, Double_t par_up = 0, Bool_t calc_deltas = true);

  // flux function to use
  TF1 *fJ;        // for fitting
  TF1 *fE3J;      // if supplied by the user

  TF1 *fJ_null;    // for null hypothesis
  TF1 *fE3J_null;  // if supplied by user

  // energy correction functions for each experiment (mapped by the experiment name),
  // if needed.  All these correction functions must use the same parameters, albeit in
  // different ways, if necessary.
  std::map<TString, TF1*> fEnCorr;

private:

  // minimizer
  TMinuit *mFIT; // the pointer can be obtained via special method by the outside code

  // collector for TCRFlux objects that have been internally created during the lifetime of the class
  TObjArray TCRFlux_Objects_Created_By_This;

ClassDef(TCRFluxFit,1)
  ;

};

#endif
