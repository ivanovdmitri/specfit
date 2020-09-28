// Dmitri Ivanov <dmiivanov@gmail.com>

// Class that describes cosmic ray flux measurements

#ifndef _TCRFlux_h_
#define _TCRFlux_h_

#include <vector>
#include <algorithm>
#include "TNamed.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"

// class that describes all the information + manipulation
// for the flux
class TCRFlux: public TNamed
{
public:

  TCRFlux() :
      log10en_min_data(0), log10en_max_data(0), nevents_min_restricted(7), fJ(0), fE3J(0), fJ_null(0), fE3J_null(0), fEnCorr(0)
  {
    ;
  }

  TCRFlux(const char *name, const char *title = "Flux");

  virtual ~TCRFlux();

  // Load data from an ASCII file with columns
  // col1: energies log10(E/eV) of the bin centers
  // col2: log10(E/eV) bin sizes
  // col3: numbers of events in each energy bin
  // col4: exposure [m^2 sr s] for each energy bin center value
  Bool_t Load(const char *ascii_file);

  // load data from C-like arrays
  Bool_t Load(Int_t nbins,                   // number of bins for the spectrum
      const Double_t *log10en_values,        // energies log10(E/eV) of the bin centers
      const Double_t *log10en_bsize_values,  // log10(E/eV) bin sizes
      const Double_t *nevents_values,        // numbers of events in each energy bin
      const Double_t *exposure_values);      // exposure [m^2 sr s] for each energy bin center value

  // load data from STL vectors
  Bool_t Load(const std::vector<Double_t> &log10en_values,  // energies log10(E/eV) of the bin centers
      const std::vector<Double_t> &log10en_bsize_values,    // log10(E/eV) bin sizes
      const std::vector<Double_t> &nevents_values,          // numbers of events in each energy bin
      const std::vector<Double_t> &exposure_values);        // exposure [m^2 sr s] for each energy bin center value

  // determine the energy range of the spectrum measurement
  void find_min_max_log10en();

  // this selects the data only within the desirable energy range
  void SelectEnergyRange(Double_t log10en_min = 17.0, Double_t log10en_max = 21.0);

  // Set the pointer to the functions that describes the fitted flux versus log10(E/eV)
  // E^3 J function is optional, it's mainly used for plotting the results
  void SetFluxFun(TF1 *fJ_set, TF1 *fE3J_set = 0)
  {
    fJ = fJ_set;
    fE3J = fE3J_set;
  }

  // Add the energy scale correction function as a function of log10(E/eV) and its parameters if applicable
  void SetEncorr(TF1 *fEnCorr_set = 0)
  {
    fEnCorr = fEnCorr_set;
  }

  // Set the minimum number of events per bin for calculating the restricted log likelihood
  void SetNeventsMinRestricted(Int_t nevents_min_restricted_value = 7)
  {
    nevents_min_restricted = nevents_min_restricted_value;
  }

  // contribution to the log likelihood function from this instance
  // first member of the pair is the log likelihood, second member of the pair
  // is the number of bins that are contributing
  void CalcLogLikelihood(Double_t log10en_min, Double_t log10en_max);

  // Set functions that are to be used in evaluating the null hypothesis
  void SetNullFun(TF1 *fJ_null_set, TF1 *fE3J_null_set = 0)
  {
    fJ_null = fJ_null_set;
    fE3J_null = fE3J_null_set;
  }

  // count the number of events between the minimum and maximum energies and (if the null hypothesis flux function is provided)
  // return the number of events expected from the flux function and the number of events observed in the data
  std::pair<Double_t, Double_t> EvalNull();

  // contribution to the log likelihood function from this instance
  // first member of the pair is the log likelihood, second member of the pair
  // is the number of bins that are contributing
  const std::pair<Double_t, Double_t>& GetLogLikelihood(Double_t log10en_min, Double_t log10en_max)
  {
    CalcLogLikelihood(log10en_min, log10en_max);
    return log_likelihood;
  }

  const std::pair<Double_t, Double_t>& GetLogLikelihood()
  {
    return log_likelihood;
  }
  const std::pair<Double_t, Double_t>& GetLogLikelihoodNonzero()
  {
    return log_likelihood_nonzero;
  }
  const std::pair<Double_t, Double_t>& GetLogLikelihoodRestricted()
  {
    return log_likelihood_restricted;
  }

  // graphs for the Flux, E^3 x Flux, numbers of events, fit prediction numbers of events, and the numbers
  // for the null hypothesis.
  TGraphAsymmErrors* GetJ();
  TGraphErrors *GetJ_simple_errors();
  TGraphAsymmErrors* GetE3J();
  TGraphErrors *GetE3J_simple_errors();
  TGraphAsymmErrors* GetNevents();
  TGraphErrors *GetNevents_simple_errors();
  TGraph* GetNeventsFit();
  TGraph* GetNeventsNull();

  // plotting functions
  void Plot(const char *what = "e3j", const char *draw_opt = "a,e1p");
  void Draw(Option_t *opt = "e3j,a,e1p");

  // print flux data to ASCII file
  void to_ascii_file(const char* ascii_file_name = 0);

  // To rescale the exposure (after the data has been loaded) 
  // by some constant factor. Sometimes
  // useful for displaying purposes.
  void RescaleExposure(Double_t c = 1.0);

  Double_t log10en_min_data;            // minimum energy available from the loaded data
  Double_t log10en_max_data;            // maximum energy available from the loaded data

  std::vector<Double_t> log10en;        // log10(E/eV)
  std::vector<Double_t> log10en_bsize;  // log10(E/eV) bin size
  std::vector<Double_t> nevents;        // number of events
  std::vector<Double_t> exposure;       // exposure [m^2 sr s]
  std::vector<Double_t> nevents_fit;    // fit prediction for the numbers of events

  std::vector<Int_t> bins_null;         // energy bin indices for the bins that participate in null hypothesis evaluation
  std::vector<Double_t> nevents_null;   // numbers of events expected by the null hypothesis function for those energy bins

  std::pair<Double_t, Double_t> log_likelihood;             // log likelihood and number of fit bins over all available bins
  std::pair<Double_t, Double_t> log_likelihood_nonzero;     // log likelihood and the number of fit bins for non-zero event bins
  Int_t nevents_min_restricted;                             // minimum number of events for calculating restricted log likelihood
  std::pair<Double_t, Double_t> log_likelihood_restricted;  // log likelihood and the number of fit bins for bins that meet the smallest number of events

  // Flux versus log10(E/eV) function.  If the pointer to this function is zero then zeros are returned for the log likelihood calculations
  TF1 *fJ;   // flux function that's used for the fits and displaying the results
  TF1 *fE3J; // flux function that's used for displaying the results

  TF1 *fJ_null; // flux function that's used as the null hypothesis
  TF1 *fE3J_null; // E^3 x flux function that's used as the null hypothesis

  // energy correction function that's applicable for this instance of the flux measurement
  // if this pointer is other than zero then it will be used in calculation of the log likelihood
  // contribution of the flux as well as in correcting the flux when plotting the results
  // the arguments to this function are assumed to be
  // function of log10(E/eV)
  TF1 *fEnCorr;

  // for the class dictionary generation
ClassDef(TCRFlux,1)
  ;

};

#endif
