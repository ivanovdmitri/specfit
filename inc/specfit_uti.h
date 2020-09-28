// Dmitri Ivanov <dmiivanov@gmail.com>
#ifndef _specfit_uti_h_
#define _specfit_uti_h_

#include "TObject.h"
#include "TString.h"
#include "TF1.h"
#include <algorithm>

namespace specfit_uti
{
  // To get a unique name for an object
  TString get_unique_object_name(const char *basename);

  // to get Feldman-Cousins error bars
  std::pair<Double_t, Double_t> get_fc_errors(Double_t n);

  // to get the lower error bar
  Double_t get_fc_error_low(Double_t n);

  // to get the upper error bar
  Double_t get_fc_error_high(Double_t n);

  // linear bin size
  Double_t GetLinBinSize(Double_t log10en, Double_t log10en_bsize);

  // get the significance in sigma units
  // from the chance probability
  Double_t pchance2sigma(Double_t pchance, Bool_t pwarning = true);

  // express chance probability in sigma units
  Double_t Sigma2Pchance(Double_t pchange_in_sigma);

  // get the chance probability of a Poisson fluctuation
  Double_t PoissonPchance(Int_t nobserved, Double_t nexpected, Bool_t in_sigma_units = true);

  // to obtain E^{3}J function from J if J was constructed using formula
  TF1* get_e3j_from_j(TF1 *f_J);

}
#endif
