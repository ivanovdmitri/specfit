// Dmitri Ivanov <dmiivanov@gmail.com>

#include "TCRFluxFit.h"
#include <cstdio>
#include <cstdlib>
#include "TTree.h"


ClassImp(TCRFluxFit);

TCRFluxFit::~TCRFluxFit()
{
  // since the Minuit minimizer pointer is frequently created and destroyed,
  // we keep track of it separately
  if(mFIT)
    delete mFIT;

  // be sure the clean up any TCRFlux objects that were created by this class
  if(TCRFlux_Objects_Created_By_This.GetEntries())
    {
      TIterator *itr = TCRFlux_Objects_Created_By_This.MakeIterator();
      for (TObject *obj = 0; (obj = itr->Next()); delete obj)
	;
    }
  TCRFlux_Objects_Created_By_This.Clear();
  Fluxes.clear();
  Fluxes_ordered.clear();
}

// Set the pointer to the functions that describes the fitted flux versus log10(E/eV)
// E^3 J function is optional, it's mainly used for plotting the results
void TCRFluxFit::SetFluxFun(TF1 *fJ_set, TF1 *fE3J_set)
{
  fJ = fJ_set;
  // if E^3 J function was not given then attempt to construct it from
  // the J function
  if(!fE3J_set && fJ)
    fE3J = specfit_uti::get_e3j_from_j(fJ);
  for (std::map<TString, TCRFlux*>::iterator iflux = Fluxes.begin(); iflux != Fluxes.end(); iflux++)
    {
      TCRFlux &flux = *iflux->second;
      flux.SetFluxFun(fJ, fE3J);
    }
}

// Set the pointer to the functions that describes the null hypothesis flux versus log10(E/eV)
// E^3 J function is optional, it's mainly used for plotting the results
void TCRFluxFit::SetNullFun(TF1 *fJ_null_set, TF1 *fE3J_null_set)
{
  fJ_null = fJ_null_set;
  // if E^3 J function was not given then attempt to construct it from
  // the J function
  if(!fE3J_null_set && fJ_null)
    fE3J_null = specfit_uti::get_e3j_from_j(fJ_null);
  for (std::map<TString, TCRFlux*>::iterator iflux = Fluxes.begin(); iflux != Fluxes.end(); iflux++)
    {
      TCRFlux &flux = *iflux->second;
      flux.SetNullFun(fJ_null, fE3J_null);
    }
}

// add the measured flux result plus the energy correction function that corresponds to that flux
// this function assumes that the flux has been initialized elsewhere and the class will not attempt
// to clean it up
Bool_t TCRFluxFit::Add(TCRFlux *flux, TF1 *fEnCorr_set)
{
  if(!flux)
    {
      fprintf(stderr, "ERROR: attempting to add flux pointer that hasn't been initialized!\n");
      return false;
    }
  TString name = flux->GetName();
  if(Fluxes.find(name) != Fluxes.end())
    {
      fprintf(stderr, "WARNING: flux result '%s' has been already added\n", name.Data());
      return false;
    }
  TCRFlux &added_flux = *(Fluxes[flux->GetName()] = flux);
  added_flux.SetFluxFun(fJ, fE3J);
  added_flux.SetEncorr(fEnCorr_set);
  if(fEnCorr_set)
    fEnCorr[added_flux.GetName()] = fEnCorr_set; // keep track of pointers to all non-zero energy correction functions
  Fluxes_ordered.push_back(&added_flux);         // store pointers to fluxes ordered by addition
  return true;
}

// nbins: number of bins
// log10en_values: energies log10(E/eV) of the bin centers
// log10en_bsize_values: log10(E/eV) bin sizes
// nevents_values: numbers of events in each energy bin
// exposure_values: exposure [m^2 sr s] for each energy bin center value
Bool_t TCRFluxFit::Add(const char *name, const char *title, Int_t nbins, const Double_t *log10en_values, const Double_t *log10en_bsize_values, const Double_t *nevents_values,
    const Double_t *exposure_values, TF1 *fEnCorr_set)
{
  if(Fluxes.find(name) != Fluxes.end())
    {
      fprintf(stderr, "WARNING: flux result named '%s' has been already added\n", name);
      return false;
    }
  TCRFlux &flux = *(Fluxes[name] = new TCRFlux(name, title));
  flux.Load(nbins, log10en_values, log10en_bsize_values, nevents_values, exposure_values);
  flux.SetFluxFun(fJ, fE3J);
  flux.SetEncorr(fEnCorr_set);
  if(fEnCorr_set)
    fEnCorr[flux.GetName()] = fEnCorr_set; // keep track of pointers to all non-zero energy correction functions
  Fluxes_ordered.push_back(&flux);         // store pointers to fluxes ordered by addition
  // keep track of all flux objects that we are creating
  TCRFlux_Objects_Created_By_This.Add(&flux);
  return true;
}

// add a flux result from an ASCII file that has:
// col1: energies log10(E/eV) of the bin centers
// col2: log10(E/eV) bin sizes
// col3: numbers of events in each energy bin
// col4: exposure [m^2 sr s] for each energy bin center value
Bool_t TCRFluxFit::Add(const char *name, const char *title, const char *ascii_file, TF1 *fEnCorr_set)
{
  TTree *t = new TTree("t", "");
  if(!t->ReadFile(ascii_file, "log10en/D:log10en_bsize/D:nevents/D:exposure/D"))
    {
      delete t;
      return false;
    }
  t->Draw("log10en:log10en_bsize:nevents:exposure", "", "goff");
  if(!Add(name, title, t->GetSelectedRows(), t->GetV1(), t->GetV2(), t->GetV3(), t->GetV4(), fEnCorr_set))
    {
      delete t;
      return false;
    }
  delete t;
  return true;
}

void TCRFluxFit::SelectEnergyRange(Double_t log10en_min, Double_t log10en_max)
{
  for (std::map<TString, TCRFlux*>::iterator iflux = Fluxes.begin(); iflux != Fluxes.end(); iflux++)
    {
      TCRFlux &flux = *iflux->second;
      flux.SelectEnergyRange(log10en_min, log10en_max);
    }
}

// return stored flux result
TCRFlux* TCRFluxFit::GetFlux(const char *name)
{
  std::map<TString, TCRFlux*>::iterator iflux = Fluxes.find(name);
  if(iflux == Fluxes.end())
    {
      fprintf(stderr, "ERROR: flux named '%s' not found!\n", name);
      return 0;
    }
  return iflux->second;
}

// return stored flux result by index
TCRFlux* TCRFluxFit::GetFlux(Int_t iflux)
{
  if(iflux < 0 || iflux > (Int_t) Fluxes_ordered.size() - 1)
    {
      fprintf(stderr, "error: iflux must be in 0 to %d range\n", (Int_t) Fluxes_ordered.size() - 1);
      return 0;
    }
  return Fluxes_ordered[iflux];
}

// Set the parameters for the flux functions and optionally errors on parameters
// (errors may be useful for displaying the results)
void TCRFluxFit::SetFluxPar(const Double_t *params, const Double_t *parerrors)
{
  if(fJ)
    fJ->SetParameters(params);
  if(fE3J)
    fE3J->SetParameters(params);
  if(parerrors)
    {
      if(fJ)
	fJ->SetParErrors(parerrors);
      if(fE3J)
	fE3J->SetParErrors(parerrors);
    }
}

// All correction functions must use the same parameters (but use them in different ways).
// Set the parameters for the correction functions and optionally errors on parameters
// (errors may be useful for displaying the results)
void TCRFluxFit::SetEncorrPar(const Double_t *params, const Double_t *parerrors)
{
  if(!fEnCorr.size())
    return;
  for (std::map<TString, TF1*>::iterator i = fEnCorr.begin(); i != fEnCorr.end(); i++)
    {
      TF1 *&f = i->second;
      if(f)
	{
	  f->SetParameters(params);
	  if(parerrors)
	    f->SetParErrors(parerrors);
	}
    }
}

// To set the parameters while fitting.
// This depends on whether just the flux or the flux and
// the nonlinearity correction are being fitted.  Parameters for the nonlinearity
// correction follow those of the flux in the "par" array
// that's being adjusted by Minuit minimizer while finding the minimum.
void TCRFluxFit::SetParameters(const Double_t *par)
{
  if(nfluxpar)
    SetFluxPar(par);
  if(nencorrpar)
    SetEncorrPar(par + nfluxpar);
}

// calculate the overall log likelihood
void TCRFluxFit::CalcLogLikelihood()
{
  log_likelihood = std::make_pair(0, 0);
  log_likelihood_nonzero = std::make_pair(0, 0);
  log_likelihood_restricted = std::make_pair(0, 0);
  for (std::map<TString, TCRFlux*>::iterator iflux = Fluxes.begin(); iflux != Fluxes.end(); iflux++)
    {
      TCRFlux &flux = *iflux->second;
      flux.CalcLogLikelihood(log10en_min, log10en_max);
      log_likelihood.first += flux.log_likelihood.first;
      log_likelihood.second += flux.log_likelihood.second;
      log_likelihood_nonzero.first += flux.log_likelihood_nonzero.first;
      log_likelihood_nonzero.second += flux.log_likelihood_nonzero.second;
      log_likelihood_restricted.first += flux.log_likelihood_restricted.first;
      log_likelihood_restricted.second += flux.log_likelihood_restricted.second;
    }
}

// calculate the numbers of events within a certain energy range that are expected from some null hypothesis flux function
// and that are actually observed in the data
std::pair<Double_t, Double_t> TCRFluxFit::EvalNull()
{
  std::pair<Double_t, Double_t> nexpected_nobserved = std::make_pair(0, 0);
  for (std::map<TString, TCRFlux*>::iterator iflux = Fluxes.begin(); iflux != Fluxes.end(); iflux++)
    {
      TCRFlux &flux = *iflux->second;
      std::pair<Double_t, Double_t> x = flux.EvalNull();
      nexpected_nobserved.first += x.first;
      nexpected_nobserved.second += x.second;
    }
  return nexpected_nobserved;
}

static TCRFluxFit *pointer_to_global_instance_of_TCRFluxFit = 0;
static void fcn_for_mFIT(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  (void) (npar);
  (void) (gin);
  (void) (iflag);
  f = 0;
  if(pointer_to_global_instance_of_TCRFluxFit)
    {
      TCRFluxFit &fit = (*pointer_to_global_instance_of_TCRFluxFit);
      fit.SetParameters(par);
      f = fit.GetLogLikelihood().first;
    }
}

Bool_t TCRFluxFit::Fit()
{
  if(!Fluxes.size())
    {
      fprintf(stderr, "ERROR: add some flux results before fitting!\n");
      return false;
    }
  if(!fJ)
    {
      fprintf(stderr, "ERROR: set the flux function before fitting!\n");
      return false;
    }

  nfluxpar = fJ->GetNpar(); // contribution to the number of fit parameters from the flux fit function

  // additional fit parameters if there are energy correction functions involved
  // note: if the energy correction functions are involved then they must use the same
  // parameters
  TF1 *fEnCorr_first = (fEnCorr.size() ? fEnCorr.begin()->second : 0);
  nencorrpar = (fEnCorr_first ? fEnCorr_first->GetNpar() : 0);

  // Initialize the Minuit minimizer
  nfitpar = nfluxpar + nencorrpar;
  if(mFIT)
    delete mFIT;
  mFIT = new TMinuit(nfitpar);

  // declare the parameters
  for (Int_t i = 0; i < nfluxpar; i++)
    {
      Double_t parmin = 0, parmax = 0;
      fJ->GetParLimits(i, parmin, parmax);
      mFIT->DefineParameter(i, fJ->GetParName(i), fJ->GetParameter(i), fJ->GetParError(i), parmin, parmax);
    }
  // declare the parameters that correspond to the energy correction function
  for (Int_t i = 0; i < nencorrpar; i++)
    {
      Double_t parmin = 0, parmax = 0;
      fEnCorr_first->GetParLimits(i, parmin, parmax);
      mFIT->DefineParameter(nfluxpar + i, fEnCorr_first->GetParName(i), fEnCorr_first->GetParameter(i), fEnCorr_first->GetParError(i), parmin, parmax);
    }

  // function to minimize
  pointer_to_global_instance_of_TCRFluxFit = this;
  mFIT->SetFCN(fcn_for_mFIT);

  // We expect that the change of -2 *  log (likelihood) by 1 will correspond to 1 sigma errors
  mFIT->SetErrorDef(1.0);

  // Perform minimization
  mFIT->Migrad();

  // Get the best fit parameters
  std::vector<Double_t> params(nfitpar, 0), parerrors(nfitpar, 0);
  for (Int_t i = 0; i < nfitpar; i++)
    mFIT->GetParameter(i, params[i], parerrors[i]);

  // set the best fit parameters to the corresponding functions
  SetFluxPar(&params[0], &parerrors[0]); // flux fit function
  // energy correction function, if correction parameters are fitted
  if(nencorrpar)
    SetEncorrPar(&params[nfluxpar], &parerrors[nfluxpar]);

  // calculate the overall log likelihood again, using the best fit parameters
  CalcLogLikelihood();

  // calculate the chi2 = (normalized log likelihood) and the
  // number of degrees of freedom = (number of fitted bins) - (total number of fit parameters)
  chi2 = log_likelihood.first;
  ndof = log_likelihood.second - (Double_t) nfitpar;

  // return success
  return true;
}

TMinuit* TCRFluxFit::GetMinuit()
{
  // To obtain the Minuit pointer for whatever reason.
  // In order for it to behave correctly, the global FCN must be
  // pointed to use this instance of the class.
  pointer_to_global_instance_of_TCRFluxFit = this;
  return mFIT;
}

