// Dmitri Ivanov <dmiivanov@gmail.com>

#include <cstdio>
#include <cstdlib>
#include "TCRFlux.h"
#include "TMath.h"
#include "TTree.h"
#include "TROOT.h"
#include "specfit_uti.h"
#include "TAxis.h"

// for the class dictionary generation
ClassImp(TCRFlux);

TCRFlux::TCRFlux(const char *name, const char *title) :
    log10en_min_data(0), log10en_max_data(0), nevents_min_restricted(7), fJ(0), fE3J(0), fJ_null(0), fE3J_null(0), fEnCorr(0)
{
  SetName(name);
  SetTitle(title);
  gROOT->GetListOfSpecials()->Add(this); // add the object to the list of ROOT's specials
}

TCRFlux::~TCRFlux()
{
  gROOT->GetListOfSpecials()->Remove(this);
}

// Load data from an ASCII file with columns
// col1: energies log10(E/eV) of the bin centers
// col2: log10(E/eV) bin sizes
// col3: numbers of events in each energy bin
// col4: exposure [m^2 sr s] for each energy bin center value
Bool_t TCRFlux::Load(const char *ascii_file)
{
  TTree *t = new TTree("t", "");
  if(!t->ReadFile(ascii_file, "log10en/D:log10en_bsize/D:nevents/D:exposure/D"))
    {
      delete t;
      return false;
    }
  t->Draw("log10en:log10en_bsize:nevents:exposure", "", "goff");
  if(!Load(t->GetSelectedRows(), t->GetV1(), t->GetV2(), t->GetV3(), t->GetV4()))
    {
      delete t;
      return false;
    }
  delete t;
  return true;
}

// load data from C-like arrays
Bool_t TCRFlux::Load(Int_t nbins,                   // number of bins for the spectrum
    const Double_t *log10en_values,        // energies log10(E/eV) of the bin centers
    const Double_t *log10en_bsize_values,  // log10(E/eV) bin sizes
    const Double_t *nevents_values,        // numbers of events in each energy bin
    const Double_t *exposure_values)       // exposure [m^2 sr s] for each energy bin center value
{
  log10en = std::vector<Double_t>(log10en_values, log10en_values + nbins);
  log10en_bsize = std::vector<Double_t>(log10en_bsize_values, log10en_bsize_values + nbins);
  nevents = std::vector<Double_t>(nevents_values, nevents_values + nbins);
  exposure = std::vector<Double_t>(exposure_values, exposure_values + nbins);
  nevents_fit = std::vector<Double_t>(nevents.size(), 0);
  find_min_max_log10en();
  return true;
}

// load data from STL vectors
Bool_t TCRFlux::Load(const std::vector<Double_t> &log10en_values,  // energies log10(E/eV) of the bin centers
    const std::vector<Double_t> &log10en_bsize_values,    // log10(E/eV) bin sizes
    const std::vector<Double_t> &nevents_values,          // numbers of events in each energy bin
    const std::vector<Double_t> &exposure_values)         // exposure [m^2 sr s] for each energy bin center value
{
  if(log10en_values.size() != log10en_bsize_values.size() || log10en_values.size() != nevents_values.size() || log10en_values.size() != exposure_values.size())
    {
      fprintf(stderr, "ERROR: Load: sizes of the arrays are not the same!\n");
      return false;
    }
  log10en = log10en_values;
  log10en_bsize = log10en_bsize_values;
  nevents = nevents_values;
  exposure = exposure_values;
  nevents_fit = std::vector<Double_t>(nevents.size(), 0);
  find_min_max_log10en();
  return true;
}

// determine the energy range of the spectrum measurement
void TCRFlux::find_min_max_log10en()
{
  if(!log10en.size())
    {
      log10en_min_data = 0;
      log10en_max_data = 0;
    }
  Int_t imin = 0, imax = 0;
  for (Int_t i = 0; i < (Int_t) log10en.size(); i++)
    {
      if(log10en[i] < log10en[imin])
	imin = i;
      if(log10en[i] > log10en[imax])
	imax = i;
    }
  log10en_min_data = log10en[imin] - log10en_bsize[imin] / 2.0;
  log10en_max_data = log10en[imax] + log10en_bsize[imax] / 2.0;
}

// this selects the data only within the desirable energy range
void TCRFlux::SelectEnergyRange(Double_t log10en_min, Double_t log10en_max)
{
  for (Int_t i = 0; i < (Int_t) log10en.size(); i++)
    {
      if(log10en[i] < log10en_min || log10en[i] > log10en_max)
	{
	  log10en.erase(log10en.begin() + i);
	  log10en_bsize.erase(log10en_bsize.begin() + i);
	  nevents.erase(nevents.begin() + i);
	  exposure.erase(exposure.begin() + i);
	  nevents_fit.erase(nevents_fit.begin() + i);
	  i--;
	}
    }
  find_min_max_log10en();
}

// contribution to the log likelihood function from this instance
// first member of the pair is the log likelihood, second member of the pair
// is the number of bins that are contributing
void TCRFlux::CalcLogLikelihood(Double_t log10en_min, Double_t log10en_max)
{

  log_likelihood = std::make_pair(0, 0);
  log_likelihood_nonzero = std::make_pair(0, 0);
  log_likelihood_restricted = std::make_pair(0, 0);

  for (Int_t i = 0; i < (Int_t) log10en.size(); i++)
    {
      // apply the energy limits
      if(log10en[i] < log10en_min || log10en[i] > log10en_max)
	continue;
      // number of events from the fit function
      Double_t bsize = specfit_uti::GetLinBinSize(log10en[i], log10en_bsize[i]);
      // correct the fit predictions appropriately if the energy correction function is being applied
      Double_t encorr = (fEnCorr ? fEnCorr->Eval(log10en[i]) : 1.0);
      Double_t log10en_corr = log10en[i] + TMath::Log10(encorr);
      Double_t lgl = 0; // contribution to log likelihood from the bin
      if(fJ)
	{
	  nevents_fit[i] = fJ->Eval(log10en_corr) * (encorr * bsize) * exposure[i];
	  // log likelihood formula when number of events is zero
	  lgl = 2.0 * nevents_fit[i];
	  // log likelihood formula when number of events is not zero
	  if(nevents[i] > 1e-3)
	    lgl = 2.0 * ((nevents_fit[i] - nevents[i]) + nevents[i] * TMath::Log(nevents[i] / nevents_fit[i]));
	}
      // if the flux function was never given then set the fit prediction numbers of events and
      // log likelihoods to zeros
      else
	{
	  nevents_fit[i] = 0;
	  lgl = 0;
	}

      log_likelihood.first += lgl;
      log_likelihood.second++;

      if(nevents[i] > 0)
	{
	  log_likelihood_nonzero.first += lgl;
	  log_likelihood_nonzero.second++;
	}
      if(nevents[i] >= nevents_min_restricted)
	{
	  log_likelihood_restricted.first += lgl;
	  log_likelihood_restricted.second++;
	}
    }
}
// count the number of events between the minimum and maximum energies and (if the null hypothesis flux function is provided)
// return the number of events expected from the flux function and the number of events observed in the data
std::pair<Double_t, Double_t> TCRFlux::EvalNull()
{
  std::pair<Double_t, Double_t> nexpect_nobserve = std::make_pair(0, 0);
  bins_null.clear();
  nevents_null.clear();
  if(!fJ_null)
    return nexpect_nobserve;
  // determine the energy scale correction for the experiment and translate the
  // energies to experiment's energy scale
  Double_t encorr_en_min = (fEnCorr ? fEnCorr->Eval(fJ_null->GetXmin()) : 1.0);
  Double_t encorr_en_max = (fEnCorr ? fEnCorr->Eval(fJ_null->GetXmax()) : 1.0);
  Double_t log10en_min_corr = fJ_null->GetXmin() + TMath::Log10(encorr_en_min);
  Double_t log10en_max_corr = fJ_null->GetXmax() + TMath::Log10(encorr_en_max);
  for (Int_t i = 0; i < (Int_t) log10en.size(); i++)
    {
      // apply the energy limits with energy correction that's appropriate for the experiment
      if(log10en[i] < log10en_min_corr || log10en[i] > log10en_max_corr)
	continue;
      Double_t bsize = specfit_uti::GetLinBinSize(log10en[i], log10en_bsize[i]);
      // correct the predictions appropriately if the energy correction function is being applied
      Double_t encorr = (fEnCorr ? fEnCorr->Eval(log10en[i]) : 1.0);
      Double_t log10en_corr = log10en[i] + TMath::Log10(encorr);
      // number of events expected from the given flux function
      Double_t nexpect = fJ_null->Eval(log10en_corr) * (encorr * bsize) * exposure[i];
      bins_null.push_back(i);
      nevents_null.push_back(nexpect);
      nexpect_nobserve.first += nexpect;
      // count the data events that are within the energy range
      nexpect_nobserve.second += nevents[i];
    }
  // return the answer
  return nexpect_nobserve;
}

TGraphAsymmErrors* TCRFlux::GetJ()
{
  TGraphAsymmErrors *g = new TGraphAsymmErrors(nevents.size());
  g->SetName(specfit_uti::get_unique_object_name(TString("g") + TString(GetName()) + "_J"));
  g->SetTitle(TString::Format("%s;log_{10}(E/eV);J [ eV^{-1} m^{-2} sr^{-1} s^{-1} ]", GetTitle()));
  g->GetXaxis()->SetTitleSize(0.055);
  g->GetYaxis()->SetTitleSize(0.055);
  g->SetMarkerStyle(20);
  Double_t ylow = 1e256, yhigh = -1;
  for (Int_t i = 0; i < (Int_t) nevents.size(); i++)
    {
      Double_t bsize = specfit_uti::GetLinBinSize(log10en[i], log10en_bsize[i]);
      Double_t encorr = (fEnCorr ? fEnCorr->Eval(log10en[i]) : 1.0);
      Double_t j = 1.0 / encorr / bsize / exposure[i] * nevents[i];
      Double_t j_e1 = 1.0 / encorr / bsize / exposure[i] * specfit_uti::get_fc_error_low(nevents[i]);
      Double_t j_e2 = 1.0 / encorr / bsize / exposure[i] * specfit_uti::get_fc_error_high(nevents[i]);
      g->SetPoint(i, log10en[i] + TMath::Log10(encorr), j);
      g->SetPointError(i, 0, 0, j_e1, j_e2);
      if(nevents[i] > 0.5)
	{
	  if(0.9 * (j - j_e1) < ylow)
	    ylow = 0.9 * (j - j_e1);
	  if(1.1 * (j + j_e2) > yhigh)
	    yhigh = 1.1 * (j + j_e2);
	}
    }
  if(fJ)
    {
      Double_t y[2] =
      { 0.9 * fJ->GetMinimum(), 1.1 * fJ->GetMaximum() };
      if(ylow > y[0])
	ylow = y[0];
      if(yhigh < y[1])
	yhigh = y[1];
    }
  if(fJ_null)
    {
      Double_t y[2] =
      { 0.9 * fJ_null->GetMinimum(), 1.1 * fJ_null->GetMaximum() };
      if(ylow > y[0])
	ylow = y[0];
      if(yhigh < y[1])
	yhigh = y[1];
    }
  g->GetYaxis()->SetRangeUser(ylow, yhigh);
  return g;
}

TGraphErrors* TCRFlux::GetJ_simple_errors()
{
  TGraphErrors *g = new TGraphErrors(nevents.size());
  g->SetName(specfit_uti::get_unique_object_name(TString("g") + TString(GetName()) + "_J"));
  g->SetTitle(TString::Format("%s;log_{10}(E/eV);J [ eV^{-1} m^{-2} sr^{-1} s^{-1} ]", GetTitle()));
  g->GetXaxis()->SetTitleSize(0.055);
  g->GetYaxis()->SetTitleSize(0.055);
  g->SetMarkerStyle(20);
  Double_t ylow = 1e256, yhigh = -1;
  for (Int_t i = 0; i < (Int_t) nevents.size(); i++)
    {
      Double_t bsize = specfit_uti::GetLinBinSize(log10en[i], log10en_bsize[i]);
      Double_t encorr = (fEnCorr ? fEnCorr->Eval(log10en[i]) : 1.0);
      Double_t j = 1.0 / encorr / bsize / exposure[i] * nevents[i];
      Double_t j_e1 = 1.0 / encorr / bsize / exposure[i] * TMath::Sqrt(nevents[i]);
      g->SetPoint(i, log10en[i] + TMath::Log10(encorr), j);
      g->SetPointError(i, 0, j_e1);
      if(nevents[i] > 0.5)
	{
	  if(0.9 * (j - j_e1) < ylow)
	    ylow = 0.9 * (j - j_e1);
	  if(1.1 * (j + j_e1) > yhigh)
	    yhigh = 1.1 * (j + j_e1);
	}
    }
  if(fJ)
    {
      Double_t y[2] =
      { 0.9 * fJ->GetMinimum(), 1.1 * fJ->GetMaximum() };
      if(ylow > y[0])
	ylow = y[0];
      if(yhigh < y[1])
	yhigh = y[1];
    }
  if(fJ_null)
    {
      Double_t y[2] =
      { 0.9 * fJ_null->GetMinimum(), 1.1 * fJ_null->GetMaximum() };
      if(ylow > y[0])
	ylow = y[0];
      if(yhigh < y[1])
	yhigh = y[1];
    }
  g->GetYaxis()->SetRangeUser(ylow, yhigh);
  return g;
}

TGraphAsymmErrors* TCRFlux::GetE3J()
{
  TGraphAsymmErrors *g = new TGraphAsymmErrors(nevents.size());
  g->SetName(specfit_uti::get_unique_object_name(TString("g") + TString(GetName()) + "_E3J"));
  g->SetTitle(TString::Format("%s;log_{10}(E/eV);E^{3} J [ eV^{-2} m^{-2} sr^{-1} s^{-1} ]", GetTitle()));
  g->GetXaxis()->SetTitleSize(0.055);
  g->GetYaxis()->SetTitleSize(0.055);
  g->SetMarkerStyle(20);
  Double_t ylow = 1e256, yhigh = -1;
  for (Int_t i = 0; i < (Int_t) nevents.size(); i++)
    {
      Double_t bsize = specfit_uti::GetLinBinSize(log10en[i], log10en_bsize[i]);
      Double_t e3 = TMath::Power(10.0, 3.0 * log10en[i]);
      Double_t encorr = (fEnCorr ? fEnCorr->Eval(log10en[i]) : 1.0);
      Double_t e3j = encorr * encorr * e3 / bsize / exposure[i] * nevents[i];
      Double_t e3j_e1 = encorr * encorr * e3 / bsize / exposure[i] * specfit_uti::get_fc_error_low(nevents[i]);
      Double_t e3j_e2 = encorr * encorr * e3 / bsize / exposure[i] * specfit_uti::get_fc_error_high(nevents[i]);
      g->SetPoint(i, log10en[i] + TMath::Log10(encorr), e3j);
      g->SetPointError(i, 0, 0, e3j_e1, e3j_e2);
      if(nevents[i] > 0.5)
	{
	  if(0.9 * (e3j - e3j_e1) < ylow)
	    ylow = 0.9 * (e3j - e3j_e1);
	  if(1.1 * (e3j + e3j_e2) > yhigh)
	    yhigh = 1.1 * (e3j + e3j_e2);
	}
    }
  if(fE3J)
    {
      Double_t y[2] =
      { 0.9 * fE3J->GetMinimum(), 1.1 * fE3J->GetMaximum() };
      if(ylow > y[0])
	ylow = y[0];
      if(yhigh < y[1])
	yhigh = y[1];
    }
  if(fE3J_null)
    {
      Double_t y[2] =
      { 0.9 * fE3J_null->GetMinimum(), 1.1 * fE3J_null->GetMaximum() };
      if(ylow > y[0])
	ylow = y[0];
      if(yhigh < y[1])
	yhigh = y[1];
    }
  g->GetYaxis()->SetRangeUser(ylow, yhigh);
  return g;
}

TGraphErrors* TCRFlux::GetE3J_simple_errors()
{
  TGraphErrors *g = new TGraphErrors(nevents.size());
  g->SetName(specfit_uti::get_unique_object_name(TString("g") + TString(GetName()) + "_E3J"));
  g->SetTitle(TString::Format("%s;log_{10}(E/eV);E^{3} J [ eV^{-2} m^{-2} sr^{-1} s^{-1} ]", GetTitle()));
  g->GetXaxis()->SetTitleSize(0.055);
  g->GetYaxis()->SetTitleSize(0.055);
  g->SetMarkerStyle(20);
  Double_t ylow = 1e256, yhigh = -1;
  for (Int_t i = 0; i < (Int_t) nevents.size(); i++)
    {
      Double_t bsize = specfit_uti::GetLinBinSize(log10en[i], log10en_bsize[i]);
      Double_t e3 = TMath::Power(10.0, 3.0 * log10en[i]);
      Double_t encorr = (fEnCorr ? fEnCorr->Eval(log10en[i]) : 1.0);
      Double_t e3j = encorr * encorr * e3 / bsize / exposure[i] * nevents[i];
      Double_t e3j_e1 = encorr * encorr * e3 / bsize / exposure[i] * TMath::Sqrt(nevents[i]);
      g->SetPoint(i, log10en[i] + TMath::Log10(encorr), e3j);
      g->SetPointError(i, 0, e3j_e1);
      if(nevents[i] > 0.5)
	{
	  if(0.9 * (e3j - e3j_e1) < ylow)
	    ylow = 0.9 * (e3j - e3j_e1);
	  if(1.1 * (e3j + e3j_e1) > yhigh)
	    yhigh = 1.1 * (e3j + e3j_e1);
	}
    }
  if(fE3J)
    {
      Double_t y[2] =
      { 0.9 * fE3J->GetMinimum(), 1.1 * fE3J->GetMaximum() };
      if(ylow > y[0])
	ylow = y[0];
      if(yhigh < y[1])
	yhigh = y[1];
    }
  if(fE3J_null)
    {
      Double_t y[2] =
      { 0.9 * fE3J_null->GetMinimum(), 1.1 * fE3J_null->GetMaximum() };
      if(ylow > y[0])
	ylow = y[0];
      if(yhigh < y[1])
	yhigh = y[1];
    }
  g->GetYaxis()->SetRangeUser(ylow, yhigh);
  return g;
}


TGraphAsymmErrors* TCRFlux::GetNevents()
{
  TGraphAsymmErrors *g = new TGraphAsymmErrors(nevents.size());
  g->SetName(specfit_uti::get_unique_object_name(TString("g") + TString(GetName()) + "_N"));
  g->SetTitle(TString::Format("%s;log_{10}(E/eV);N_{EVENTS} / BIN", GetTitle()));
  g->GetXaxis()->SetTitleSize(0.055);
  g->GetYaxis()->SetTitleSize(0.055);
  g->SetMarkerStyle(20);
  Double_t ylow = 1e256, yhigh = -1;
  for (Int_t i = 0; i < (Int_t) nevents.size(); i++)
    {
      Double_t n_e1 = specfit_uti::get_fc_error_low(nevents[i]);
      Double_t n_e2 = specfit_uti::get_fc_error_high(nevents[i]);
      Double_t encorr = (fEnCorr ? fEnCorr->Eval(log10en[i]) : 1.0);
      g->SetPoint(i, log10en[i] + TMath::Log10(encorr), nevents[i]);
      g->SetPointError(i, log10en_bsize[i] / 2.0, log10en_bsize[i] / 2.0, n_e1, n_e2);
      if(nevents[i] > 0.5)
	{
	  if(0.9 * (nevents[i] - n_e1) < ylow)
	    ylow = 0.9 * (nevents[i] - n_e1);
	  if(1.1 * (nevents[i] + n_e2) > yhigh)
	    yhigh = 1.1 * (nevents[i] + n_e2);
	}
    }
  g->GetYaxis()->SetRangeUser(ylow, yhigh);
  return g;
}

TGraphErrors* TCRFlux::GetNevents_simple_errors()
{
  TGraphErrors *g = new TGraphErrors(nevents.size());
  g->SetName(specfit_uti::get_unique_object_name(TString("g") + TString(GetName()) + "_N"));
  g->SetTitle(TString::Format("%s;log_{10}(E/eV);N_{EVENTS} / BIN", GetTitle()));
  g->GetXaxis()->SetTitleSize(0.055);
  g->GetYaxis()->SetTitleSize(0.055);
  g->SetMarkerStyle(20);
  Double_t ylow = 1e256, yhigh = -1;
  for (Int_t i = 0; i < (Int_t) nevents.size(); i++)
    {
      Double_t n_e1 = TMath::Sqrt(nevents[i]);
      Double_t encorr = (fEnCorr ? fEnCorr->Eval(log10en[i]) : 1.0);
      g->SetPoint(i, log10en[i] + TMath::Log10(encorr), nevents[i]);
      g->SetPointError(i, log10en_bsize[i] / 2.0, n_e1);
      if(nevents[i] > 0.5)
	{
	  if(0.9 * (nevents[i] - n_e1) < ylow)
	    ylow = 0.9 * (nevents[i] - n_e1);
	  if(1.1 * (nevents[i] + n_e1) > yhigh)
	    yhigh = 1.1 * (nevents[i] + n_e1);
	}
    }
  g->GetYaxis()->SetRangeUser(ylow, yhigh);
  return g;
}

TGraph* TCRFlux::GetNeventsFit()
{
  TGraph *g = new TGraph(nevents_fit.size());
  g->SetName(specfit_uti::get_unique_object_name(TString("g") + TString(GetName()) + "_N_fit"));
  g->SetTitle(TString::Format("%s;log_{10}(E/eV);N_{EVENTS}^{FIT} / BIN", GetTitle()));
  g->GetXaxis()->SetTitleSize(0.055);
  g->GetYaxis()->SetTitleSize(0.055);
  g->SetMarkerStyle(20);
  for (Int_t i = 0; i < (Int_t) nevents_fit.size(); i++)
    {
      Double_t encorr = (fEnCorr ? fEnCorr->Eval(log10en[i]) : 1.0);
      g->SetPoint(i, log10en[i] + TMath::Log10(encorr), nevents_fit[i]);
    }
  if(fJ)
    {
      g->SetLineStyle(fJ->GetLineStyle());
      g->SetLineColor(fJ->GetLineColor());
    }
  return g;
}

TGraph* TCRFlux::GetNeventsNull()
{
  TGraph *g = new TGraph(bins_null.size());
  g->SetName(specfit_uti::get_unique_object_name(TString("g") + TString(GetName()) + "_N_null"));
  g->SetTitle(TString::Format("%s;log_{10}(E/eV);N_{EVENTS}^{NULL} / BIN", GetTitle()));
  g->GetXaxis()->SetTitleSize(0.055);
  g->GetYaxis()->SetTitleSize(0.055);
  g->SetMarkerStyle(20);
  for (Int_t i = 0; i < (Int_t) bins_null.size(); i++)
    {
      Double_t encorr = (fEnCorr ? fEnCorr->Eval(log10en[bins_null[i]]) : 1.0);
      g->SetPoint(i, log10en[bins_null[i]] + TMath::Log10(encorr), nevents_null[i]);
    }
  if(fJ_null)
    {
      g->SetLineStyle(fJ_null->GetLineStyle());
      g->SetLineColor(fJ_null->GetLineColor());
    }
  return g;
}

void TCRFlux::Plot(const char *what, const char *draw_opt)
{
  TString s_what(what);
  s_what.ToLower();
  if(s_what == "e3j")
    {
      GetE3J()->Draw(draw_opt);
      if(fE3J)
	fE3J->Draw("same");
      if(fE3J_null)
	fE3J_null->Draw("same");
    }
  else if(s_what == "" || s_what == "j")
    {
      GetJ()->Draw(draw_opt);
      if(fJ)
	fJ->Draw("same");
      if(fJ_null)
	fJ_null->Draw("same");
    }
  else if(s_what == "n" || s_what == "nevent" || s_what == "nevents")
    {
      if(nevents.size())
	GetNevents()->Draw(draw_opt);
      if(nevents_fit.size())
	GetNeventsFit()->Draw("L,same");
      if(nevents_null.size())
	GetNeventsNull()->Draw("L,same");
    }
  else
    {
      fprintf(stderr, "error: what = '%s' not recognized\n", what);
      return;
    }
}

void TCRFlux::Draw(Option_t *opt)
{
  TString s_opt(opt);
  s_opt.ToLower();
  TString what = "";
  if(s_opt.Contains("e3j"))
    what = "e3j";
  else if(!s_opt.Contains("e3j") && s_opt.Contains("j"))
    what = "j";
  else if(!s_opt.Contains("j") && (s_opt.Contains("n") || s_opt.Contains("nevent")))
    what = "n";
  else
    {
      fprintf(stderr, "option string must include \"e3j\", or \"j\", or \"n\" in it\n");
      return;
    }
  TString draw_opt = s_opt;
  while (draw_opt.Contains("e3j"))
    draw_opt.Remove(draw_opt.First("e3j"), TString("e3j").Length());
  while (draw_opt.Contains("e3j"))
    draw_opt.Remove(draw_opt.First("e3j"), TString("e3j").Length());
  while (draw_opt.Contains("nevent"))
    draw_opt.Remove(draw_opt.First("nevent"), TString("nevent").Length());
  while (draw_opt.Contains("n"))
    draw_opt.Remove(draw_opt.First("n"), TString("n").Length());
  while (draw_opt.BeginsWith(","))
    draw_opt.Remove(0, 1);
  Plot(what, draw_opt);
}

void TCRFlux::RescaleExposure(Double_t c)
{
  // To rescale the exposure (after the data has been loaded) 
  // by some constant factor. Sometimes
  // useful for displaying purposes.
  for (std::vector<Double_t>::iterator it = exposure.begin(); it != exposure.end(); it++)
    (*it) *= c;
}

void TCRFlux::to_ascii_file(const char* ascii_file_name)
{
  // print flux data to ASCII file
  TString s_ascii_file_name = ascii_file_name ? TString(ascii_file_name) : TString(GetName()) + ".txt";
  FILE *fp = fopen(s_ascii_file_name,"w");
  if(!fp)
    {
      fprintf(stderr,"ERROR: failed to open %s for writing!\n", s_ascii_file_name.Data());
      return;
    }
  fprintf(fp,"%s %12s %13s %15s\n",
          "#log10en",
          "log10en_bsize",
          "nevents",
          "exposure");
  for(Int_t ibin=0; ibin < (Int_t)log10en.size(); ibin++)
    {
      fprintf(fp,"%6.2f %11.2f %19.5e %15.5e\n",
              log10en[ibin],
              log10en_bsize[ibin],
              nevents[ibin], 
              exposure[ibin]);
    }
  fflush(fp);
  fclose(fp);
}

