// Dmitri Ivanov <dmiivanov@gmail.com>

#include <cstdio>
#include <cstdlib>
#include <vector>
#include "specfit_uti.h"
#include "TF1.h"
#include "TGraph.h"
#include "TROOT.h"
#include "TFeldmanCousins.h"
#include "TMath.h"
#if __cplusplus >= 201103L
#include <cstdint>
#else
#include <stdint.h>
#endif

using namespace std;

#ifndef SIZE_MAX
#define SIZE_MAX ((size_t)(-1))
#endif

// To get a unique name for an object
TString specfit_uti::get_unique_object_name(const char *basename)
{
  for (size_t i = 0; i < (size_t) -1; i++)
    {
      TString name = (TString(basename) + "_") + TString::Format("%zu", i);
      // skip if already present
      if(gROOT->FindObject(name))
	continue;
      return name;
    }
  fprintf(stderr, "warning: ran out of object names for prefix %s\n", basename);
  return TString("");
}

// to get Feldman-Cousins error bars
std::pair<Double_t, Double_t> specfit_uti::get_fc_errors(Double_t n)
{
  static TGraph *gFCerr[2] =
  { 0 };
  if(!gFCerr[0])
    {
      gFCerr[0] = new TGraph(20);
      gFCerr[1] = new TGraph(20);
      TFeldmanCousins fldc_err(0.683);
      for (Int_t i = 0; i < gFCerr[0]->GetN(); i++)
	{
	  gFCerr[0]->SetPoint(i, (Double_t) i, (Double_t) i - fldc_err.CalculateLowerLimit((Double_t) i, 0.0));
	  gFCerr[1]->SetPoint(i, (Double_t) i, fldc_err.CalculateUpperLimit((Double_t) i, 0.0) - (Double_t) i);
	}
    }
  if(n < (Double_t) gFCerr[0]->GetN())
    return std::pair<Double_t, Double_t>(gFCerr[0]->Eval(n), gFCerr[1]->Eval(n));
  return std::pair<Double_t, Double_t>(TMath::Sqrt(n), TMath::Sqrt(n));
}

// to get the lower error bar
Double_t specfit_uti::get_fc_error_low(Double_t n)
{
  return get_fc_errors(n).first;
}

// to get the upper error bar
Double_t specfit_uti::get_fc_error_high(Double_t n)
{
  return get_fc_errors(n).second;
}

// linear bin size
Double_t specfit_uti::GetLinBinSize(Double_t log10en, Double_t log10en_bsize)
{
  return TMath::Power(10.0, log10en + log10en_bsize / 2.0) - TMath::Power(10.0, log10en - log10en_bsize / 2.0);
}

// get the significance in sigma units
// from the chance probability
Double_t specfit_uti::pchance2sigma(Double_t pchance, Bool_t pwarning)
{
  if(pchance <= 0.0 || pchance > 1.0)
    {
      if(pwarning)
	{
	  fprintf(stderr, "warning: pchance2sigma: ");
	  fprintf(stderr, "invalid chance probability, ");
	  fprintf(stderr, "must be in (0.0-1.0] range\n");
	}
      return 0.0;
    }
  if(pchance >= 0.5)
    return 0.0;
  return sqrt(2.0) * TMath::ErfcInverse(2.0 * pchance);
}

// express chance probability in sigma units
Double_t specfit_uti::Sigma2Pchance(Double_t pchange_in_sigma)
{
  return 0.5 * (1 - TMath::Erf(pchange_in_sigma / sqrt(2.0)));
}

// get the chance probability of a Poisson fluctuation
Double_t specfit_uti::PoissonPchance(Int_t nobserved, Double_t nexpected, Bool_t in_sigma_units)
{
  Double_t pchance = 0.0;
  if((Double_t) nobserved <= nexpected)
    {
      for (Int_t i = nobserved; i >= 0; i--)
	{
	  Double_t d_pchance = TMath::PoissonI(i, nexpected);
	  pchance += d_pchance;
	  if(d_pchance < 1e-300)
	    break;
	}
    }
  else
    {
      for (Int_t i = nobserved; i <= 0x7FFFFFFF; i++)
	{
	  Double_t d_pchance = TMath::PoissonI(i, nexpected);
	  pchance += d_pchance;
	  if(d_pchance < 1e-300)
	    break;
	}
    }
  return (in_sigma_units ? pchance2sigma(pchance) : pchance);
}

// to obtain E^{3}J function from J if J was constructed using formula
TF1* specfit_uti::get_e3j_from_j(TF1 *f_J)
{
  if(!f_J)
    return 0;
  TString frm = f_J->GetExpFormula();
  for (Int_t i = 0; i < f_J->GetNpar(); i++)
    (void) (frm.ReplaceAll(TString::Format("[%s]", f_J->GetParName(i)), TString::Format("[%d]", i)));
  if(!frm.Length())
    return 0;
  frm = TString("10^(3.0*x)*") + frm;
  TString name = f_J->GetName();
  if(name.BeginsWith("fJ"))
    {
      name.Remove(0, 2);
      name = TString("fE3J") + name;
    }
  else
    name += "_E3";
  TF1 *f = new TF1(specfit_uti::get_unique_object_name(name), frm, f_J->GetXmin(), f_J->GetXmax());
  for (Int_t i = 0; i < f_J->GetNpar(); i++)
    {
      f->SetParName(i, f_J->GetParName(i));
      f->SetParameter(i, f_J->GetParameter(i));
      f->SetParError(i, f_J->GetParError(i));
    }
  f->SetTitle(";log_{10}(E/eV);E^{3}J");
  f->SetLineStyle(f_J->GetLineStyle());
  f->SetLineColor(f_J->GetLineColor());
  return f;
}

NamespaceImp(specfit_uti);
