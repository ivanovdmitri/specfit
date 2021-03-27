// Dmitri Ivanov <dmiivanov@gmail.com>

#include <cstdio>
#include <cstdlib>
#include <vector>
#include "specfit_uti.h"
#include "TF1.h"
#include "TAxis.h"
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
      if(gROOT->FindObject(name)) // @suppress("Function cannot be resolved") // @suppress("Method cannot be resolved")
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

void specfit_uti::interpolate_flux_point(const TGraphErrors* g, Double_t log10en, Double_t &flux, Double_t &eflux)
{
  // to interpolate flux

  flux = 0; // interpolated flux
  eflux = 0; // interpolating error bars
  if(!g)
    {
      fprintf(stderr, "error: interpolate_flux_point: TGraphErrors *g pointer is NULL\n");
      return;
    }
  if(g->GetN() < 1)
    return;
  Double_t *x = g->GetX(); // array of values of energy
  Double_t *y = g->GetY(); // array of values of flux
  if(g->GetN() == 1)
    {
      flux = y[0];
      eflux = g->GetErrorY(0);
      return;
    }
  // points from which we do interpolation
  Double_t log10en1 = 0.0, log10en2 = 0.0;
  Double_t flux1 = 0.0, eflux1 = 0.0, flux2 = 0.0, eflux2 = 0.0;
  // Determine the flux start and stop points (flux must be non-zero)
  // non-zero means that the flux should be more than 1E-30 times
  // the average value of the flux in the range under consideration
  Double_t flux_average = 0.0;
  for (Int_t i = 0; i < g->GetN(); i++)
    flux_average += y[i];
  flux_average /= (Double_t) g->GetN();
  Int_t i_start = g->GetN() - 1;
  Int_t i_stop = 0;
  for (Int_t i = 0; i < g->GetN(); i++)
    {
      // If we landed on an actual data point then use that data point
      // provided that the flux is not zero there
      if(TMath::Abs(x[i] - log10en) < 1e-12 && TMath::Abs(y[i]) > 1e-30 * flux_average)
	{
	  flux = y[i];
	  eflux = g->GetErrorY(i);
	  return;
	}
      if(TMath::Abs(y[i]) > 1e-30 * flux_average)
	{
	  if(i < i_start)
	    i_start = i;
	  if(i > i_stop)
	    i_stop = i;
	}
    }
  // if there's only one non-zero flux point then
  // we can only use that
  if(i_start == i_stop)
    {
      flux = y[i_start];
      eflux = g->GetErrorY(i_start);
      return;
    }
  // Determine the closest lower and upper energy bounds, and
  // so that flux is non-zero
  Int_t i1 = i_start;
  Int_t i2 = i_stop;
  for (Int_t i = i_start; i <= i_stop; i++)
    {
      if(TMath::Abs(y[i]) < 1e-30 * flux_average)
	continue; // use only non-zero flux points
      // closest lower energy bound
      if(x[i] < log10en)
	{
	  if(x[i] > x[i1])
	    i1 = i;
	}
      // closest upper energy bound
      if(x[i] > log10en)
	{
	  if(x[i] < x[i2])
	    i2 = i;
	}
    }
  // if lower bound correspond to the upper end point of the flux
  // then use it as the upper bound and find the next
  // closest non-zero lower bound
  if(i1 == i_stop)
    {
      i2 = i_stop;
      for (Int_t i = i_stop - 1; i >= i_start; i--)
	{
	  if(TMath::Abs(y[i]) < 1e-30 * flux_average)
	    continue;
	  i1 = i;
	  break;
	}
    }
  // if the upper bound corresponds to the lower end point of the
  // flux then use it as the lower bound and find the next
  // closest non-zero upper bound
  if(i2 == i_start)
    {
      i1 = i_start;
      for (Int_t i = i_start + 1; i <= i_stop; i++)
	{
	  if(TMath::Abs(y[i]) < 1e-30 * flux_average)
	    continue;
	  i2 = i;
	  break;
	}
    }
  // do logarithmic interpolation
  g->GetPoint(i1, log10en1, flux1);
  g->GetPoint(i2, log10en2, flux2);
  flux = exp(log(flux1) + (log(flux2) - log(flux1)) / (log10en2 - log10en1) * (log10en - log10en1));
  // lower error bar obtained by interpolating the lower error bars
  eflux1 = g->GetErrorY(i1);
  eflux2 = g->GetErrorY(i2);
  eflux = flux / TMath::Abs(log10en2 - log10en1)
      * sqrt(eflux1 * eflux1 / flux1 / flux1 * (log10en2 - log10en) * (log10en2 - log10en) + eflux2 * eflux2 / flux2 / flux2 * (log10en - log10en1) * (log10en - log10en1));
}

// 1st divided by the 2nd
TGraphErrors* specfit_uti::FluxRatio(const TGraphErrors* flux1, const TGraphErrors* flux2, Bool_t ok_to_extrapolate)
{
  TGraphErrors *gRat = new TGraphErrors(0);
  gRat->SetMarkerStyle(21);
  gRat->SetLineWidth(3);
  if(!flux1)
    {
      fprintf(stderr, "error: TGraphErrors *flux1 is NULL!\n");
      return gRat;
    }
  if(!flux2)
    {
      fprintf(stderr, "error: TGraphErrors *flux2 is NULL!\n");
      return gRat;
    }
  Int_t irat_point = 0;
  for (Int_t i = 0; i < flux1->GetN(); i++)
    {
      Double_t log10en = 0, f1 = 0;
      flux1->GetPoint(i, log10en, f1);
      if(!ok_to_extrapolate)
	{
	  if(log10en > flux2->GetX()[flux1->GetN() - 1])
	    continue;
	  if(log10en < flux2->GetX()[0])
	    continue;
	}
      Double_t ef1 = flux1->GetErrorY(i);
      Double_t f2 = 0, ef2 = 0;
      interpolate_flux_point(flux2, log10en, f2, ef2);
      if(TMath::Abs(f1) < 1e-30 * TMath::Abs(flux1->GetY()[flux1->GetN() / 2]))
	continue;
      if(TMath::Abs(f2) < 1e-30 * TMath::Abs(flux2->GetY()[flux2->GetN() / 2]))
	continue;
      Double_t r = f1 / f2;
      Double_t er = r * sqrt(ef1 * ef1 / f1 / f1 + ef2 * ef2 / f2 / f2);
      gRat->SetPoint(irat_point, log10en, r);
      gRat->SetPointError(irat_point, 0, er);
      irat_point++;
    }
  gRat->GetXaxis()->CenterTitle();
  gRat->GetYaxis()->CenterTitle();
  gRat->GetXaxis()->SetTitleSize(0.055);
  gRat->GetYaxis()->SetTitleSize(0.055);
  gRat->SetMarkerStyle(24);
  gRat->SetMarkerSize(1.5);
  return gRat;
}

TGraphErrors* specfit_uti::FluxRatio_Energy_Bins(const TGraphErrors* flux1, const TGraphErrors* flux2, Int_t nebins, Double_t log10en_lo, Double_t log10en_up, Bool_t ok_to_extrapolate)
{
  // Ratio of the fluxes using binning

  TGraphErrors *gRat = new TGraphErrors(0);
  gRat->SetMarkerStyle(21);
  gRat->SetLineWidth(3);
  if(!flux1)
    {
      fprintf(stderr, "error: TGraphErrors *flux1 is NULL!\n");
      return gRat;
    }
  if(!flux2)
    {
      fprintf(stderr, "error: TGraphErrors *flux2 is NULL!\n");
      return gRat;
    }
  Int_t irat_point = 0;
  Double_t log10en_bsize = (log10en_up - log10en_lo) / (Double_t) nebins;

  Double_t log10en2_nonzero_lo = flux2->GetX()[flux2->GetN() - 1];
  Double_t log10en2_nonzero_up = flux2->GetX()[0];
  for (Int_t i = 0; i < flux2->GetN(); i++)
    {
      if(TMath::Abs(flux2->GetY()[flux2->GetN() - 1 - i]) > 1e-30 * TMath::Abs(flux2->GetY()[flux2->GetN() / 2]))
	log10en2_nonzero_lo = flux2->GetX()[flux2->GetN() - 1 - i] - 1e-3;
      if(TMath::Abs(flux2->GetY()[i]) > 1e-30 * TMath::Abs(flux2->GetY()[flux2->GetN() / 2]))
 	log10en2_nonzero_up = flux2->GetX()[i] + 1e-3;
    }
  for (Int_t i = 0; i < nebins; i++)
    {
      Double_t log10en = log10en_lo + log10en_bsize * ((Double_t) i + 0.5);
      Double_t f1 = 0, ef1 = 0;
      interpolate_flux_point(flux1, log10en, f1, ef1);
      Double_t f2 = 0, ef2 = 0;
      interpolate_flux_point(flux2, log10en, f2, ef2);
      // don't compute the ratio after the first flux ends in all cases
      if(log10en > flux1->GetX()[flux1->GetN() - 1])
	continue;
      if(log10en < flux1->GetX()[0])
	continue;
      if(!ok_to_extrapolate)
	{
	  if(log10en > log10en2_nonzero_up)
	    continue;
	  if(log10en < log10en2_nonzero_lo)
	    continue;
	}
      if(TMath::Abs(f1) < 1e-30 * TMath::Abs(flux1->GetY()[flux1->GetN() / 2]))
	continue;
      if(TMath::Abs(f2) < 1e-30 * TMath::Abs(flux2->GetY()[flux2->GetN() / 2]))
	continue;
      Double_t r = f1 / f2;
      Double_t er = r * sqrt(ef1 * ef1 / f1 / f1 + ef2 * ef2 / f2 / f2);
      gRat->SetPoint(irat_point, log10en, r);
      gRat->SetPointError(irat_point, 0, er);
      irat_point++;
    }
  gRat->GetXaxis()->CenterTitle();
  gRat->GetYaxis()->CenterTitle();
  gRat->GetXaxis()->SetTitleSize(0.055);
  gRat->GetYaxis()->SetTitleSize(0.055);
  gRat->SetMarkerStyle(24);
  gRat->SetMarkerSize(1.5);
  return gRat;
}

NamespaceImp(specfit_uti);
