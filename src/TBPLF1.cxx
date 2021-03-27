#include "TBPLF1.h"
#include "specfit_uti.h"
#include "TMath.h"

ClassImp(TBPLF1);


// make flux formulas based on the function type and how many break points
TString TBPLF1::make_formula(Int_t nbreaks, const char *ftype, Double_t scalefactor, Double_t log10en_min)
{
  if(nbreaks < 0)
    nbreaks = 0;
  const TString xdiff = TString::Format("(x-%f)", log10en_min);
  TString s_ftype = TString(ftype);
  s_ftype.ToUpper();
  TString frm = "";
  if(s_ftype == "J" || s_ftype == "E3J" || s_ftype == "EJ")
    {
      frm = TString::Format("%e*[0]*", scalefactor);
      if(s_ftype == "E3J")
	frm += "10^(3.0*x)*";
      if(s_ftype == "EJ")
	frm += "10^(x)*";
      frm += "(";
      if(nbreaks == 0)
	frm += TString::Format("10^([1]*%s)", xdiff.Data());
      else
	{
	  // power law before the first break
	  frm += TString::Format("(x<[%d])*10^([1]*%s)", 2 + nbreaks, xdiff.Data());
	  // power laws between the breaks
	  for (Int_t ibreak = 0; ibreak < nbreaks - 1; ibreak++)
	    {
	      Int_t ipar_p_after_break = ibreak + 2;
	      Int_t ipar_break = nbreaks + 2 + ibreak;
	      Int_t ipar_next_break = nbreaks + 3 + ibreak;
	      frm += TString::Format("+([%d]<=x)*(x<[%d])*10^(%s+[%d]*%s)", ipar_break, ipar_next_break, get_pcf(ibreak, nbreaks, log10en_min).Data(), ipar_p_after_break,
		  xdiff.Data());
	    }
	  // power law after the last break
	  frm += TString::Format("+([%d]<=x)*10^(%s+[%d]*%s)", 2 * nbreaks + 1, get_pcf(nbreaks - 1, nbreaks, log10en_min).Data(), nbreaks + 1, xdiff.Data());
	}
      frm += ")";
    }
  else if(s_ftype == "J>" || s_ftype == "E2J>")
    {
      frm = TString::Format("%e*[0]*", scalefactor);
      if(s_ftype == "E2J>")
	frm += "10^(2.0*x)*";
      frm += TString::Format("10^(%f)*", log10en_min);
      frm += "(";
      if(nbreaks == 0)
	frm += TString::Format("-(x<%f)/(1+[1])-(%f<=x)*10^((1+[1])*%s)/(1+[1])", log10en_min, log10en_min, xdiff.Data());
      else
	{
	  TString lbdiff = TString::Format("([%d]-%f)", nbreaks + 2, log10en_min);
	  frm += TString::Format("(x<=[%d])*((10^((1+[1])*%s)-10^((1+[1])*%s))/(1+[1])", nbreaks + 2, lbdiff.Data(), xdiff.Data());
	  for (Int_t jbreak = 0; jbreak < nbreaks; jbreak++)
	    {
	      TString jp2 = TString::Format("(1.0+[%d])", 2 + jbreak);
	      TString jlb1 = TString::Format("[%d]", nbreaks + 2 + jbreak);
	      TString jlb1diff = TString::Format("([%d]-%f)", nbreaks + 2 + jbreak, log10en_min);
	      TString jpcf = get_pcf(jbreak, nbreaks, log10en_min);
	      if(jbreak < nbreaks - 1)
		{
		  TString jlb2diff = TString::Format("([%d]-%f)", nbreaks + 3 + jbreak, log10en_min);
		  frm += TString::Format("+10^(%s)*(10^(%s*%s)-10^(%s*%s))/%s", jpcf.Data(), jp2.Data(), jlb2diff.Data(), jp2.Data(), jlb1diff.Data(), jp2.Data());
		}
	      else
		frm += TString::Format("-10^(%s+%s*%s)/%s", jpcf.Data(), jp2.Data(), jlb1diff.Data(), jp2.Data());
	    }
	  frm += ")";
	  for (Int_t ibreak = 0; ibreak < nbreaks; ibreak++)
	    {
	      TString ip2 = TString::Format("(1+[%d])", 2 + ibreak);
	      TString ilb1 = TString::Format("[%d]", nbreaks + 2 + ibreak);
	      TString ilb1diff = TString::Format("([%d]-%f)", nbreaks + 2 + ibreak, log10en_min);
	      TString ipcf = get_pcf(ibreak, nbreaks, log10en_min);
	      if(ibreak < nbreaks - 1)
		{
		  TString ilb2 = TString::Format("[%d]", nbreaks + 2 + ibreak + 1);
		  TString ilb2diff = TString::Format("([%d]-%f)", nbreaks + 2 + ibreak + 1, log10en_min);
		  frm += TString::Format("+(%s<x)*(x<=%s)*(", ilb1.Data(), ilb2.Data());
		  frm += TString::Format("10^(%s)*(10^(%s*%s)-10^(%s*%s))/%s", ipcf.Data(), ip2.Data(), ilb2diff.Data(), ip2.Data(), xdiff.Data(), ip2.Data());
		  for (Int_t jbreak = ibreak + 1; jbreak < nbreaks; jbreak++)
		    {
		      TString jp1 = TString::Format("(1+[%d])", 2 + jbreak);
		      TString jlb1diff = TString::Format("([%d]-%f)", nbreaks + 2 + jbreak, log10en_min);
		      TString jpcf = get_pcf(jbreak, nbreaks, log10en_min);
		      if(jbreak < nbreaks - 1)
			{
			  TString jlb2diff = TString::Format("([%d]-%f)", nbreaks + 2 + jbreak + 1, log10en_min);
			  frm += TString::Format("+10^(%s)*(10^(%s*%s)-10^(%s*%s))/%s", jpcf.Data(), jp1.Data(), jlb2diff.Data(), jp1.Data(), jlb1diff.Data(), jp1.Data());
			}
		      else
			frm += TString::Format("-10^(%s+%s*%s)/%s", jpcf.Data(), jp1.Data(), jlb1diff.Data(), jp1.Data());
		    }
		  frm += ")";
		}
	      else
		frm += TString::Format("-(%s<x)*10^(%s+%s*%s)/%s", ilb1.Data(), ipcf.Data(), ip2.Data(), xdiff.Data(), ip2.Data());
	    }
	}
      frm += ")";
    }
  else
    {
      std::cerr << "ERROR: make_formula: function type '" << ftype << "' not understood; use 'J', 'J>', 'EJ', 'E3J', or 'E2J>'" << std::endl;
      return TString("");
    }
  return frm;
}

// get the power coefficient that applies after the break point
// this takes into account all the prior break points
TString TBPLF1::get_pcf(Int_t ibreak, Int_t nbreaaks, Double_t log10en_min)
{
  TString frm = "";
  if(ibreak < 0 || ibreak > nbreaaks - 1 || nbreaaks < 1)
    {
	frm = "0.0";
	return frm;
    }
  for (Int_t jbreak = 0; jbreak <= ibreak; jbreak++)
    {
	Int_t jpar_p_before_break = jbreak + 1;
	Int_t jpar_p_after_break = jbreak + 2;
	Int_t jpar_break = nbreaaks + 2 + jbreak;
	if(frm.Length())
	  frm += "+";
	frm += TString::Format("([%d]-[%d])*([%d]-%f)", jpar_p_before_break, jpar_p_after_break, jpar_break, log10en_min);
    }
  return frm;
}


Double_t TBPLF1::MultiplyAndIntegrate_dE(const TF1* f, Double_t log10en_start, Double_t log10en_end, Double_t esprel) const
{
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
#define _Integral_(_a_,_b_) Integral(_a_,_b_,esprel)
#else
#define _Integral_(_a_,_b_) Integral(_a_,_b_,(const Double_t*)0,esprel)
#endif
  if(log10en_end < log10en_start)
    {
      fprintf(stderr,"ERROR: TBPLF1::Integrate_dE: upper integration limit is smaller than the lower!\n");
      return 0.0;
    }
  TSPECFITF1* fdE = new TSPECFITF1(specfit_uti::get_unique_object_name("f_dE"),"10^(x) * TMath::Ln10()", log10en_start, log10en_end);
  TSPECFITF1* fJoined = (f ? Multiply(specfit_uti::get_unique_object_name(TString(GetName())+TString(f->GetName())),this,f) : (TSPECFITF1*)Clone(specfit_uti::get_unique_object_name(TString(GetName()))));
  TSPECFITF1* fIntegrand0 = Multiply(specfit_uti::get_unique_object_name(TString(GetName())+"_dE"),fJoined,fdE);
  Double_t A = TMath::Abs(fIntegrand0->Eval(log10en_start));
  if(A < 1e-323)
    {
      fprintf(stderr,"ERROR: TBPLF1::Integrate_dE: integrand value is too small at log10en_start = %f for double precision!\n",log10en_start);
      return 0.0;
    }
  TSPECFITF1* fIntegrand = Add(specfit_uti::get_unique_object_name(TString(GetName())+"_dE_scaled"),fIntegrand0,1/A);
  delete fdE;
  delete fJoined;
  delete fIntegrand0;
  Double_t result = 0.0; // result of the integration
  std::vector<Double_t> xd_breaks; // break points that have been crossed during the integration
  for (Int_t ibreak = 0; ibreak < GetNbreaks(); ibreak++)
    {
      if(log10en_start < GetParameter(GetNbreaks() + 2 + ibreak) && log10en_end > GetParameter(GetNbreaks() + 2 + ibreak))
	xd_breaks.push_back(GetParameter(GetNbreaks() + 2 + ibreak));
    }
  if(xd_breaks.empty())
    result += fIntegrand->_Integral_(log10en_start, log10en_end);
  else
    {
      result += fIntegrand->_Integral_(log10en_start, xd_breaks[0]);
      for (Int_t i = 0; i < (Int_t) xd_breaks.size() - 1; i++)
	result += fIntegrand->_Integral_(xd_breaks[i], xd_breaks[i + 1]);
      result += fIntegrand->_Integral_(xd_breaks[xd_breaks.size() - 1], log10en_end);
    }
  result *= A;
#undef _Integral_
  delete fIntegrand; // clean up
  return result; // return the answer
}

// title according to the type of the function
void TBPLF1::set_default_title(const char *ftype)
{
  TString s_ftype(ftype);
  if(s_ftype == "J")
    SetTitle(";log_{10}(E/eV);J");
  if(s_ftype == "E3J")
    SetTitle(";log_{10}(E/eV);E^{3}J");
  if(s_ftype == "EJ")
    SetTitle(";log_{10}(E/eV);EJ");
  if(s_ftype == "J>")
    SetTitle(";log_{10}(E/eV);J_{>}");
  if(s_ftype == "E2J>")
    SetTitle(";log_{10}(E/eV);E^{2}J_{>}");
}
