#include "TBPLF1.h"

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
      frm = TString::Format("[0]*%e*", scalefactor);
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
      frm = TString::Format("[0]*%e*", scalefactor);
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
