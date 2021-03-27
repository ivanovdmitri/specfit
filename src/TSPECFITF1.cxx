#include "TSPECFITF1.h"

ClassImp(TSPECFITF1);

TSPECFITF1::~TSPECFITF1()
{
  ;
}


void TSPECFITF1::SetParLimitsCS(const char *cpsparlimits)
{
  std::vector<std::pair<Double_t, Double_t> > pairs = TSPECFITF1::parse_cs_doubles_pairs(cpsparlimits);
  if((Int_t)pairs.size() != GetNpar())
    {
      fprintf(stderr,"error: number of parsed limit pairs %d found in '%s' is not the same as the number of function parameters %d!\n",
	  (Int_t)pairs.size(),cpsparlimits,GetNpar());
      return;
    }
  for (Int_t ipar = 0; ipar < GetNpar(); ipar++)
      SetParLimits(ipar, pairs[ipar].first, pairs[ipar].second);
}

TString TSPECFITF1::GetParLimitsCS() const
{
  TString result = "";
  for (Int_t ipar = 0; ipar < GetNpar(); ipar++)
    {
	Double_t parmin = 0, parmax = 0;
	GetParLimits(ipar, parmin, parmax);
	if(result.Length())
	  result += ",";
	result += TString::Format("(%.9e,%.9e)", parmin, parmax);
    }
  return result;
}

TString TSPECFITF1::cs_doubles(Int_t n, const Double_t *data)
{
  TString result = "";
  for (Int_t i = 0; i < n; i++)
    {
      if(result.Length())
	result += ",";
      result += TString::Format("%.9e", data[i]);
    }
  return result;
}

TString TSPECFITF1::cs_strings(Int_t n, const TString *data)
{
  TString result = "";
  for (Int_t i = 0; i < n; i++)
    {
      if(result.Length())
	result += ",";
      result += data[i];
    }
  return result;
}

// parse comma-separated array of doubles into a vector of doubles
std::vector<std::pair<Double_t,Double_t> > TSPECFITF1::parse_cs_doubles_pairs(const char *cs_doubles_pairs)
{
  std::vector<std::pair<Double_t,Double_t> > result;
  if(!cs_doubles_pairs)
    return result;
  Ssiz_t from = 0;
  TString tok = "";
  TString s_cs_doubles_pairs = cs_doubles_pairs;
  s_cs_doubles_pairs.ReplaceAll(" ", "");
  while (s_cs_doubles_pairs.Tokenize(tok, from, "),"))
    {
      if (tok.Length() < 3)
	{
	  fprintf(stderr,"error: TSPECFITF1::parse_cs_doubles_pairs: failed to parse %s\n", cs_doubles_pairs);
	  result.clear();
	  return result;
	}

      if(tok[tok.Length() - 1] == ',')
	tok.Remove(tok.Length() - 1);
      if(tok[0] == '(')
	tok.Remove(0, 1);
      if(tok[tok.Length() - 1] == ')')
	tok.Remove(tok.Length() - 1);

      std::vector<Double_t> vp = parse_cs_doubles(tok.Data());
      if (vp.size() != 2)
       	{
       	  fprintf(stderr,"error: TSPECFITF1::parse_cs_doubles_pairs: failed to parse %s\n", cs_doubles_pairs);
       	  result.clear();
       	  return result;
       	}
      std::pair<Double_t,Double_t> p = std::make_pair(vp[0],vp[1]);
      result.push_back(p);
    }
  return result;
}

// parse comma-separated array of doubles into a vector of doubles
std::vector<Double_t> TSPECFITF1::parse_cs_doubles(const char *cs_doubles)
{
  std::vector<Double_t> result;
  if(!cs_doubles)
    return result;
  Ssiz_t from = 0;
  TString tok = "";
  TString s_cs_doubles = cs_doubles;
  s_cs_doubles.ReplaceAll(" ", "");
  while (s_cs_doubles.Tokenize(tok, from, ","))
    {
      if(tok.IsFloat())
	result.push_back(tok.Atof());
    }
  return result;
}

// parse comma-separated strings into a vector of strings
std::vector<TString> TSPECFITF1::parse_cs_strings(const char *cs_strings)
{
  std::vector<TString> result;
  if(!cs_strings)
    return result;
  Ssiz_t from = 0;
  TString tok = "";
  TString s_cs_strings = cs_strings;
  s_cs_strings.ReplaceAll(" ", "");
  while (s_cs_strings.Tokenize(tok, from, ","))
    {
      if(!tok.IsWhitespace())
	result.push_back(tok);
    }
  return result;
}

// in addition, add a method to set parameters from the array;
// a method with these arguments is missing
// in the TF1 class.
void TSPECFITF1::SetParNames(const TString *parnames)
{
  if(!parnames)
    return;
  for (Int_t i = 0; i < GetNpar(); i++)
    SetParName(i, parnames[i]);
}

void TSPECFITF1::SetParNamesCS(const char *csparnames)
{
  if(!csparnames)
    return;
  std::vector<TString> parnames = parse_cs_strings(csparnames);
  if((Int_t) parnames.size() == GetNpar())
    {
      if((Int_t) parnames.size())
	SetParNames(&parnames[0]);
    }
  else
    {
      fprintf(stderr, "ERROR: size of your comma separated parameter_names (%d) not equal to the number of parameters (%d); ", (Int_t) parnames.size(), GetNpar());
      fprintf(stderr, "not using par_names='%s'!\n", csparnames);
    }
}

void TSPECFITF1::SetParametersCS(const char *csparameters)
{
  if(!csparameters)
    return;
  std::vector<Double_t> parameters = parse_cs_doubles(csparameters);
  if((Int_t) parameters.size() == GetNpar())
    {
      if((Int_t) parameters.size())
	SetParameters(&parameters[0]);
    }
  else
    {
      fprintf(stderr, "ERROR: size of your comma separated parameters (%d) not equal to the number of parameters (%d); ", (Int_t) parameters.size(), GetNpar());
      fprintf(stderr, "not using parameters='%s'!\n", csparameters);
    }
}

void TSPECFITF1::SetParErrorsCS(const char *csparerrors)
{
  if(!csparerrors)
    return;
  std::vector<Double_t> parerrors = parse_cs_doubles(csparerrors);
  if((Int_t) parerrors.size() == GetNpar())
    {
      if((Int_t) parerrors.size())
	SetParErrors(&parerrors[0]);
    }
  else
    {
      fprintf(stderr, "ERROR: size of your comma separated parameter_errors (%d) not equal to the number of parameters (%d); ", (Int_t) parerrors.size(), GetNpar());
      fprintf(stderr, "not using parameter_errors='%s'!\n", csparerrors);
    }
}

void TSPECFITF1::GetParNames(TString *parnames) const
{
  for (Int_t i = 0; i < GetNpar(); i++)
    parnames[i] = GetParName(i);
}

 std::vector<TString> TSPECFITF1::GetParNames() const
{
  std::vector<TString> parnames = std::vector<TString>(GetNpar());
  GetParNames(&parnames[0]);
  return parnames;
}

TSPECFITF1* TSPECFITF1::Add(const char *newname, const char *frm_start, const TF1 *f1, const TF1 *f2, Double_t c1, Double_t c2)
{
  TF1 *f0 = new TF1(TString(f1->GetName()) + TString(f2->GetName()) + "_tmp", frm_start);
  TString frm_0 = GetExpFormula(f0, 0);
  TString frm_1 = GetExpFormula(f1, f0->GetNpar());
  TString frm_2 = GetExpFormula(f2, f0->GetNpar() + f1->GetNpar());
  TString frm = TString::Format("(%s + %.9e * (%s) + %.9e * (%s))", frm_0.Data(), c1, frm_1.Data(), c2, frm_2.Data());
  std::vector<TString> parnames;
  std::vector<Double_t> params;
  std::vector<Double_t> parerrors;
// for f0, put in the default values that occur after initialization
  for (Int_t i = 0; i < f0->GetNpar(); i++)
    {
      parnames.push_back(f0->GetParName(i));
      params.push_back(f0->GetParameter(i));
      parerrors.push_back(f0->GetParError(i));
    }
  delete f0; // clean up
// put in whatever values were stored in f1, f2
  for (Int_t i = 0; i < f1->GetNpar(); i++)
    {
      parnames.push_back(f1->GetParName(i));
      params.push_back(f1->GetParameter(i));
      parerrors.push_back(f1->GetParError(i));
    }
  for (Int_t i = 0; i < f2->GetNpar(); i++)
    {
      parnames.push_back(f2->GetParName(i));
      params.push_back(f2->GetParameter(i));
      parerrors.push_back(f2->GetParError(i));
    }

// domain of the function is determined by f1, f2 (f0 was used
// only for manipulating the formula
  Double_t xmin = TMath::Min(f1->GetXmin(), f2->GetXmin());
  Double_t xmax = TMath::Max(f1->GetXmax(), f2->GetXmax());
  TSPECFITF1 *f = new TSPECFITF1(newname, frm.Data(), xmin, xmax, &parnames[0], &params[0], &parerrors[0]);
  return f;
}

TSPECFITF1* TSPECFITF1::Add(const char *newname, const char *frm_start, const TF1 *f, Double_t c)
{
  TF1 *f2 = new TF1(TString(f->GetName()) + "_tmp_2", "0.0");
  TSPECFITF1 *ff = Add(newname, frm_start, f, f2, c, 1.0);
  delete f2;
  return ff;
}


TString TSPECFITF1::GetExpFormula(const TF1 *f, Int_t n_offset)
{
  TString frm = f->GetExpFormula();
  for (Int_t i = 0; i < f->GetNpar(); i++)
    {
      TString named_param = TString::Format("[%s]", f->GetParName(i));
      TString intno_param = TString::Format("[%d]", i);
      frm.ReplaceAll(named_param, intno_param); // first cast formula into a more basic form (used in older ROOT by default)
    }
  for (Int_t i = f->GetNpar() - 1; i >= 0; i--)
    {
      TString intno_param = TString::Format("[%d]", i);
      TString intno_param_plus_offset = TString::Format("[%d]", i + n_offset);
      frm.ReplaceAll(intno_param, intno_param_plus_offset); // apply offset, if requested
    }
  return frm;
}


TSPECFITF1* TSPECFITF1::Multiply(const char *newname, const TF1 *f1, const TF1 *f2)
{
  TString frm_1 = GetExpFormula(f1, 0);
  TString frm_2 = GetExpFormula(f2, f1->GetNpar());
  TString frm = TString::Format("(%s) * (%s)", frm_1.Data(), frm_2.Data());
  std::vector<TString> parnames;
  std::vector<Double_t> params;
  std::vector<Double_t> parerrors;
  // put in whatever values were stored in f1, f2
  for (Int_t i = 0; i < f1->GetNpar(); i++)
    {
      parnames.push_back(f1->GetParName(i));
      params.push_back(f1->GetParameter(i));
      parerrors.push_back(f1->GetParError(i));
    }
  for (Int_t i = 0; i < f2->GetNpar(); i++)
    {
      parnames.push_back(f2->GetParName(i));
      params.push_back(f2->GetParameter(i));
      parerrors.push_back(f2->GetParError(i));
    }
  // domain of the function is determined by f1, f2
  Double_t xmin = TMath::Min(f1->GetXmin(), f2->GetXmin());
  Double_t xmax = TMath::Max(f1->GetXmax(), f2->GetXmax());
  TSPECFITF1 *f = new TSPECFITF1(newname, frm.Data(), xmin, xmax, &parnames[0], &params[0], &parerrors[0]);
  return f;
}

void TSPECFITF1::Scale(Double_t c)
{
  std::vector<TString>  parnames = GetParNames();
  std::vector<Double_t> parameters(GetParameters(), GetParameters() + GetNpar());
  std::vector<Double_t> parerrors(GetParErrors(), GetParErrors() + GetNpar());
  TString frm_0 = GetExpFormula(this, 0);
  TString frm = TString::Format("%e * (%s)", c, frm_0.Data());
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
  GetFormula()->SetTitle(frm); // In ROOT6 and above TFormula is embedded
  GetFormula()->Compile(frm);
#else
  SetTitle(frm); // Before ROOT6, TF1 inherited from TFormula
  Compile(frm);
#endif
  SetParNames(&parnames.front());
  SetParameters(&parameters.front());
  SetParErrors(&parerrors.front());
}
