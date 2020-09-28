{
  THtml h;
  TString SPECFIT= (gSystem->Getenv("SPECFIT") ? 
                    TString(gSystem->Getenv("SPECFIT")) : 
                    TString(gSystem->DirName(__FILE__)));
  gROOT->Macro(SPECFIT+TString("/load_specfit_cpplib.C"));
  h.SetOutputDir(SPECFIT+"/htmldoc");
  h.SetInputDir(SPECFIT+TString("/inc")+
		TString(":")+
		SPECFIT+TString("/src")+
		TString(":")+
		TString(gSystem->Getenv("ROOTSYS"))+
		TString(":")+
		TString(gSystem->Getenv("ROOTSYS"))+"/include");
  h.SetProductName("SPECFIT");
  h.MakeIndex("TCRFlux*|TSPECFITF1|TBPLF1|specfit_*");
  h.MakeAll(1,"TCRFlux*|TSPECFITF1|TBPLF1|specfit_*");
  h.MakeTree("TCRFlux");
  h.MakeTree("TCRFluxFit");
  h.MakeTree("TSPECFITF1");
  h.MakeTree("TBPLF1");
  h.CreateAuxiliaryFiles();
}
