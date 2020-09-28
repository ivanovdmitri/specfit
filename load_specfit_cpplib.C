// Dmitri Ivanov <dmiivanov@gmail.com>

// Loads the shared library with the cosmic ray spectrum fit classes into ROOT. 

{
  TString SPECFIT= (gSystem->Getenv("SPECFIT") ? 
		    TString(gSystem->Getenv("SPECFIT")) : 
		    TString(gSystem->DirName(__FILE__)));
  if(gSystem->Load(SPECFIT+"/lib/libspecfit") < 0)
    {
      fprintf(stderr,"\n\n");
      fprintf(stderr,"POSSIBLE CAUSES: specfit environment not");
      fprintf(stderr," set or specfit not compiled. \n");
      fprintf(stderr,"SUGGESTIONS:\n");
      fprintf(stderr,"source /full/path/to/specfit/specfit_env.sh\n");
      fprintf(stderr,"(add this line to %s/.bashrc",gSystem->Getenv("HOME"));
      fprintf(stderr," for the future)\n");
      fprintf(stderr,"And/or re-install specfit:\n");
      fprintf(stderr,"cd /full/path/to/specfit\n");
      fprintf(stderr,"make clean; make -j3\n");
      fprintf(stderr,"\n");
      exit(2);
    }
  
  // ROOT 6 needs to have the header file included
  if(gROOT->GetVersionInt() / 10000 >= 6)
    gInterpreter->ProcessLine(TString::Format("#include \"%s/inc/specfit.h\"",SPECFIT.Data()));
  
  // making sure that vector include is available when dealing with older ROOT versions.
  if(gROOT->GetVersionInt() / 10000 < 6)
    gInterpreter->ProcessLine("#include <vector>");
}
