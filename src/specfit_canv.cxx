// Dmitri Ivanov <dmiivanov@gmail.com>

#include "TCanvas.h"
#include "TVirtualPad.h"
#include "TCollection.h"
#include "TString.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCollection.h"
#include "TObjArray.h"
#include "specfit_canv.h"


// Collection of routines for manipulating canvases
TObjArray specfit_canv::AllCanvases;

Int_t specfit_canv::_arrangement_factor_ = 8;
  
Int_t specfit_canv::count_subpads(TVirtualPad *pad)
{
  if(!pad)
    return 0;
  Int_t npads = 0;
  TIter next(pad->GetListOfPrimitives());
  for (TObject *obj = 0; (obj=next()); )
    {
      if(obj->InheritsFrom(TVirtualPad::Class()))
	npads ++;
    }
  return npads;
}


TCanvas* specfit_canv::get_canvas(Int_t canvas_number)
{ 
  if(canvas_number < 1 || canvas_number > (Int_t) AllCanvases.GetEntries())
    {
      fprintf(stderr,"error: zoomin: canvas_number must be in 1-%d range\n",
	      (Int_t) AllCanvases.GetEntries());
      return 0;
    }
  return (TCanvas *) AllCanvases[canvas_number - 1];
}

Int_t specfit_canv::get_ncanvases()
{
  return (Int_t) AllCanvases.GetEntries();
}

Int_t specfit_canv::get_npads(Int_t canvas_number)
{
  TCanvas *canv = get_canvas(canvas_number);
  return count_subpads(canv);
}


void specfit_canv::Update(TVirtualPad *p)
{
  p->Modified();
  p->Update();
  Int_t nsubpads = count_subpads(p);
  if(nsubpads)
    {
      for (Int_t i = 1; i <= nsubpads; i++)
	{
	  if(p->cd(i))
	    {
	      p->cd(i)->Modified();
	      p->cd(i)->Update();
	    }
	}
    }
}
  
void specfit_canv::Update(Int_t canvas_number)
{
  TCanvas *canv = get_canvas(canvas_number);
  if(!canv)
    return;
  Update(canv);
}
  
void specfit_canv::Update()
{
  for (Int_t icanvas=0; icanvas < (Int_t)AllCanvases.GetEntries(); icanvas++)
    {
      TCanvas* canv= (TCanvas *)AllCanvases[icanvas];
      Update(canv);
    }
}

TVirtualPad* specfit_canv::cd(Int_t canvas_number, Int_t subpadnumber)
{
  TCanvas *canv = get_canvas(canvas_number);
  if(!canv)
    return 0;
  return canv->cd(subpadnumber);
}

void specfit_canv::zoomin(TCanvas* canv, Int_t xsize, Int_t ysize)
{  
  canv->SetWindowSize(xsize,ysize);
  canv->SetWindowPosition(1,1);
  canv->SetWindowPosition(0,0);
  Update(canv);
}

void specfit_canv::zoomin(Int_t canvas_number, Int_t xsize, Int_t ysize)
{ 
  TCanvas *canv = get_canvas(canvas_number);
  if(!canv)
    return;
  zoomin(canv,xsize,ysize);
}

void specfit_canv::zoomin(TCanvas* canv, Int_t xysize)
{ 
  zoomin(canv,xysize,xysize); 
}

void specfit_canv::zoomin(Int_t canvas_number, Int_t xysize)
{ 
  zoomin(canvas_number,xysize,xysize); 
}

void specfit_canv::zoom_out(TCanvas *canv)
{
  zoomin(canv,300);
}

void specfit_canv::zoom_out(Int_t canvas_number)
{
  TCanvas *canv = get_canvas(canvas_number);
  if(!canv)
    return;
  zoom_out(canv);
}

void specfit_canv::zoom_out()
{
  for (Int_t icanvas = 0; icanvas < (Int_t) AllCanvases.GetEntries(); icanvas++)
    {
      TCanvas *canv = (TCanvas *)AllCanvases[icanvas];
      zoom_out(canv);
    }
}

void specfit_canv::bigTitle()
{
  gStyle->SetTitleFontSize(0.1);
  Update();
}


void specfit_canv::smallTitle()
{
  gStyle->SetTitleFontSize(0.05);
  Update();
}
  
void specfit_canv::Clear(TCanvas *canv)
{
  canv->Clear();
  Update(canv);
}

void specfit_canv::Clear(Int_t canvas_number)
{
  TCanvas *canv = get_canvas(canvas_number);
  if(!canv)
    return;
  Clear(canv);
}
  
void specfit_canv::Clear()
{
  for (Int_t icanvas=0; icanvas < (Int_t) AllCanvases.GetEntries(); icanvas++)
    {
      TCanvas* canv=(TCanvas *)AllCanvases[icanvas];
      Clear(canv);
    }
}

bool specfit_canv::save_canvas(TCanvas* canv, const char* fname, Int_t xysize)
{
  if(!canv)
    return false;
  canv->cd(); // so that it works correctly in batch mode
  if(xysize)
    zoomin(canv,xysize);
  if(strlen(fname))
    canv->SaveAs(fname);
  else
    {
      TString plot_fname = canv->GetName();
      plot_fname += ".png";
      canv->SaveAs(plot_fname);
    }
  return true;
}


void specfit_canv::set_glob_style()
{
  gROOT->SetStyle("Plain");
  gStyle->SetTitleSize(0.055,"X");
  gStyle->SetTitleSize(0.055,"Y");
  gStyle->SetTitleOffset(1.2,"X");
  gStyle->SetTitleOffset(1.2,"Y");
  gStyle->SetLineWidth(3);
  gStyle->SetPalette(1,0);
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);
  gStyle->SetTitleStyle(1111);
  gStyle->SetTitleX(0.1);
  gStyle->SetTitleW(0.8);
  gStyle->SetTitleBorderSize(0);
}


void specfit_canv::set_tick(TVirtualPad *p, Int_t x, Int_t y)
{
  p->SetTickx(x);
  p->SetTicky(y);
  Update(p);
}
  
void specfit_canv::set_tick(Int_t x, Int_t y)
{
  for (Int_t icanvas = 0; icanvas < (Int_t) AllCanvases.GetEntries(); icanvas++)
    {
      TCanvas *canv = (TCanvas *)AllCanvases[icanvas];
      set_tick(canv,x,y);
    }
}

void specfit_canv::set_grid(TVirtualPad* p, Int_t x, Int_t y)
{
  p->SetGridx(x);
  p->SetGridy(y);
  Update(p);
}
  
void specfit_canv::set_grid(Int_t x, Int_t y)
{
  for (Int_t icanvas = 0; icanvas < (Int_t) AllCanvases.GetEntries(); icanvas++)
    {
      TCanvas *canv = (TCanvas *)AllCanvases[icanvas];
      set_grid(canv,x,y);
    }
}

void specfit_canv::set_log(TVirtualPad* p, Int_t x, Int_t y, Int_t z)
{
  p->SetLogx(x);
  p->SetLogy(y);
  p->SetLogz(z);
  Update(p);
}
  
void specfit_canv::set_log(Int_t x, Int_t y, Int_t z)
{
  for (Int_t icanvas = 0; icanvas < (Int_t) AllCanvases.GetEntries(); icanvas++)
    {
      TCanvas *canv = (TCanvas *)AllCanvases[icanvas];
      set_log(canv,x,y,z);
    }
}
  
void specfit_canv::adjust_margins(TVirtualPad *p,
				  Double_t left_margin, 
				  Double_t bottom_margin, 
				  Double_t right_margin, 
				  Double_t top_margin)
{
  p->SetLeftMargin(left_margin);
  p->SetBottomMargin(bottom_margin);
  p->SetRightMargin(right_margin);
  p->SetTopMargin(top_margin);
  Update(p);
}
  
void specfit_canv::adjust_margins(Double_t left_margin, 
				  Double_t bottom_margin, 
				  Double_t right_margin, 
				  Double_t top_margin)
{
  for (Int_t icanvas = 0; icanvas < (Int_t) AllCanvases.GetEntries(); icanvas++)
    {
      TCanvas *canv = (TCanvas *)AllCanvases[icanvas];
      adjust_margins(canv, left_margin, bottom_margin, right_margin,top_margin);
    }
}

void specfit_canv::arrange_canvases(Int_t arrangement_factor)
{

  if(arrangement_factor > 8)
    {
      fprintf(stderr,"maximum arrangement factor is 8\n");
      arrangement_factor = 8;
    }

  if(arrangement_factor < 1)
    {
      fprintf(stderr,"minimum arrangement factor is 1\n");
      arrangement_factor = 1;
    }

  _arrangement_factor_ = arrangement_factor;
  

  // make windows smaller
  zoom_out();
    
  // arrange window positions
  for (Int_t icanvas = 0; icanvas < AllCanvases.GetEntries(); icanvas++)
    {
      // specific arrangement of windows
      Int_t icanvas_mod_arrangement = icanvas % arrangement_factor;
      if(icanvas_mod_arrangement == 0)
	((TCanvas *)AllCanvases[icanvas])->SetWindowPosition(0,0);
      if(icanvas_mod_arrangement == 1)
	((TCanvas *)AllCanvases[icanvas])->SetWindowPosition(309,0); 
      if(icanvas_mod_arrangement == 2)
	((TCanvas *)AllCanvases[icanvas])->SetWindowPosition(618,0);
      if(icanvas_mod_arrangement == 3)
	((TCanvas *)AllCanvases[icanvas])->SetWindowPosition(0,350);
      if(icanvas_mod_arrangement == 4)
	((TCanvas *)AllCanvases[icanvas])->SetWindowPosition(309,350);
      if(icanvas_mod_arrangement == 5)
	((TCanvas *)AllCanvases[icanvas])->SetWindowPosition(618,350);
      if(icanvas_mod_arrangement == 6)
	((TCanvas *)AllCanvases[icanvas])->SetWindowPosition(927,0);
      if(icanvas_mod_arrangement == 7)
	((TCanvas *)AllCanvases[icanvas])->SetWindowPosition(927,350);
      Update((TCanvas *)AllCanvases[icanvas]);
    }
}

void specfit_canv::arrange_canvases()
{
  arrange_canvases(_arrangement_factor_);
}

void specfit_canv::Iconify(TCanvas *canv)
{
  canv->Iconify();
  Update(canv);
}
  
void specfit_canv::Iconify(Int_t canvas_number)
{
  TCanvas *canv = get_canvas(canvas_number);
  if(!canv)
    return;
  Iconify(canv);
}

void specfit_canv::Iconify()
{ 
  for (Int_t icanvas=0; icanvas < (Int_t) AllCanvases.GetEntries(); icanvas++)
    {
      TCanvas *canv = (TCanvas *)AllCanvases[icanvas];
      Iconify(canv);
    }
}

void specfit_canv::Show(TCanvas *canv)
{
  canv->Show();
  Update(canv);
}


void specfit_canv::Show(Int_t canvas_number)
{
  TCanvas *canv = get_canvas(canvas_number);
  if(!canv)
    return;
  Show(canv);
}
  
void specfit_canv::Show()
{ 
  for (Int_t icanvas=0; icanvas < (Int_t) AllCanvases.GetEntries(); icanvas++)
    {
      TCanvas *canv = (TCanvas *)AllCanvases[icanvas];
      Show(canv);
    }
}
  
bool specfit_canv::save_canvas(Int_t canvas_number, const char* fname, Int_t xysize)
{  
  TCanvas *canv = get_canvas(canvas_number);
  if(!canv)
    return false;
  return save_canvas(canv,fname,xysize);
}

void specfit_canv::save_plots(const char* basename, 
			      const char* fExt,
			      Int_t xysize)
{
  for (Int_t icanvas = 0; icanvas < (Int_t) AllCanvases.GetEntries(); icanvas++)
    {
      TCanvas *canv = (TCanvas *)AllCanvases[icanvas];
      if(xysize)
	zoomin(canv,xysize);
      canv->cd(); // so that it works correctly in batch mode
      TString fname = basename;
      fname += canv->GetName();
      fname+=fExt;
      canv->SaveAs(fname);
    }
  arrange_canvases();
}


void specfit_canv::save_plot(TCanvas *canv,
			     const char* fname,
			     Int_t xysize)
{
  save_canvas(canv,fname,xysize);
  arrange_canvases();
}


void specfit_canv::init_canvases(Int_t ncanvases, Int_t arrangement_factor, Bool_t tick, Bool_t grid)
{

  Int_t icanvas_start   = (AllCanvases.GetEntries() ? (Int_t)AllCanvases.GetEntries() : 0);
  Int_t ncanvases_total = (AllCanvases.GetEntries() ? AllCanvases.GetEntries() + ncanvases : ncanvases);
 
  // apply gStyle settings
  set_glob_style();
    
  // allocate all canvases
  for (Int_t icanvas = icanvas_start; icanvas < ncanvases_total; icanvas ++)
    {
      TString canvName = TString::Format("c%d",icanvas+1);
      TCanvas *canv = new TCanvas(canvName,canvName,0,10,700,700);
      if(tick)
	set_tick(canv,1,1);
      if(grid)
	set_grid(canv,1,1);
      AllCanvases.Add(canv);
    }
    
  // make sure canvases have reasonable margins so that the titles of the axes
  // are displayed properly
  adjust_margins();
    
  // arrange canvases on the screen
  arrange_canvases(arrangement_factor);
    
}
  

NamespaceImp(specfit_canv);
