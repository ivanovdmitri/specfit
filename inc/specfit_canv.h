// Dmitri Ivanov <dmiivanov@gmail.com>

#ifndef _specfit_canv_h_
#define _specfit_canv_h_

#include "TCanvas.h"
#include "TVirtualPad.h"
#include "TCollection.h"
#include "TString.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCollection.h"
#include "TObjArray.h"

namespace specfit_canv
{
  // Collection of routines for manipulating canvases
  extern TObjArray AllCanvases;
  extern Int_t _arrangement_factor_;
  Int_t count_subpads(TVirtualPad *pad);
  TCanvas* get_canvas(Int_t canvas_number);
  Int_t get_ncanvases();
  Int_t get_npads(Int_t canvas_number);
  void Update(TVirtualPad *p);
  void Update(Int_t canvas_number);  
  void Update();
  TVirtualPad* cd(Int_t canvas_number, Int_t subpadnumber = 0);
  void zoomin(TCanvas* canv, Int_t xsize, Int_t ysize);
  void zoomin(Int_t canvas_number, Int_t xsize, Int_t ysize);
  void zoomin(TCanvas* canv, Int_t xysize=600);
  void zoomin(Int_t canvas_number, Int_t xysize=600);
  void zoom_out(TCanvas *canv);
  void zoom_out(Int_t canvas_number);
  void zoom_out();
  void bigTitle();
  void smallTitle(); 
  void Clear(TCanvas *canv);
  void Clear(Int_t canvas_number); 
  void Clear();
  bool save_canvas(TCanvas* canv, const char* fname = "", Int_t xysize=0);
  void set_glob_style();
  void set_tick(TVirtualPad *p, Int_t x, Int_t y);
  void set_tick(Int_t x=1, Int_t y=1);
  void set_grid(TVirtualPad* p, Int_t x, Int_t y);
  void set_grid(Int_t x=1, Int_t y=1);
  void set_log(TVirtualPad* p, Int_t x, Int_t y, Int_t z);
  void set_log(Int_t x=1, Int_t y=1, Int_t z=1); 
  void adjust_margins(TVirtualPad *p,
		      Double_t left_margin = 0.15, 
		      Double_t bottom_margin = 0.15, 
		      Double_t right_margin = 0.1, 
		      Double_t top_margin = 0.1);  
  void adjust_margins(Double_t left_margin = 0.15, 
		      Double_t bottom_margin = 0.15, 
		      Double_t right_margin = 0.1, 
		      Double_t top_margin = 0.1);
  void arrange_canvases(Int_t arrangement_factor);
  void arrange_canvases();
  void Iconify(TCanvas *canv);
  void Iconify(Int_t canvas_number);
  void Iconify();
  void Show(TCanvas *canv);
  void Show(Int_t canvas_number);
  void Show();
  bool save_canvas(Int_t canvas_number=1, const char* fname = "", Int_t xysize=0);
  void save_plots(const char* basename="plots_", 
		  const char* fExt=".png",
		  Int_t xysize = 800);
  void save_plot(TCanvas *canv,
		 const char* fname = "",
		 Int_t xysize = 800);  
  void init_canvases(Int_t ncanvases = 1, Int_t arrangement_factor = 8, Bool_t tick=true, Bool_t grid=false); 
}

#endif
