#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TString.h>
#include "arrlist.h"
#include <TStyle.h>

/*	
	Author: Teemu Kallio
	Version: 1.0
	Date: 19.10.2022

	Program info:
	Correlation distributions
*/


TH2D* h2D_same[nbins_mult][nbins_pt]; // for input data scatter plots
TH2D* h2D_mixed[nbins_mult][nbins_pt]; // for input data scatter plots
TH1D* hDeltaphi_trig[nbins_mult]; // for input data scatter plots
TH2D* h2d_corr_final[nbins_mult][nbins_pt];



void loadhistos(TString inputname); // loads 2D histo from given .root file
void makeLongRangeCorrDist(TH2D* h2d_same, TH2D* h2d_mixed, TH1D* hdeltaphi_trig, int ic, int iptt);
void writeToRoot(TString outname);


// Program starts here
void longrangecorrdist() {

	TString inputname = "input/outfile_hist_pp13TeV_set00_grp000_pT_try000.root"; // PYTHIA data
	loadhistos(inputname);
	TString outname = "output/fout_corr_dist_hist_pp13TeV_set00_grp000_pT_try000.root"; // PYTHIA
	writeToRoot(outname);

	inputname = "input/outfile_hist_AMPT_pp_13TeV_grp003_pT_try000.root"; // AMPT data
	loadhistos(inputname);
	outname = "output/fout_corr_dist_hist_AMPT_pp13TeV_grp003_pT_try000.root"; // AMPT
	writeToRoot(outname);
	
	// Sample draw
	gStyle->SetPalette(55);
	h2d_corr_final[2][9]->Draw("[cut1,-cut2],SURF2Z");

}

// Load histograms and call correlation distribution member function
void loadhistos(TString inputname){

	TFile *fIn = new TFile (inputname, "read"); // input file

		for (int ic = 0; ic < nbins_mult; ic++){
			for (int iptt = 0; iptt < nbins_pt; iptt++){
				
				h2D_same[ic][iptt] = (TH2D*) fIn->Get(Form("h2d_same_dphi_deta_fmda_h_mult%02d_pt%02d", ic, iptt));
				h2D_mixed[ic][iptt] = (TH2D*) fIn->Get(Form("h2d_mixed_dphi_deta_fmda_h_mult%02d_pt%02d", ic, iptt));
				hDeltaphi_trig[ic] = (TH1D*) fIn->Get(Form("hntrig_same_fmda_h_mult%02d", ic));	

				makeLongRangeCorrDist(h2D_same[ic][iptt],h2D_mixed[ic][iptt],hDeltaphi_trig[ic],ic,iptt); 

			} // iptt
		} // ic

}

// Correlation distribution
void makeLongRangeCorrDist(TH2D* h2d_same, TH2D* h2d_mixed, TH1D* hdeltaphi_trig, int ic, int iptt) {

	Double_t trig_content = hdeltaphi_trig->GetBinContent(1);
	Int_t scale = trig_content; 

	h2d_same->Scale(1/trig_content);
	Double_t hB_Norm = h2d_mixed->GetMaximum();
	h2d_mixed->Scale(1/hB_Norm);

	TH2D* h2d_Corr = (TH2D*) h2d_same->Clone();

	Int_t binamount = h2d_Corr->GetXaxis()->GetNbins();

	/*
		Dividing (1/N_trig)*S with B/B(0,0) to get C function 
	*/
	for (int ix = 0; ix <= h2d_Corr->GetXaxis()->GetNbins(); ix++){
		for (int iy = 0; iy <= h2d_Corr->GetYaxis()->GetNbins(); iy++){

			Double_t mixed_binval = h2d_mixed->GetBinContent(ix,iy);
			Double_t same_binval = h2d_same->GetBinContent(ix,iy);
			Double_t corr_binval = same_binval / mixed_binval;
			h2d_Corr->SetBinContent(ix, iy, corr_binval);

		}
	}

	h2d_corr_final[ic][iptt] = (TH2D*) h2d_Corr->Clone();
}

// Saving
void writeToRoot(TString outname){

	TFile *fOut = new TFile (outname, "recreate"); // input file

		for (int ic = 0; ic < nbins_mult; ic++){
			for (int iptt = 0; iptt < nbins_pt; iptt++){
				
				h2d_corr_final[ic][iptt]->SetTitle(Form("h2d_corr_dist_dphi_deta_fmda_h_mult%02d_pt%02d", ic, iptt));
				h2d_corr_final[ic][iptt]->Write();
			} // iptt
		} // ic
}