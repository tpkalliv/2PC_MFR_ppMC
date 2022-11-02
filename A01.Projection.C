#include "arrlist.h"
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TString.h>

/*	
	Author: Teemu Kallio
	Version: 1.0
	Date: 19.10.2022

	Program info:
	Long-range correlation projections
*/

TH2D *h2D_corr_dist[nbins_mult][nbins_pt];




void load2DHistos(TString inputname); // loads 2D histo from given .root file
void makeLongRangeCorr();
void writeToRoot(TString outname);

// Program starts here
void projection() {

	//TString inputname = "input/outfile_hist_pp13TeV_set00_grp000_pT_try000.root"; // PYTHIA data
	TString inputname = "output/fout_corr_dist_hist_AMPT_pp13TeV_grp003_pT_try000.root"; // AMPT data
	load2DHistos(inputname);
	makeLongRangeCorr();
	gStyle->SetPalette(55);
	//TString outname = "output/fout_hist_proj_pp13TeV_set00_grp000_pT_try000.root"; // PYTHIA
	TString outname = "output/fout_corr_dist_proj_hist_AMPT_pp13TeV_grp003_pT_try000.root"; // AMPT
	writeToRoot(outname);

}

/* 
	Loops over .root files for same, mixed and trigger distributions and loads files to 2d histos
*/
void load2DHistos(TString inputname){

	TFile *fIn = new TFile (inputname, "read"); // input file

		for (int ic = 0; ic < nbins_mult; ic++){
			for (int iptt = 0; iptt < nbins_pt; iptt++){
				h2D_corr_dist[ic][iptt] = (TH2D*) fIn->Get(Form("h2D_corr_dist_proj_dphi_deta_fmda_h_mult%02d_pt%02d", ic, iptt));
			} // iptt
		} // ic
	}


void makeLongRangeCorr() {

	for (int ic = 0; ic < nbins_mult; ic++){
		for (int iptt = 0; iptt < nbins_pt-1; iptt++){

			hDeltaphi_corr[ic][iptt] = (TH1D*) h2D_corr_dist[ic][iptt]->ProjectionX(Form("h2dsameETAG%02dC%02dPTT%02d", 00, ic, iptt), 
				h2D_corr_dist[ic][iptt]->GetYaxis()->FindBin(1.9), h2D_corr_dist[ic][iptt]->GetYaxis()->FindBin(4.8));

		} // iptt
	} // ic
} 


void writeToRoot(TString outname = "output/fout_corr_dist_hist_proj_pp13TeV_set00_grp000_pT_try000.root") {

	TFile* fOut = new TFile (outname, "recreate");
	//TString htitle;
	//htitle = Form("%.0f< N_{ch}<%.0f, %.3f< p_{T} <%.3f, %.1f< |#Delta#eta| <%.1f",Mult_HMT_min[ic],Mult_HMT_max[ic],pTMin[iiptt],pTMax[iiptt], etamin[ig], etamax);
	//hDeltaphi->SetTitle(htitle);

	for (int ic = 0; ic < nbins_mult; ic++){
		for (int iptt = 0; iptt < nbins_pt-1; iptt++){

			hDeltaphi_corr[ic][iptt]->Write(Form("h2dCorrC%02dPTT%02d", ic, iptt));

		} // iptt
	} // ic
}

