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

TH2D *h2D_same[nbins_mult][nbins_pt];
TH1D *hDeltaphi_same[nbins_mult][nbins_pt];
TH2D *h2D_mixed[nbins_mult][nbins_pt];
TH1D *hDeltaphi_mixed[nbins_mult][nbins_pt];
TH1D *hDeltaphi_trig[nbins_mult];



void load2DHistos(TString inputname); // loads 2D histo from given .root file
void makeLongRangeCorr();
void writeToRoot(TString outname);

// Program starts here
void projection() {

	TString inputname = "input/outfile_hist_pp13TeV_set00_grp000_pT_try000.root";
	load2DHistos(inputname);
	makeLongRangeCorr();
	TString outname = "output/fout_hist_proj_pp13TeV_set00_grp000_pT_try000.root";
	writeToRoot(outname);

}

/* 
	Loops over .root files for same, mixed and trigger distributions and loads files to 2d histos
*/
void load2DHistos(TString inputname = "input/outfile_hist_pp13TeV_set00_grp000_pT_try000.root"){

	TFile *fIn = new TFile (inputname, "read"); // input file

		for (int ic = 0; ic < nbins_mult; ic++){
			for (int iptt = 0; iptt < nbins_pt; iptt++){
				h2D_same[ic][iptt] = (TH2D*) fIn->Get(Form("h2d_same_dphi_deta_fmda_h_mult%02d_pt%02d", ic, iptt));
				h2D_mixed[ic][iptt] = (TH2D*) fIn->Get(Form("h2d_mixed_dphi_deta_fmda_h_mult%02d_pt%02d", ic, iptt));
				hDeltaphi_trig[ic] = (TH1D*) fIn->Get(Form("hntrig_same_fmda_h_mult%02d", ic));	

				cout << Form("pt%02d", iptt) << endl;

			} // iptt
		} // ic
	}


void makeLongRangeCorr() {

	for (int ic = 0; ic < nbins_mult; ic++){
		for (int iptt = 0; iptt < nbins_pt-1; iptt++){

			hDeltaphi_same[ic][iptt] = (TH1D*) h2D_same[ic][iptt]->ProjectionX(Form("h2dsameETAG%02dC%02dPTT%02d", 00, ic, iptt), h2D_same[ic][iptt]->GetYaxis()->FindBin(1.9), h2D_same[ic][iptt]->GetYaxis()->FindBin(4.8));
			//TH1D *hTemp_same = (TH1D*) h2D_same->ProjectionX(Form("h2dsameETAG%02dC%02diptt%02d", 00, 21, 00), h2D_same->GetYaxis()->FindBin(-3.1), h2D_same->GetYaxis()->FindBin(-1.9));
			//hDeltaphi_same->Add(hTemp_same, 1);
			
			hDeltaphi_mixed[ic][iptt] = (TH1D*) h2D_mixed[ic][iptt]->ProjectionX(Form("h2dmixETAG%02dC%02dPTT%02d", 00, ic, iptt), h2D_mixed[ic][iptt]->GetYaxis()->FindBin(1.9), h2D_mixed[ic][iptt]->GetYaxis()->FindBin(4.8));
			//TH1D *hTemp_mixed = (TH1D*) h2D_mixed->ProjectionX(Form("h2dmixETAG%02dC%02diptt%02d", 00, 21, 00), h2D_mixed->GetYaxis()->FindBin(-3.1), h2D_mixed->GetYaxis()->FindBin(-1.9));
			//hDeltaphi_mixed->Add(hTemp_mixed, 1);

		} // iptt
	} // ic
} 


void writeToRoot(TString outname = "output/fout_hist_proj_pp13TeV_set00_grp000_pT_try000.root") {

	TFile* fOut = new TFile (outname, "recreate");
	//TString htitle;
	//htitle = Form("%.0f< N_{ch}<%.0f, %.3f< p_{T} <%.3f, %.1f< |#Delta#eta| <%.1f",Mult_HMT_min[ic],Mult_HMT_max[ic],pTMin[iiptt],pTMax[iiptt], etamin[ig], etamax);
	//hDeltaphi->SetTitle(htitle);

	for (int ic = 0; ic < nbins_mult; ic++){
		for (int iptt = 0; iptt < nbins_pt-1; iptt++){

			hDeltaphi_same[ic][iptt]->Write(Form("h2dsameETAG%02dC%02dPTT%02d", 00, ic, iptt));
			hDeltaphi_mixed[ic][iptt]->Write(Form("h2dmixETAG%02dC%02dPTT%02d", 00, ic, iptt));
			hDeltaphi_trig[ic]->Write(Form("h2dtrigETAG%02dC%02d", 00, ic));

		} // iptt
	} // ic
}

