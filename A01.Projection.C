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
	Long-range correlation projections for AMPT and PYTHIA model data
*/

TH2D *h2D_corr_dist[nbins_mult][nbins_pt];
TH1D *hDeltaphi_corr[nbins_mult][nbins_pt];




void load2DHistos(TString inputname); // loads 2D histo from given .root file
void makeLongRangeCorr(TString inputname);
void writeToRoot(TString outputname);

// Program starts here
void projection() {

	TString inputname = "output/fout_corr_dist_hist_AMPT_pp13TeV_grp003_pT_try000.root"; // AMPT data
	load2DHistos(inputname);
	makeLongRangeCorr(inputname);
	TString outname = "2.output_projections/fout_corr_dist_proj_hist_AMPT_pp13TeV_grp003_pT_try000.root"; // AMPT
	writeToRoot(outname);


	inputname = "output/fout_corr_dist_hist_pp13TeV_set00_grp000_pT_try000.root"; // PYTHIA data
	load2DHistos(inputname);
	makeLongRangeCorr(inputname);
	outname = "2.output_projections/fout_hist_proj_pp13TeV_set00_grp000_pT_try000.root";
	writeToRoot(outname);
}

/* 
	Loops over .root files for same, mixed and trigger distributions and loads files to 2d histos
*/
void load2DHistos(TString inputname){

	TFile *fIn = new TFile (inputname, "read"); // input file

	if (inputname == "output/fout_corr_dist_hist_pp13TeV_set00_grp000_pT_try000.root") {

		for (int ic = 0; ic < nbins_mult; ic++){
			for (int iptt = 0; iptt < nbins_pt; iptt++){

				h2D_corr_dist[ic][iptt] = (TH2D*) fIn->Get(Form("h2d_pythia_corr_dist_dphi_deta_fmda_h_mult%02d_pt%02d", ic, iptt));

			} // iptt
		} // ic

	}

	if (inputname == "output/fout_corr_dist_hist_AMPT_pp13TeV_grp003_pT_try000.root") {

		for (int ic = 0; ic < nbins_mult; ic++){
			for (int iptt = 0; iptt < nbins_pt; iptt++){

				h2D_corr_dist[ic][iptt] = (TH2D*) fIn->Get(Form("h2d_ampt_corr_dist_dphi_deta_fmda_h_mult%02d_pt%02d", ic, iptt));

			} // iptt
		} // ic
	}
}

void makeLongRangeCorr(TString inputname) {

	if (inputname == "output/fout_corr_dist_hist_pp13TeV_set00_grp000_pT_try000.root") {

		for (int ic = 0; ic < nbins_mult; ic++){
			for (int iptt = 0; iptt < nbins_pt; iptt++){

				cout << "ic " << ic << " iptt " << iptt << endl;
				hDeltaphi_corr[ic][iptt] = (TH1D*) h2D_corr_dist[ic][iptt]->ProjectionX(Form("h2d_pythia_corr_proj_hist_C%02dPTT%02d", ic, iptt), 
					h2D_corr_dist[ic][iptt]->GetYaxis()->FindBin(1.9), h2D_corr_dist[ic][iptt]->GetYaxis()->FindBin(4.8));
				

			} // iptt
		} // ic
	}

	if (inputname == "output/fout_corr_dist_hist_AMPT_pp13TeV_grp003_pT_try000.root") {

		for (int ic = 0; ic < nbins_mult; ic++){
			for (int iptt = 0; iptt < nbins_pt; iptt++){

				cout << "ic " << ic << " iptt " << iptt << endl;
				hDeltaphi_corr[ic][iptt] = (TH1D*) h2D_corr_dist[ic][iptt]->ProjectionX(Form("h2d_ampt_corr_proj_hist_C%02dPTT%02d", ic, iptt), 
					h2D_corr_dist[ic][iptt]->GetYaxis()->FindBin(1.9), h2D_corr_dist[ic][iptt]->GetYaxis()->FindBin(4.8));
				

			} // iptt
		} // ic
	}
} 


void writeToRoot(TString outputname) {

	TFile* fOut = new TFile (outputname, "recreate");

	if (outputname == "2.output_projections/fout_hist_proj_pp13TeV_set00_grp000_pT_try000.root") {

		for (int ic = 0; ic < nbins_mult; ic++){
			for (int iptt = 0; iptt < nbins_pt; iptt++){

				hDeltaphi_corr[ic][iptt]->Write(Form("h2dPythiaCorrProjC%02dPTT%02d", ic, iptt));

			} // iptt
		} // ic
	}

	if (outputname == "2.output_projections/fout_corr_dist_proj_hist_AMPT_pp13TeV_grp003_pT_try000.root") {

		for (int ic = 0; ic < nbins_mult; ic++){
			for (int iptt = 0; iptt < nbins_pt; iptt++){

				hDeltaphi_corr[ic][iptt]->Write(Form("h2dAMPTCorrProjC%02dPTT%02d", ic, iptt));

			} // iptt
		} // ic


	}


}

