#include "2pc_MFR.C"
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TString.h>
#include "arrlist.h"

/*	
	Author: Teemu Kallio
	Version: 1.0
	Date: 03.11.2022

	Program info:

*/

TH1D* hDeltaphi[nbins_mult][nbins_pt];

// 


void loadingData() {

	TFile *fIn = new TFile ("2.output_projections/fout_corr_dist_proj_hist_AMPT_pp13TeV_grp003_pT_try000.root", "read"); // input file

	for (int ic = 0; ic < nbins_mult; ic++) {
		for (int iptt = 0; iptt < nbins_pt; iptt++) {

				hDeltaphi[ic][iptt] = (TH1D*) fIn->Get(Form("h2dCorrProjC%02dPTT%02d", ic, iptt));

			
		} // iptt
	} // ic

	fIn = new TFile ("2.output_projections/fout_hist_proj_pp13TeV_set00_grp000_pT_try000.root", "read"); // input file

	for (int ic = 0; ic < nbins_mult; ic++) {
		for (int iptt = 0; iptt < nbins_pt; iptt++) {

				hDeltaphi[ic][iptt] = (TH1D*) fIn->Get(Form("h2dCorrProjC%02dPTT%02d", ic, iptt));

			
		} // iptt
	} // ic

}

void AnaAllCentPtBins() {

	loadingData();

		for (int ic = 0; ic < nbins_mult; ic++) {
			for (int iptt = 0; iptt < nbins_pt; iptt++) {
				cout << "-------------------------------------------------" << endl;
				cout << Form("ic:%02d_iptt:%02d", ic, iptt) << endl;
				h2dLMTempFitOne(hDeltaphi[ic][iptt],hDeltaphi[0][iptt], ic, iptt); 
				cout << "-------------------------------------------------" << endl;
			}
		}
	

}