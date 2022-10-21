#include "arrlist.h"
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TString.h>

/*	
	Author: Teemu Kallio
	Version: 1.0
	Date: 21.10.2022

	Program info:
	To produce unbiased two-particle correlation distribution using same, mixed and trigger distributions.
*/

TH1D* hDphiproj_same[nbins_mult][nbins_pt];
TH1D* hDphiproj_mixed[nbins_mult][nbins_pt];
TH1D* hDphi_trig[nbins_mult];
TH1D* dist[nbins_mult][nbins_pt];

void twoParticleDist();
void toWriteRoot(TString outname);
void loadProjections(TString inputname);

void twoPCDist(){

	//TString inputname = "output/fout_hist_proj_pp13TeV_set00_grp000_pT_try000.root"; // PYTHIA
	TString inputname = "output/fout_hist_proj_pp13TeV_set00_grp000_pT_try000.root"; // AMPT
	loadProjections(inputname);
	twoParticleDist();
	//TString outname = "output/fout_2pacdist_hist_pp13TeV.root"; // PYTHIA
	TString outname = "output/fout_2pacdist_hist_AMPT_pp13TeV.root"; // AMPT
	toWriteRoot(outname);

}


void loadProjections(TString inputname){

	TFile *fIn = new TFile (inputname, "read"); // input file

		for (int ic = 0; ic < nbins_mult-1; ic++){
			for (int iptt = 0; iptt < nbins_pt-1; iptt++){

				hDphiproj_same[ic][iptt] = (TH1D*) fIn->Get(Form("h2dsameETAG%02dC%02dPTT%02d", 00, ic, iptt));
				hDphiproj_mixed[ic][iptt] = (TH1D*) fIn->Get(Form("h2dmixETAG%02dC%02dPTT%02d", 00, ic, iptt));
				hDphi_trig[ic] = (TH1D*) fIn->Get(Form("h2dtrigETAG%02dC%02d", 00, ic));

			} // iptt
		} // ic
	}

void twoParticleDist(){ 

	for (int ic = 0; ic < nbins_mult-1; ic++){
		for (int iptt = 0; iptt < nbins_pt-1; iptt++){

			Double_t trig_part_val = 0;

			for (int k = 0; k <= hDphi_trig[ic]->GetNbinsX(); k++) trig_part_val += hDphi_trig[ic]->GetBinContent(k); 

			dist[ic][iptt] = (TH1D*) hDphiproj_same[ic][iptt]->Clone();
			
			for (int bin = 0; bin <= dist[ic][iptt]->GetNbinsX(); bin++){ // Same amount of bins in same and mixed events

				Double_t same_val = dist[ic][iptt]->GetBinContent(bin);
				Double_t mixed_val = hDphiproj_mixed[ic][iptt]->GetBinContent(bin);
				Double_t divided_val = same_val / (mixed_val / trig_part_val);
				dist[ic][iptt]->SetBinContent(bin, divided_val);
			}

		} // iptt
	} // ic

}

void toWriteRoot(TString outname){

	TFile* fOut = new TFile(outname, "recreate");

	for (int ic = 0; ic < nbins_mult-1; ic++){
		for (int iptt = 0; iptt < nbins_pt-1; iptt++){

			dist[ic][iptt]->Write(Form("h2pac_corrected_dist_ETAG%02dC%02dPTT%02d", 00, ic, iptt));

		}
	} 

}