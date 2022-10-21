#include "TROOT.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TString.h"
#include "TLegend.h"
#include <iostream>
#include <fstream>
#include "bininfo.h"


/*
	This code is for 2 particle correlation for mid-forward rapidity in pp collisions from EPOS and AMTP models
*/


void h2dLMTempFitOne(TH1D* hY ,TH1D* hY_MB, int ic, int iptt, int ig); 

Double_t Chi2(TH1D *hY_a, TF1 *fFit);

// Initializations and constants

const int numbOfFVar = 500; // Number of F values
Double_t factorF[numbOfFVar];
double F_min = 0.9;
double F_max = 2;
TString errNames[] = {"fit_G_err","fit_V2_err ","fit_V3_err "};
TString paramNames[] = {"G", "v22", "v33"};

//---------------------------------------
// Test with ALICE data
//---------------------------------------
void h2dLMTempFit(){

 	// Loading data
	TFile *fIn = new TFile ("output/fout_2pacdist_hist_proj_pp13TeV.root", "read");
	TH1D* hY = (TH1D*) fIn->Get("h2pac_corrected_dist_ETAG00C02PTT00"); // HM
	TH1D* hY_MB = (TH1D*) fIn->Get("h2pac_corrected_dist_ETAG00C00PTT00"); // LM
	h2dLMTempFitOne(hY, hY_MB, 0, 0, 0);

}
//---------------------------------------
// Two inputs : HM and LM-template, need ic and iptt for multiplicity bin and pt bins
//---------------------------------------
void h2dLMTempFitOne(TH1D* hY ,TH1D* hY_MB, int ic, int iptt, int ig) {
	// F factor values
	Double_t stepsize = (F_max - F_min)/(double) numbOfFVar;
	for (int i = 0; i <= numbOfFVar; i++) factorF[i] = F_min + (i*stepsize);
	cout << Form("%.2f<F<%.2f, NF=%d, step=%.4f\n",F_min,F_max,numbOfFVar,stepsize) << endl;

 	TH1D* hY_a;
 	Double_t chi2val;
	TF1 *fitvn_s[NH];
	TH1D* hY_a_G;
	Double_t vn[NH];
	Double_t vnError[NH];
	Double_t params[sizeof(paramNames)/sizeof(*paramNames)]; // Reserve an array of same size as paramNames
	TH1D* hFitTotal;


	TH1D *hchiq2 = new TH1D("hchiq2","chi2 of fits",numbOfFVar,0,5);

	//	FIT FUNCTION FOR PARAMETERS
 	string cosine = "[0]*(1";
	for (int i = 0; i < NH; i++) {
		ostringstream app;
		app << "+2*[" << i+1 << "]*TMath::Cos(" << i+2 << "*x)"; // Coefficients are Vn(delta)phi = Vn^2
		string append = app.str();
		cosine = cosine + append;
	}
	cosine = cosine + ")";
	cout << cosine << endl;
	const char* cos = cosine.c_str();
  

	TF1* fFit[numbOfFVar]; // Array for diff fit functions, size is numbOfFVar

 	Double_t chiq_min = -1.;
 	int index = -1;

	for (int j = 0; j < numbOfFVar; j++){ 
		fFit[j]= new TF1(Form("fFit%d",j), cos, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
		fFit[j]->SetParameter(0, 1);
		for (int i = 0; i < NH; i++) fFit[j]->SetParName(i+1, paramNames[i]); 
		for (int i = 0; i < NH; i++) fFit[j]->SetParameter(i+1, TMath::Power(1.0 - (i*0.06),2)); // Initial Vn values are Vn(delta)phi = Vn^2

 		chi2val = 0;
 		hY_a = (TH1D*) hY->Clone(); 
 		hY_a->Add(hY_MB, -factorF[j]);

 		// Fit options: Q = don't print, N = don't draw  
 		hY_a->Fit(fFit[j], "QN", "", -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0); // Adds numbOfFVar diff fits to an array
    	chi2val = fFit[j]->GetChisquare() / fFit[j]->GetNDF(); // Chi-square test ( X^2 / NDF ) 
 		hchiq2->Fill(chi2val);		
 		cout << Form("F[%d] = %.4f, chi2=%0.4f",j,factorF[j],chi2val) << endl;

   		if(chi2val < chiq_min || chiq_min < 0.){
        	chiq_min = chi2val;
        	index = j;
   		}	
	}

	cout << Form("F value=%.4f, chi2_min=%0.4f, index=%d/%d\n",factorF[index], chiq_min, index, numbOfFVar)<< endl;

	TF1* fBestFit = (TF1*) fFit[index]->Clone();

	delete *fFit;
	delete hY_a;
  
	params[0] = fBestFit->GetParameter(0);
	for (int l = 1; l < NH; l++) params[l] = fBestFit->GetParameter(l); // Saving Vn^2 values

    cout << "Saving results in to the ROOT file" << endl;

 	hFitTotal = (TH1D*) hY_MB->Clone();
 	for (int k = 1; k <= hY_MB->GetNbinsX(); k++) 
 	{
 		Double_t ylm = hY_MB->GetBinContent(k); // Taking k'th bin value
 		Double_t x = hY_MB->GetXaxis()->GetBinCenter(k);
 		Double_t tot =  fBestFit->Eval(x) + (factorF[index]*ylm); // Adding all up
 		hFitTotal->SetBinContent(k, tot);
 	}

	// F*Y_LM + G DISTRIBUTION
    hY_a_G = (TH1D*) hY_MB->Clone(); 
	for (int k = 1; k <= hY_a_G->GetNbinsX(); k++) {
			double value = hY_a_G->GetBinContent(k);
			value = value*factorF[index] + params[0];
			hY_a_G->SetBinContent(k, value);
	}
	hY_a_G->SetMarkerStyle(24);



	// SAVINGS (Signal, Fit, F*Y_LM+G)
	TFile *fOut = new TFile (Form("2.output_LMfits/out_LMtemplate_C%02dPTT%02dETAG%02d.root",ic,iptt, ig), "recreate");
	hY->Write("hDphiHM"); // SIGNAL
	fBestFit->Write("fFit_best"); 
	hY_a_G->Write("hY_a_G"); // F*Y_LM+G
	hFitTotal->Write("hFitTotal");
	hchiq2->Write();
	TString strF = Form("%.2f",factorF[index]);
	//strF->Write();

	// PRODUCING V2 AND V3 HARMONICS AND SAVING 
	Double_t Y_LM_min = hY_MB->GetMinimum(0);
	Double_t ScaleFYmin = factorF[index]*Y_LM_min;
	// Saving vn results to text file
	fstream file;
	TString outtextname = Form("2.output_LMfits/out_LMtemplate_C%02dPTT%02dETAG%02d.txt",ic,iptt, ig);
    file.open(outtextname.Data(), ios_base::out);

	for (Int_t n=0; n<NH; n++)
	{
		TString formula = Form("[0]*(1 + 2*[1]*TMath::Cos(%d*x)) + [2]",n+2);					
		fitvn_s[n]= new TF1(Form("fit_s_v%d", n+2),formula, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
		vn[n] = fFit[index]->GetParameter(n+1); // v2 v3 correctly
		vnError[n] = fFit[index]->GetParError(n+1);															
		fitvn_s[n]->SetParameter(1, vn[n]);
		fitvn_s[n]->SetParameter(0, params[0]);
		fitvn_s[n]->SetParameter(2, ScaleFYmin);
		fitvn_s[n]->Write();
	}
	
	// Writes ig, ic, iptt, vns, vn error to a txt file for z03.makegraphs.C to read
	TString texttmp = Form("%d %d %d %.5f %.5f %.5f %.5f\n", ig, ic, iptt, vn[0], vnError[0], vn[1], vnError[1]);
	file << texttmp.Data();
	file.close();

	TCanvas *c2 = new TCanvas ("hFit", "hHM", 1);
	hY->SetMinimum(hY->GetMinimum()*0.98);
	hY->Draw();
	hY_a_G->SetLineColor(2);
	hY_a_G->Draw("same");
	fitvn_s[0]->Draw("same");
	fitvn_s[1]->Draw("same");
	hFitTotal->SetLineColor(4);
	hFitTotal->Draw("same");

	//gROOT->ProcessLine("TBrowser r"); // Opens TBrowser when done
	
} // PROGRAM ENDS HERE

/* 

STEPS IN THE ALGORITHM:

1. 	Loads two input histos
2. 	Creates G(fourier) fit function
3. 	Gives fit some initial values
4. 	Multiplies hY_MB histo with F_i value
5. 	Substracts hY_MB from hY histo to create hY' histo
6. 	Fits G(fourier) to hY'
7. 	Calculates Chi2 value
8. 	Compares Chi2 value to the previous best Chi2 value
9. 	Saves histos and fits to .root file
10.	Outputs best chi2 value and associated parameters on screen

*/
