import numpy as np
import ROOT

import scipy
from scipy import interpolate

import sys
sys.path.append("JPyPlotRatio");


import JPyPlotRatio

'''
For figs
'''

file = sys.argv[1];
# data set is manual
datatitle = "AMPT" 
outfigname = file.replace(".root","_fit.png")

f = ROOT.TFile(file,"read");
dataTypePlotParams = [
        {'plotType':'data','color':'k','fmt':'o','markersize':5.0},
	{'plotType':'theory','facecolor':'C0','edgecolor':'C0','alpha':0.5,'linestyle':'dashdot','linecolor':'C0'},
	{'plotType':'data','color':'r','fmt':'s','markersize':3.0},
        {'plotType':'theory','facecolor':'C1','edgecolor':'C1','alpha':0.5,'linestyle':'solid','linecolor':'C1'},
        {'plotType':'theory','facecolor':'C2','edgecolor':'C2','alpha':0.5,'linestyle':'dashed','linecolor':'C2'},
        {'plotType':'theory','facecolor':'C3','edgecolor':'C3','alpha':0.5,'linestyle':'dashdot','linecolor':'C3'}
];


# define panel/xaxis limits/titles
nrow = 1;
ncol = 1;
xlimits = {0:(-1.4,4.5)};
ylimits = {0:(1.7,2.0)};
rlimits = [(0.99,1.01)];


histnames = ["hDphiHM","hFitTotal","hY_a_G", "fit_s_v2","fit_s_v3" ];
#ls 2.output_LMfits/*.root | awk '{com=sprinf("python3 LMTempFitFigs2.py %s",$1); print(com)}'
# add here the histogram names for each pad
# add labels for each pad
plables = [ ""
	#		"$1.0 < p_\\mathrm{T,trigg(assoc)} < 2.0$","$2.0 < p_\\mathrm{T,trigg(assoc)} < 3.0$","$3.0 < p_\\mathrm{T,trigg(assoc)} < 4.0$",
	#		"$1.0 < p_\\mathrm{T,trigg(assoc)} < 2.0$","$2.0 < p_\\mathrm{T,trigg(assoc)} < 3.0$","$3.0 < p_\\mathrm{T,trigg(assoc)} < 4.0$"
		 ];
# model names : for histonames in ROOT file
#modelStrInROOT = ["pythiadefault","stringshoving","eposlhc"]; # for data reading
modelStr = ["AMPT near-side","AMPT away-side"]; # for legend

#TypeName =["Signal", "Fit", "$F$LM + $G$","$G$(1 + $v_{2,2}$)","$G$(1 + $v_{3,3}$)" ];
TypeName =["Signal (0--0.1\%)", "Fit","$FY_{\\mathrm{LM}} + G$", "$G(1+2v_{2,2}cos(2\\Delta\\varphi))$ \n $+ FY_{\\mathrm{LM,min}}$","$G(1+2v_{3,3}cos(3\\Delta\\varphi))$ \n $+ FY_{\\mathrm{LM,min}}$"];

#xtitle = ["$p_\\mathrm{T,trig(assoc)} (\\mathrm{GeV}/c)$"];
xtitle = ["$\\Delta\\varphi (\\mathrm{rad})$"];
ytitle = ["$\\frac{1}{N_{\\mathrm{trig}}}\\frac{\\mathrm{d}N^{\\mathrm{pair}}}{\\mathrm{d}\\Delta\\varphi}$"];

# Following two must be added
toptitle = "pp $\\sqrt{s}$ = 13 TeV"; # need to add on the top


# getting result informations
# htitle =  15< N_{ch}<35, 4.000< p_{T} <6.000, 1.6< |#Delta#eta| <1.8
hist = f.Get("hDphiHM")
htitle = hist.GetTitle()
print("htitle = ", htitle)
hinfo = htitle.split(",")
TypeName[0] = "$"+hinfo[0]+"$"
hinfo[1] = hinfo[1].replace(".000",".0")
dataDetail = "$"+hinfo[1]+" (\\mathrm{GeV}/c)$"
hinfo[2] = hinfo[2].replace("#","\\")
dataDetailEta = "$"+hinfo[2]+"$"


plot = JPyPlotRatio.JPyPlotRatio(panels=(nrow,ncol),
	panelsize=(12,12), # change the size
	#rowBounds=ylimits,  # for nrow
	colBounds=xlimits,  # for ncol
	panelLabel=plables,  # nrowxncol
	ratioBounds=rlimits,# for nrow
#	disableRatio=[0],
	ratioPlot=True,
	panelLabelLoc=(0.1,0.87),panelLabelSize=10,panelLabelAlign="left",
	legendPanel=0,
	legendLoc=(0.16,0.75),legendSize=11,xlabel=xtitle[0],ylabel=ytitle[0]);



#plot.EnableLatex(True);

plotMatrix = np.empty((nrow,ncol),dtype=int);

#plot.GetAxes(0).set_xticks([0,1,2,3,4,5]);


for d in range(0,5):
	plot.GetAxes(0).xaxis.set_ticks_position('both');
	plot.GetAxes(0).yaxis.set_ticks_position('both');
	gr = f.Get("{}".format(histnames[d]));
	data = plot.Add(0,gr,**dataTypePlotParams[d],label=TypeName[d]);
	if(d==0):
		grFit = f.Get("{}".format(histnames[1]));
		fit = plot.Add(0,grFit,**dataTypePlotParams[1],label=TypeName[1]);
		plot.Ratio(data, fit );


f.Close();

plot.GetPlot().text(0.38,0.83,datatitle,fontsize=11);
plot.GetPlot().text(0.38,0.79,toptitle,fontsize=11);
plot.GetPlot().text(0.38,0.75,dataDetail,fontsize=11);
plot.GetPlot().text(0.38,0.71,dataDetailEta,fontsize=11);
plot.GetRatioAxes(0).xaxis.set_ticks_position('both');
plot.GetRatioAxes(0).yaxis.set_ticks_position('both');
#plot.GetAxes(0).xticks(rotation=45)
plot.Plot();


plot.Save(outfigname);
#plot.Show();