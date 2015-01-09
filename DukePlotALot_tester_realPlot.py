#!/bin/env python

from lib.DukePlotALot import *
from lib.plotlib import HistSorage,getColorList,getDictValue
import matplotlib.pyplot as plt
from lib.configobj import ConfigObj
try:
    from collections import OrderedDict
except ImportError:
    from lib.ordered import OrderedDict

def main():

    basedir="/user/padeken/out/output2014_12_16_14_49/merged"
    lumi=19712

    xs= ConfigObj("/home/home1/institut_3a/padeken/Analysis/SirPlotAlot/xsv100.cfg")
    bghists=HistSorage(xs,lumi,basedir)
    bghists.setDataDriven("dataDrivenQCD")

    bglist=OrderedDict()
    bglist["Diboson"]=['WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1SIM',
     'ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1SIM',
     'WZtoAnything_ptmin500_TuneZ2Star_8TeV-pythia6-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1SIM',
     'ZZtoAnything_ptmin500_TuneZ2Star_8TeV-pythia6-tauola_Summer12_DR53X-PU_S10_START53_V7A-v2SIM',
     'WWtoAnything_ptmin500_TuneZ2Star_8TeV-pythia6-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1SIM',
     'WZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1SIM']
    bglist["DY"]=['DYToTauTau_M_10To20_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1SIM',
     'DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1SIM',
     'DYToTauTau_M-100to200_TuneZ2Star_8TeV-pythia6-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1SIM',
     'DYToTauTau_M-200to400_TuneZ2Star_8TeV-pythia6-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1SIM',
     'DYToTauTau_M-400to800_TuneZ2Star_8TeV-pythia6-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1SIM',
     'DYJetsToLL_PtZ-50To70_TuneZ2star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1SIM',
     'DYJetsToLL_PtZ-70To100_TuneZ2star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2SIM',
     'DYJetsToLL_PtZ-100_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2SIM',
     'DYJetsToLL_PtZ-180_TuneZ2star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7C-v1SIM']
    bglist["Top"]=['Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1SIM',
     'T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1SIM',
     'T_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1SIM',
     'TT_CT10_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v2SIM',
     'Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1SIM',
     'Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1SIM',
     'T_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1SIM']
    bglist["QCD jet"]=["dataDrivenQCD"]
    bglist["W"]=['WToTauNu_ptmin500_TuneZ2Star_8TeV-pythia6-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1SIM',
     'WToTauNu_ptmin100_ptmax500_TuneZ2Star_8TeV-pythia6-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1SIM',
     'WJetsToLNu_PtW-70To100_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1SIM',
     'WToENu_ptmin500_TuneZ2Star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1SIM',
     'WJetsToLNu_PtW-180_TuneZ2star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7C-v1SIM',
     'WToMuNu_ptmin500_TuneZ2Star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1SIM',
     'WToMuNu_ptmin100_ptmax500_TuneZ2Star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1SIM',
     'WToENu_ptmin100_ptmax500_TuneZ2Star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1SIM',
     'WJetsToLNu_PtW-100_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1SIM',
     'WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2SIM',
     'WJetsToLNu_PtW-50To70_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1SIM']

    colorList={}
    colorList["W"]="lightblue"
    colorList["QCD jet"]="darkblue"
    colorList["Top"]="pink"
    colorList["Diboson"]="green"
    colorList["DY"]="red"

    #print bglist
    bghists.addFileList(bglist)



    dat_hist=HistSorage(xs,lumi,basedir,isData=True)
    dat_hist.addFile("allDataMET")

    hists=["byLooseCombinedIsolationDeltaBetaCorr3Hits/h1_5_byLooseCombinedIsolationDeltaBetaCorr3Hits_MT",
    "byLooseCombinedIsolationDeltaBetaCorr3Hits/h1_5_byLooseCombinedIsolationDeltaBetaCorr3Hits_tau_pt",
    "byLooseCombinedIsolationDeltaBetaCorr3Hits/h1_5_byLooseCombinedIsolationDeltaBetaCorr3Hits_met_et",
    ]

    binning={
            "_pt":10,
            "_MT":20,
            "_met_et":30,
    }

    xranges={
            "_pt":[0,900],
            "_MT":[200,1500],
            #"_met_et":[0,900],
    }
    for hist in hists:
        print hist
        bghists.clearHists()
        dat_hist.clearHists()
        bghists.getHist(hist)
        dat_hist.getHist(hist)
        binf=getDictValue(hist,binning)
        if binf is not None:
            dat_hist.rebin(width=binf)
            bghists.rebin(width=binf)
        bghists.setStyle(colors=colorList)
        dat_hist.getHistList()[0].SetTitle("data")

        test = plotter(hist=bghists.getHistList(),style='CMS')
        test.Add_data(dat_hist.getHistList()[0])
        test.Add_plot('Signi',pos=1, height=15)
        #test.Add_plot('DiffRatio',pos=1, height=15)
        #test.Add_plot('Diff',pos=2, height=15)
        #test.Add_plot('Ratio',pos=2, height=15)
        #test.Add_error_hist([sys_hist_2,sys_hist], band_center = 'ref')
        mxrange=getDictValue(hist,xranges)
        if mxrange is not None:
            test.Set_axis(xmin=mxrange[0],xmax=mxrange[1])
        name=hist.replace("/","")
        test.make_plot('%s.pdf'%(name))
    return 42



main()
