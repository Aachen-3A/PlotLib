#!/bin/env python

from lib.DukePlotALot import *
from lib.plotlib import *
from lib.configobj import ConfigObj

def main():

    basedir="/user/padeken/out/output2014_12_16_14_49/merged"
    lumi=19712
    bghists=HistSorage(basedir)
    bghists.addAllFiles(veto=["Wprime","Run2012","QCD","data"])
    #bghists.addFile("WToTauNu_ptmin500_TuneZ2Star_8TeV-pythia6-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1SIM")
    print bghists.files.keys()

    dat_hist=HistSorage(basedir)
    dat_hist.addFile("allDataMET")

    hist="byLooseCombinedIsolationDeltaBetaCorr3Hits/h1_5_byLooseCombinedIsolationDeltaBetaCorr3Hits_MT"
    bghists.getHist(hist)
    dat_hist.getHist(hist)

    xs= ConfigObj("/home/home1/institut_3a/padeken/Analysis/SirPlotAlot/xsv100.cfg")

    dat_hist.rebin(width=10)
    bghists.rebin(width=10)
    bghists.scale_cfg(lumi,xs)

    bghists.join("DY","DY")
    bghists.join("WJetsToLNu","W_1")
    bghists.join("WTo","W_2")
    bghists.join("Summer12","Rest")
    bghists.join("W_","W")

    colors=getColorList(len(bghists.hists))
    bghists.setStyle(colors=colors)
    dat_hist.getHistList()[0].SetTitle("data")



    #bag_hist, sig_hist, dat_hist, sys_hist = create_test_histos()

    #bag_hist.fillstyle = 'solid'
    #bag_hist.fillcolor = 'green'
    #bag_hist.linecolor = 'green'
    #bag_hist.linewidth = 0
    #bag_hist.yaxis.SetTitle('Events')
    #bag_hist.xaxis.SetTitle('Mass (GeV)')

    #sig_hist.fillstyle = '0'
    #sig_hist.fillcolor = 'red'
    #sig_hist.linecolor = 'red'
    #sig_hist.linewidth = 1
    #sig_hist.yaxis.SetTitle('Events')
    #sig_hist.xaxis.SetTitle('Mass (GeV)')

    #bag_hist.Scale(0.5)
    #bag_hist_2 = bag_hist.Clone(title='Background 2')
    #bag_hist_2.fillcolor = 'y'
    #bag_hist_2.linecolor = 'y'

    #sys_hist_2 = sys_hist.Clone(title='sys2')
    #sys_hist_2.Scale(0.5)

    test = plotter(hist=bghists.getHistList(),style='CMS')
    test.Add_data(dat_hist.getHistList()[0])
    #test.Add_plot('Signi',pos=0, height=15)
    #test.Add_plot('Ratio',pos=1, height=15)
    #test.Add_plot('DiffRatio',pos=2, height=15)
    #test.Add_error_hist([sys_hist_2,sys_hist], band_center = 'ref')
    test.Set_axis(xmin=200,xmax=1500)
    test.make_plot('bla_plt_Wprime.pdf')
    return 42



main()
