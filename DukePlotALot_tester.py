#!/bin/env python

from lib.DukePlotALot import *

def main():
    bag_hist, sig_hist, dat_hist = create_test_histos()

    bag_hist.fillstyle = 'solid'
    bag_hist.fillcolor = 'green'
    bag_hist.linecolor = 'green'
    bag_hist.linewidth = 0
    bag_hist.yaxis.SetTitle('Events')
    bag_hist.xaxis.SetTitle('Mass (GeV)')
    
    sig_hist.fillstyle = '0'
    sig_hist.fillcolor = 'red'
    sig_hist.linecolor = 'red'
    sig_hist.linewidth = 1
    sig_hist.yaxis.SetTitle('Events')
    sig_hist.xaxis.SetTitle('Mass (GeV)')

    test = plotter(hist=[bag_hist], sig = [sig_hist],style='CMS')
    test.Add_data(dat_hist)
    test.Add_plot('Signi',pos=0, height=15)
    test.Add_plot('Ratio',pos=1, height=15)
    test.Add_plot('DiffRatio',pos=2, height=15)
    test.make_plot('bla_plt.pdf')
    return 42

def create_test_histos():
    # set the random seed
    ROOT.gRandom.SetSeed(42)
    np.random.seed(42)

    # signal distribution
    signal = 126 + 10 * np.random.randn(1000)
    signal_obs = 126 + 10 * np.random.randn(1000)

    # create histograms
    h1 = Hist(30, 40, 200, title='Background', markersize=0)
    h2 = h1.Clone(title='Signal')
    h3 = h1.Clone(title='Data')

    # fill the histograms with our distributions
    h1.FillRandom('landau', 10000)
    map(h2.Fill, signal)
    h3.FillRandom('landau', 10000)
    map(h3.Fill, signal_obs)

    return h1, h2, h3

main()
