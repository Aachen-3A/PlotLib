#!/bin/env python

import ROOT
import numpy as np
from rootpy.plotting import Hist, HistStack, Legend, Canvas, Graph
import matplotlib
from matplotlib import rc
import matplotlib.ticker as mticker
import rootpy.plotting.root2matplotlib as rplt
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

def main():
    bag_hist, sig_hist, dat_hist = create_test_histos()

    bag_hist.fillstyle = 'solid'
    bag_hist.fillcolor = 'green'
    bag_hist.linecolor = 'green'
    bag_hist.linewidth = 0
    
    sig_hist.fillstyle = 'solid'
    sig_hist.fillcolor = 'red'
    sig_hist.linecolor = 'red'
    sig_hist.linewidth = 0

    test = plotter(hist=[bag_hist, sig_hist])
    test.Add_data(dat_hist)
    test.Add_signi(pos=0, height=15)
    test.Add_ratio(pos=1, height=15)
    test.Add_diff(pos=2, height=15)
    test.Compiler()
    test.checker()
    test.Draw()
    test.SavePlot()
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

class plotter():
    def __init__(self, style = 'CMS', hist = [], data_hist = None, data = False, doRatio = False, doSigni = False, doDiff = False):
        self.style = style
        self.hist = hist
        self.hist_height = 100
        self.hist_start = 0
        self.data = data
        self.data_hist = data_hist
        self.ratio = doRatio
        self.ratio_height = 20
        self.ratio_pos = 1
        self.signi = doSigni
        self.signi_height = 20
        self.signi_pos = 0
        self.diff = doDiff
        self.diff_height = 20
        self.diff_pos = 2

    def Add_data(self, data_hist):
        self.data = True
        self.data_hist = data_hist

    def Draw_data(self, doData = True):
        self.data = doData

    def Add_histo(self, histo):
        self.hist.append(histo)

    def Add_ratio(self, height = 20, pos = 1):
        self.ratio = True
        self.ratio_height = height
        self.ratio_pos = pos

    def Add_signi(self, height = 20, pos = 0):
        self.signi = True
        self.signi_height = height
        self.signi_pos = pos

    def Add_diff(self, height = 20, pos = 2):
        self.diff = True
        self.diff_height = height
        self.diff_pos = pos

    def Check_Consistency(self):
        if self.ratio and self.signi and not self.diff:
            if self.ratio_pos == self.signi_pos:
                print('Can not put Ratio and Significance plot at the same position, sorry!')
                return False
        if self.ratio and not self.signi and self.diff:
            if self.ratio_pos == self.diff_pos:
                print('Can not put Ratio and Difference plot at the same position, sorry!')
                return False
        if not self.ratio and self.signi and self.diff:
            if self.diff_pos == self.signi_pos:
                print('Can not put Difference and Significance plot at the same position, sorry!')
                return False
        if self.ratio and self.signi and self.diff:
            if self.ratio_pos == self.signi_pos:
                print('Can not put Ratio and Significance plot at the same position, sorry!')
                return False
            if self.ratio_pos == self.diff_pos:
                print('Can not put Ratio and Difference plot at the same position, sorry!')
                return False
            if self.diff_pos == self.signi_pos:
                print('Can not put Difference and Significance plot at the same position, sorry!')
                return False
        return True

    def Compiler(self):
        consi = self.Check_Consistency()
        if consi:
            if self.ratio and self.ratio_pos == 0:
                self.hist_start = self.ratio_height
                self.hist_height -= self.ratio_height
            elif self.diff and self.diff_pos == 0:
                self.hist_start = self.diff_height
                self.hist_height -= self.diff_height
            elif self.signi and self.signi_pos == 0:
                self.hist_start = self.signi_height
                self.hist_height -= self.signi_height
            if self.ratio and self.ratio_pos == 1:
                self.hist_height -= self.ratio_height
            elif self.diff and self.diff_pos == 1:
                self.hist_height -= self.diff_height
            elif self.signi and self.signi_pos == 1:
                self.hist_height -= self.signi_height
            if self.ratio and self.ratio_pos == 2:
                self.hist_height -= self.ratio_height
            elif self.diff and self.diff_pos == 2:
                self.hist_height -= self.diff_height
            elif self.signi and self.signi_pos == 2:
                self.hist_height -= self.signi_height
        else:
            print('Can not cpmpile the plot: see previous error messages')

    def Calc_ratio(self):
        sum_hist = self.hist[0].Clone('sum_hist')
        for i in range(1,len(self.hist)):
            sum_hist.Add(self.hist[i])
        ratio = self.data_hist.Clone('ratio')
        ratio.Divide(sum_hist)
        return ratio

    def Calc_diff(self):
        diff = self.data_hist.Clone('diff')
        for item in self.hist:
            diff.Add(item,-1)
        return diff

    def Calc_signi(self):
        sum_hist = self.hist[0].Clone('sum_hist')
        for i in range(1,len(self.hist)):
            sum_hist.Add(self.hist[i])
        signi = self.data_hist.Clone('signi')
        for i in range(signi.GetNbinsX()+1):
            value = float(self.data_hist.GetBinContent(i) - sum_hist.GetBinContent(i))/np.sqrt(float(pow(self.data_hist.GetBinError(i),2) + pow(sum_hist.GetBinError(i),2)))
            signi.SetBinContent(i,value)
            signi.SetBinError(i,1)
        return signi

    def Set_Style(self):
        matplotlib.rcParams.update({'font.size': 10})
        matplotlib.rcParams.update({'text': True})
        #rc('text', usetex=True)

    def Draw(self):
        self.Set_Style()
        self.fig = plt.figure(figsize=(6, 6), dpi=100, facecolor='w')

        ## Plot the main distribution on axis 1
        ax1 = plt.subplot2grid((100,1), (self.hist_start,0), rowspan=self.hist_height, colspan=1)
        if len(self.hist) == 1:
            rplt.bar(self.hist[0], stacked=False, axes=ax1)
            if self.data:
                rplt.errorbar(self.data_hist, xerr=False, emptybins=False, axes=ax1, markersize=4)
        else:
            stack = HistStack()
            for item in self.hist:
                stack.Add(item)
            rplt.bar(stack, stacked=True, axes=ax1)
            if self.data:
                rplt.errorbar(self.data_hist, xerr=False, emptybins=False, axes=ax1, markersize=4)
        ax1.yaxis.set_major_locator(mticker.MaxNLocator(prune='lower'))
        ## Plot a derived distribution on top of the main distribution on axis 0
        if (self.ratio and self.ratio_pos == 0) or (self.diff and self.diff_pos == 0) or (self.signi and self.signi_pos == 0):
            ax0 = plt.subplot2grid((100,1), (0,0), rowspan=self.hist_start, colspan=1, sharex = ax1)
            if (self.ratio and self.ratio_pos == 0):
                ratio_hist = self.Calc_ratio()
                rplt.errorbar(ratio_hist, xerr=False, emptybins=False, axes=ax0, markersize=4, label='Data/MC')
                ax0.axhline(1, color='blue')
            if (self.diff and self.diff_pos == 0):
                diff_hist = self.Calc_diff()
                rplt.errorbar(diff_hist, xerr=False, emptybins=False, axes=ax0, markersize=4, label='Data - MC')
                ax0.axhline(0, color='blue')
            if (self.signi and self.signi_pos == 0):
                signi_hist = self.Calc_signi()
                rplt.errorbar(signi_hist, xerr=False, emptybins=False, axes=ax0, markersize=4, label='Significance')
                ax0.axhline(0, color='blue')
            ax0.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5, prune='lower'))
            plt.setp(ax0.get_xticklabels(), visible=False)
        ## Plot a derived distribution below the main distribution on axis 2
        if (self.ratio and self.ratio_pos == 1) or (self.diff and self.diff_pos == 1) or (self.signi and self.signi_pos == 1):
            if (self.ratio and self.ratio_pos == 1):
                ax2 = plt.subplot2grid((100,1), (self.hist_start+self.hist_height,0), rowspan=self.ratio_height, colspan=1, sharex = ax1)
                ratio_hist = self.Calc_ratio()
                rplt.errorbar(ratio_hist, xerr=False, emptybins=False, axes=ax2, markersize=4, label='Data/MC')
                ax2.axhline(1, color='blue')
            if (self.diff and self.diff_pos == 1):
                ax2 = plt.subplot2grid((100,1), (self.hist_start+self.hist_height,0), rowspan=self.diff_height, colspan=1, sharex = ax1)
                diff_hist = self.Calc_diff()
                rplt.errorbar(diff_hist, xerr=False, emptybins=False, axes=ax2, markersize=4, label='Data - MC')
                ax2.axhline(0, color='blue')
            if (self.signi and self.signi_pos == 1):
                ax2 = plt.subplot2grid((100,1), (self.hist_start+self.hist_height,0), rowspan=self.signi_height, colspan=1, sharex = ax1)
                signi_hist = self.Calc_signi()
                rplt.errorbar(signi_hist, xerr=False, emptybins=False, axes=ax2, markersize=4, label='Significance')
                ax2.axhline(0, color='blue')
            plt.setp(ax1.get_xticklabels(), visible=False)
            if (self.ratio and self.ratio_pos == 2) or (self.diff and self.diff_pos == 2) or (self.signi and self.signi_pos == 2):
                plt.setp(ax2.get_xticklabels(), visible=False)
                ax2.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5, prune='both'))
            else:
                ax2.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5, prune='upper'))
        ## Plot a derived distribution at the very bottom of the main distribution on axis 3
        if (self.ratio and self.ratio_pos == 2) or (self.diff and self.diff_pos == 2) or (self.signi and self.signi_pos == 2):
            if (self.ratio and self.ratio_pos == 2):
                ax3 = plt.subplot2grid((100,1), (100-self.ratio_height,0), rowspan=self.ratio_height, colspan=1, sharex = ax1)
                ratio_hist = self.Calc_ratio()
                rplt.errorbar(ratio_hist, xerr=False, emptybins=False, axes=ax3, markersize=4, label='Data/MC')
                ax3.axhline(1, color='blue')
            if (self.diff and self.diff_pos == 2):
                ax3 = plt.subplot2grid((100,1), (100-self.diff_height,0), rowspan=self.diff_height, colspan=1, sharex = ax1)
                diff_hist = self.Calc_diff()
                rplt.errorbar(diff_hist, xerr=False, emptybins=False, axes=ax3, markersize=4, label='Data - MC')
                ax3.axhline(0, color='blue')
            if (self.signi and self.signi_pos == 2):
                ax3 = plt.subplot2grid((100,1), (100-self.signi_height,0), rowspan=self.signi_height, colspan=1, sharex = ax1)
                signi_hist = self.Calc_signi()
                rplt.errorbar(signi_hist, xerr=False, emptybins=False, axes=ax3, markersize=4, label='Significance')
                ax3.axhline(0, color='blue')
            ax3.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5, prune='upper'))

        plt.subplots_adjust(left=.08, bottom=.10, right= .92, top=.95, wspace =.2, hspace=.0)

    def SavePlot(self):
        plt.savefig('bla_plt.pdf',facecolor=self.fig.get_facecolor())

    def checker(self):
        try:
            print('histo with name:' + self.hist[0].GetName())
        except AttributeError:
            print('No histogram added')
        print('with height: ' + str(self.hist_height) + ' and start: ' + str(self.hist_start))
        print('Do ratio: ' + str(self.ratio) + ' with height: ' + str(self.ratio_height) + ' and pos: ' + str(self.ratio_pos))
        print('Do signi: ' + str(self.signi) + ' with height: ' + str(self.signi_height) + ' and pos: ' + str(self.signi_pos))
        print('Do diff: ' + str(self.diff) + ' with height: ' + str(self.diff_height) + ' and pos: ' + str(self.diff_pos))

main()
