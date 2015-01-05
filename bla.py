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
    bag_hist.yaxis.SetTitle('Events')
    bag_hist.xaxis.SetTitle('Mass (GeV)')
    
    sig_hist.fillstyle = 'solid'
    sig_hist.fillcolor = 'red'
    sig_hist.linecolor = 'red'
    sig_hist.linewidth = 0
    sig_hist.yaxis.SetTitle('Events')
    sig_hist.xaxis.SetTitle('Mass (GeV)')

    test = plotter(hist=[bag_hist, sig_hist],style='CMS')
    test.Add_data(dat_hist)
    test.Add_signi(pos=0, height=15)
    test.Add_ratio(pos=1, height=15)
    test.Add_diff(pos=2, height=15)
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

class plotter():
    ## Constructor:
    def __init__(self, style = 'Plain', hist = [], data_hist = None, data = False, doRatio = False, doSigni = False, doDiff = False):
        ## style variables
        self.style                = style
        ## BG histograms
        self.hist                 = hist
        self.hist_height          = 100
        self.hist_start           = 0
        ## Data histograms
        self.data                 = data
        self.data_hist            = data_hist
        ## Ratio variables
        self.ratio                = doRatio
        self.ratio_height         = 20
        self.ratio_pos            = 1
        ## Significance variables
        self.signi                = doSigni
        self.signi_height         = 20
        self.signi_pos            = 0
        ## Difference variables
        self.diff                 = doDiff
        self.diff_height          = 20
        self.diff_pos             = 2
        self.annotations_modified = False
        self._Set_style()

    ##------------------------------------------------------------------
    ## Public functions
    ##------------------------------------------------------------------
    def make_plot(self,out_name):
        self._Compiler()
        self._checker()
        self._Draw()
        self._SavePlot(out_name)

#    def Modify_annotations(self, lumi = self.lumi_val, cms = self.cms_val, AddText = self.additional_text, posx = self.cms_text_x, posy = self.cms_text_y, AddAlign = self.cms_text_alignment):
#        self.annotations_modified = True
#        self.lumi_val           = lumi
#        self.cms_val            = cms
#        self.additional_text    = AddText
#        self.cms_text_x         = posx
#        self.cms_text_y         = posy
#        self.cms_text_alignment = AddAlign

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
        self._Set_style()

    def Add_signi(self, height = 20, pos = 0):
        self.signi = True
        self.signi_height = height
        self.signi_pos = pos
        self._Set_style()

    def Add_diff(self, height = 20, pos = 2):
        self.diff = True
        self.diff_height = height
        self.diff_pos = pos
        self._Set_style()

    ##------------------------------------------------------------------
    ## Private functions
    ##------------------------------------------------------------------
    def _Set_style(self):
        matplotlib.rcParams.update({'font.size': 10})
        rc('text', usetex=True)
        self.xaxis_title     = self.hist[0].xaxis.GetTitle()
        self.yaxis_title     = self.hist[0].yaxis.GetTitle()
        self.signi_text      = 'Significance'
        self.ratio_text      = 'Data/MC'
        self.diff_text       = 'Data - MC'
        self.lumi_val        = 42000
        self.cms_val         = 13
        self.additional_text = '$Preliminary$'
        self.y_label_offset  = -0.1
        if self.style == 'CMS':
            self.add_cms_text           = True
            self.add_lumi_text          = True
            self.label_text_color       = 'black'
            self.annotation_text_color  = 'black'
            self.bg_color               = 'w'
            self.ref_line_color         = 'blue'
            self.spine_color            = 'black'
            self.tick_color             = 'black'
            self.marker_style           = 'o'
            self.marker_size            = 3
            self.marker_color           = 'black'
            self.marker_error_cap_width = 0
            self.cms_text_alignment     = 'row'
            if not ((self.ratio and self.ratio_pos == 0) or (self.diff and self.diff_pos == 0) or (self.signi and self.signi_pos == 0)):
                self.cms_text_x         = 0.8
                self.cms_text_y         = 0.9
            else:
                self.cms_text_x         = 0.8
                if self.ratio_pos == 0:
                    self.cms_text_y     = 0.9 - (0.8 * self.ratio_height / 100.)
                elif self.diff_pos == 0:
                    self.cms_text_y     = 0.9 - (0.8 * self.diff_height / 100.)
                elif self.signi_pos == 0:
                    self.cms_text_y     = 0.9 - (0.8 * self.signi_height / 100.)
        elif self.style == 'Plain':
            self.add_cms_text           = False
            self.add_lumi_text          = False
            self.label_text_color       = 'black'
            self.annotation_text_color  = 'black'
            self.bg_color               = 'w'
            self.ref_line_color         = 'blue'
            self.spine_color            = 'black'
            self.tick_color             = 'black'
            self.marker_style           = 'o'
            self.marker_size            = 4
            self.marker_color           = 'black'
            self.marker_error_cap_width = 1
            self.cms_text_alignment     = 'row'
            if not ((self.ratio and self.ratio_pos == 0) or (self.diff and self.diff_pos == 0) or (self.signi and self.signi_pos == 0)):
                self.cms_text_x         = 0.8
                self.cms_text_y         = 0.9
            else:
                self.cms_text_x         = 0.8
                if self.ratio_pos == 0:
                    self.cms_text_y     = 0.9 - (0.8 * self.ratio_height / 100.)
                elif self.diff_pos == 0:
                    self.cms_text_y     = 0.9 - (0.8 * self.diff_height / 100.)
                elif self.signi_pos == 0:
                    self.cms_text_y     = 0.9 - (0.8 * self.dsigni_height / 100.)
        elif self.style == 'Cool':
            self.add_cms_text           = True
            self.add_lumi_text          = True
            self.label_text_color       = 'white'
            self.annotation_text_color  = 'white'
            self.bg_color               = '#07000d'
            self.ref_line_color         = 'y'
            self.spine_color            = '#5998ff'
            self.tick_color             = 'w'
            self.marker_style           = 'o'
            self.marker_size            = 3
            self.marker_color           = 'lightgray'
            self.marker_error_cap_width = 0
            self.cms_text_alignment     = 'column'
            self.cms_text_x             = 0.1
            self.cms_text_y             = 0.955

    def _Write_additional_text(self):
        if self.add_lumi_text:
            if self.lumi_val > 1000:
                self.fig.text(0.945, 0.955, '$%.1f\,\mathrm{fb^{-1}} (%.0f\,\mathrm{TeV})$'%(self.lumi_val/1000,self.cms_val), va='bottom', ha='right', color=self.annotation_text_color, size=12)
            else:
                self.fig.text(0.945, 0.955, '$%.0f\,\mathrm{pb^{-1}} (%.0f\,\mathrm{TeV})$'%(self.lumi_val,self.cms_val), va='bottom', ha='right', color=self.annotation_text_color, size=12)
        if self.add_cms_text:
            if self.cms_text_alignment == 'row':
                self.fig.text(self.cms_text_x, self.cms_text_y, 'CMS', va='bottom', ha='left', color=self.annotation_text_color, size=14, weight='bold')
                self.fig.text(self.cms_text_x, self.cms_text_y-0.03, self.additional_text, va='bottom', ha='left', color=self.annotation_text_color, size=10)
            elif self.cms_text_alignment == 'column':
                self.fig.text(self.cms_text_x, self.cms_text_y, 'CMS', va='bottom', ha='left', color=self.annotation_text_color, size=14, weight='bold')
                self.fig.text(self.cms_text_x + 0.08, self.cms_text_y, self.additional_text, va='bottom', ha='left', color=self.annotation_text_color, size=10)
            else:
                print('At the moment only ''row'' and ''column'' are allowed alignment values')

    def _Check_Consistency(self):
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

    def _Compiler(self):
        consi = self._Check_Consistency()
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

    def _Calc_ratio(self):
        sum_hist = self.hist[0].Clone('sum_hist')
        for i in range(1,len(self.hist)):
            sum_hist.Add(self.hist[i])
        ratio = self.data_hist.Clone('ratio')
        ratio.Divide(sum_hist)
        return ratio

    def _Calc_diff(self):
        diff = self.data_hist.Clone('diff')
        for item in self.hist:
            diff.Add(item,-1)
        return diff

    def _Calc_signi(self):
        sum_hist = self.hist[0].Clone('sum_hist')
        for i in range(1,len(self.hist)):
            sum_hist.Add(self.hist[i])
        signi = self.data_hist.Clone('signi')
        for i in range(signi.GetNbinsX()+1):
            value = float(self.data_hist.GetBinContent(i) - sum_hist.GetBinContent(i))/np.sqrt(float(pow(self.data_hist.GetBinError(i),2) + pow(sum_hist.GetBinError(i),2)))
            signi.SetBinContent(i,value)
            signi.SetBinError(i,1)
        return signi

    def _Draw(self):
        self.fig = plt.figure(figsize=(6, 6), dpi=100, facecolor=self.bg_color)

        ## Plot the main distribution on axis 1
        ax1 = plt.subplot2grid((100,1), (self.hist_start,0), rowspan=self.hist_height, colspan=1, axisbg = self.bg_color)
        if len(self.hist) == 1:
            rplt.hist(self.hist[0], stacked=False, axes=ax1, zorder=2)
            if self.data:
                rplt.errorbar(self.data_hist, xerr=False, emptybins=False, axes=ax1, 
                              markersize=self.marker_size,
                              marker = self.marker_style,
                              ecolor = self.marker_color,
                              markerfacecolor = self.marker_color,
                              markeredgecolor = self.marker_color,
                              capthick = self.marker_error_cap_width)
        else:
            rplt.hist(self.hist, stacked=True, axes=ax1, zorder=2)
            if self.data:
                rplt.errorbar(self.data_hist, xerr=False, emptybins=False, axes=ax1,
                              markersize=self.marker_size,
                              marker = self.marker_style,
                              ecolor = self.marker_color,
                              markerfacecolor = self.marker_color,
                              markeredgecolor = self.marker_color,
                              capthick = self.marker_error_cap_width)
        ax1.set_ylabel(self.yaxis_title, color=self.label_text_color, va='top', ha='left')
        ax1.yaxis.set_label_coords(self.y_label_offset,0.9)
        if not ((self.ratio and self.ratio_pos == 1) or (self.diff and self.diff_pos == 1) or (self.signi and self.signi_pos == 1) or (self.ratio and self.ratio_pos == 2) or (self.diff and self.diff_pos == 2) or (self.signi and self.signi_pos == 2)):
            plt.xlabel(self.xaxis_title, color=self.label_text_color, position=(1., -0.1), va='top', ha='right')
        ax1.yaxis.set_major_locator(mticker.MaxNLocator(prune='lower'))
        ax1.spines['bottom'].set_color(self.spine_color)
        ax1.spines['top'].set_color(self.spine_color)
        ax1.spines['left'].set_color(self.spine_color)
        ax1.spines['right'].set_color(self.spine_color)
        ax1.tick_params(axis='y', colors=self.tick_color)
        ax1.tick_params(axis='x', colors=self.tick_color)
        ## Plot a derived distribution on top of the main distribution on axis 0
        if (self.ratio and self.ratio_pos == 0) or (self.diff and self.diff_pos == 0) or (self.signi and self.signi_pos == 0):
            ax0 = plt.subplot2grid((100,1), (0,0), rowspan=self.hist_start, colspan=1, sharex = ax1, axisbg = self.bg_color)
            ax0.spines['bottom'].set_color(self.spine_color)
            ax0.spines['top'].set_color(self.spine_color)
            ax0.spines['left'].set_color(self.spine_color)
            ax0.spines['right'].set_color(self.spine_color)
            ax0.tick_params(axis='y', colors=self.tick_color)
            ax0.tick_params(axis='x', colors=self.tick_color)
            if (self.ratio and self.ratio_pos == 0):
                ratio_hist = self._Calc_ratio()
                rplt.errorbar(ratio_hist, xerr=False, emptybins=False, axes=ax0,
                              markersize=self.marker_size,
                              label=self.ratio_text,
                              marker = self.marker_style,
                              ecolor = self.marker_color,
                              markerfacecolor = self.marker_color,
                              markeredgecolor = self.marker_color,
                              capthick = self.marker_error_cap_width)
                ax0.axhline(1, color=self.ref_line_color)
                ax0.set_ylabel(self.ratio_text, color=self.label_text_color, va='top', ha='left')
                ax0.yaxis.set_label_coords(self.y_label_offset,1.)
            if (self.diff and self.diff_pos == 0):
                diff_hist = self._Calc_diff()
                rplt.errorbar(diff_hist, xerr=False, emptybins=False, axes=ax0,
                              markersize=self.marker_size,
                              label=self.diff_text,
                              marker = self.marker_style,
                              ecolor = self.marker_color,
                              markerfacecolor = self.marker_color,
                              markeredgecolor = self.marker_color,
                              capthick = self.marker_error_cap_width)
                ax0.axhline(0, color=self.ref_line_color)
                ax0.set_ylabel(self.diff_text, color=self.label_text_color, va='top', ha='left')
                ax0.yaxis.set_label_coords(self.y_label_offset,1.)
            if (self.signi and self.signi_pos == 0):
                signi_hist = self._Calc_signi()
                rplt.errorbar(signi_hist, xerr=False, emptybins=False, axes=ax0,
                              markersize=self.marker_size,
                              label=self.signi_text,
                              marker = self.marker_style,
                              ecolor = self.marker_color,
                              markerfacecolor = self.marker_color,
                              markeredgecolor = self.marker_color,
                              capthick = self.marker_error_cap_width)
                ax0.axhline(0, color=self.ref_line_color)
                ax0.set_ylabel(self.signi_text, color=self.label_text_color, va='top', ha='left')
                ax0.yaxis.set_label_coords(self.y_label_offset,1.)
            ax0.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5, prune='lower'))
            #plt.xlabel(self.xaxis_title, color=self.label_text_color, position=(1., -0.1), va='top', ha='right')
            plt.setp(ax0.get_xticklabels(), visible=False)
        ## Plot a derived distribution below the main distribution on axis 2
        if (self.ratio and self.ratio_pos == 1) or (self.diff and self.diff_pos == 1) or (self.signi and self.signi_pos == 1):
            if (self.ratio and self.ratio_pos == 1):
                ax2 = plt.subplot2grid((100,1), (self.hist_start+self.hist_height,0), rowspan=self.ratio_height, colspan=1, sharex = ax1, axisbg = self.bg_color)
                ratio_hist = self._Calc_ratio()
                rplt.errorbar(ratio_hist, xerr=False, emptybins=False, axes=ax2,
                              markersize=self.marker_size,
                              label=self.ratio_text,
                              marker = self.marker_style,
                              ecolor = self.marker_color,
                              markerfacecolor = self.marker_color,
                              markeredgecolor = self.marker_color,
                              capthick = self.marker_error_cap_width)
                ax2.axhline(1, color=self.ref_line_color)
                ax2.set_ylabel(self.ratio_text, color=self.label_text_color, va='top', ha='left')
                ax2.yaxis.set_label_coords(self.y_label_offset,1.)
            if (self.diff and self.diff_pos == 1):
                ax2 = plt.subplot2grid((100,1), (self.hist_start+self.hist_height,0), rowspan=self.diff_height, colspan=1, sharex = ax1, axisbg = self.bg_color)
                diff_hist = self._Calc_diff()
                rplt.errorbar(diff_hist, xerr=False, emptybins=False, axes=ax2,
                              markersize=self.marker_size,
                              label=self.diff_text,
                              marker = self.marker_style,
                              ecolor = self.marker_color,
                              markerfacecolor = self.marker_color,
                              markeredgecolor = self.marker_color,
                              capthick = self.marker_error_cap_width)
                ax2.axhline(0, color=self.ref_line_color)
                ax2.set_ylabel(self.diff_text, color=self.label_text_color, va='top', ha='left')
                ax2.yaxis.set_label_coords(self.y_label_offset,1.)
            if (self.signi and self.signi_pos == 1):
                ax2 = plt.subplot2grid((100,1), (self.hist_start+self.hist_height,0), rowspan=self.signi_height, colspan=1, sharex = ax1, axisbg = self.bg_color)
                signi_hist = self._Calc_signi()
                rplt.errorbar(signi_hist, xerr=False, emptybins=False, axes=ax2,
                              markersize=self.marker_size,
                              label=self.signi_text,
                              marker = self.marker_style,
                              ecolor = self.marker_color,
                              markerfacecolor = self.marker_color,
                              markeredgecolor = self.marker_color,
                              capthick = self.marker_error_cap_width)
                ax2.axhline(0, color=self.ref_line_color)
                ax2.set_ylabel(self.signi_text, color=self.label_text_color, va='top', ha='left')
                ax2.yaxis.set_label_coords(self.y_label_offset,1.)
            ax2.spines['bottom'].set_color(self.spine_color)
            ax2.spines['top'].set_color(self.spine_color)
            ax2.spines['left'].set_color(self.spine_color)
            ax2.spines['right'].set_color(self.spine_color)
            ax2.tick_params(axis='y', colors=self.tick_color)
            ax2.tick_params(axis='x', colors=self.tick_color)
            if (self.ratio and self.ratio_pos == 2) or (self.diff and self.diff_pos == 2) or (self.signi and self.signi_pos == 2):
                plt.setp(ax2.get_xticklabels(), visible=False)
                ax2.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5, prune='both'))
                plt.xlabel(self.xaxis_title, color=self.label_text_color, position=(1., -0.1), va='top', ha='right')
            else:
                ax2.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5, prune='upper'))
                plt.xlabel(self.xaxis_title, color=self.label_text_color, position=(1., -0.1), va='top', ha='right')
            plt.setp(ax1.get_xticklabels(), visible=False)
        ## Plot a derived distribution at the very bottom of the main distribution on axis 3
        if (self.ratio and self.ratio_pos == 2) or (self.diff and self.diff_pos == 2) or (self.signi and self.signi_pos == 2):
            if (self.ratio and self.ratio_pos == 2):
                ax3 = plt.subplot2grid((100,1), (100-self.ratio_height,0), rowspan=self.ratio_height, colspan=1, sharex = ax1, axisbg = self.bg_color)
                ratio_hist = self._Calc_ratio()
                rplt.errorbar(ratio_hist, xerr=False, emptybins=False, axes=ax3,
                              markersize=self.marker_size,
                              label=self.ratio_text,
                              marker = self.marker_style,
                              ecolor = self.marker_color,
                              markerfacecolor = self.marker_color,
                              markeredgecolor = self.marker_color,
                              capthick = self.marker_error_cap_width)
                ax3.axhline(1, color=self.ref_line_color)
                ax3.set_ylabel(self.ratio_text, color=self.label_text_color, va='top', ha='left')
                ax3.yaxis.set_label_coords(self.y_label_offset,1.)
            if (self.diff and self.diff_pos == 2):
                ax3 = plt.subplot2grid((100,1), (100-self.diff_height,0), rowspan=self.diff_height, colspan=1, sharex = ax1, axisbg = self.bg_color)
                diff_hist = self._Calc_diff()
                rplt.errorbar(diff_hist, xerr=False, emptybins=False, axes=ax3,
                              markersize=self.marker_size,
                              label=self.diff_text,
                              marker = self.marker_style,
                              ecolor = self.marker_color,
                              markerfacecolor = self.marker_color,
                              markeredgecolor = self.marker_color,
                              capthick = self.marker_error_cap_width)
                ax3.axhline(0, color=self.ref_line_color)
                ax3.set_ylabel(self.diff_text, color=self.label_text_color, va='top', ha='left')
                ax3.yaxis.set_label_coords(self.y_label_offset,1.)
            if (self.signi and self.signi_pos == 2):
                ax3 = plt.subplot2grid((100,1), (100-self.signi_height,0), rowspan=self.signi_height, colspan=1, sharex = ax1, axisbg = self.bg_color)
                signi_hist = self._Calc_signi()
                rplt.errorbar(signi_hist, xerr=False, emptybins=False, axes=ax3,
                              markersize=self.marker_size,
                              label=self.signi_text,
                              marker = self.marker_style,
                              ecolor = self.marker_color,
                              markerfacecolor = self.marker_color,
                              markeredgecolor = self.marker_color,
                              capthick = self.marker_error_cap_width)
                ax3.axhline(0, color=self.ref_line_color)
                ax3.set_ylabel(self.signi_text, color=self.label_text_color, va='top', ha='left')
                ax3.yaxis.set_label_coords(self.y_label_offset,1.)
            ax3.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5, prune='upper'))
            ax3.spines['bottom'].set_color(self.spine_color)
            ax3.spines['top'].set_color(self.spine_color)
            ax3.spines['left'].set_color(self.spine_color)
            ax3.spines['right'].set_color(self.spine_color)
            ax3.tick_params(axis='y', colors=self.tick_color)
            ax3.tick_params(axis='x', colors=self.tick_color)
            plt.setp(ax1.get_xticklabels(), visible=False)
            plt.xlabel(self.xaxis_title, color=self.label_text_color, position=(1., -0.1), va='top', ha='right')

        plt.subplots_adjust(left=.10, bottom=.08, right= .95, top=.95, wspace =.2, hspace=.0)
        self._Write_additional_text()

    def _SavePlot(self, out_name):
        plt.savefig(out_name,facecolor=self.fig.get_facecolor())

    def _checker(self):
        try:
            print('histo with name:' + self.hist[0].GetName())
        except AttributeError:
            print('No histogram added')
        print('with height: ' + str(self.hist_height) + ' and start: ' + str(self.hist_start))
        print('Do ratio: ' + str(self.ratio) + ' with height: ' + str(self.ratio_height) + ' and pos: ' + str(self.ratio_pos))
        print('Do signi: ' + str(self.signi) + ' with height: ' + str(self.signi_height) + ' and pos: ' + str(self.signi_pos))
        print('Do diff: ' + str(self.diff) + ' with height: ' + str(self.diff_height) + ' and pos: ' + str(self.diff_pos))

main()
