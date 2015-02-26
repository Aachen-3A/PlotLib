#!/bin/env python

import sys
import matplotlib
from matplotlib import rc

matplotlib.use('Qt4Agg')

import ROOT
import subprocess
import numpy as np
from rootpy.plotting import Hist, HistStack, Legend, Canvas, Graph, Pad
import matplotlib.ticker as mticker
import rootpy.plotting.root2matplotlib as rplt
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from plotlib import duke_errorbar
from operator import methodcaller
from rounding import rounding


import style_class as sc

##@class plotter
# Class to collect matplotlib functions
#
# To use the various matplotlib functions to produce the standard
# plots (with data[optional], signal[optional], background and
# uncertainties[optional]). Also different analysis distributions
# like ratio or siginficance can be added.
#
# @TODO Include handling of overflow bins
# @TODO Include cumulative distributions
# @TODO Functionallity for rebinning/variable binning
# @TODO Handling of asymetric errors (systematics)
# @TODO Handling of the data error bars
# @TODO Include file reading functionallity
# @TODO Include hist reweighting
#
# written by Soeren Erdweg 2014-2015
class plotter():
    ## Init function
    #
    # In this function the default variables are set. Also the style can
    # be defined and the histogram input can be given.
    # @param[in] hist List of background histograms (default = [])
    # @param[in] sig List of signal histograms (default = [])
    # @param[in] data_hist Data histogram (default = None)
    # @param[in] data Bool if data should be plotted (default = False)
    # @param[in] style Style container that should be used for the plot
    # @param[in] kwargs dict of key word arguments that will be passed to style
    def __init__(self, hist = [], sig = [], hist_axis = [], data_hist = None, data = False, style = sc.style_container(), **kwargs):
        ## style variables
        # self._style                = style
        ## BG histograms
        self._hist                 = hist
        self._hist_height          = 100
        self._hist_start           = 0
        ## SG histograms
        self._sig_hist             = sig
        ## Data histograms
        self._data                 = data
        self._data_hist            = data_hist
        ## Second axis histograms
        self._hist_axis            = hist_axis
        ## Additional plots
        self._add_plots            = ['', '', '']
        self._add_plots_height     = [0, 0, 0]
        self._add_plots_labels     = ['', '', '']
        self._add_plots_ref_line   = [0, 0, 0]
        self._annotations_modified = False
        self._add_error_bands      = False
        self._error_hist           = []
        self._fig                  = None
        self._allHists=self._hist+self._sig_hist+[self._data_hist]+self._hist_axis
        self._Style_cont = style
        self._useRoot = self._Style_cont.Get_useRoot()
        self._Style_cont.AddAxisTitle(self._allHists[0])
        if self._useRoot:
            ROOT.gROOT.SetBatch()
        if len(self._hist_axis) > 0:
            self._Style_cont.AddAxisTitle_histaxis(self._hist_axis[0])
            self._Style_cont.InitStyle(histaxis = self._hist_axis)
        else:
            self._Style_cont.InitStyle()
        self.additionalPad          = []

    ## del function
    #
    # This deletes the main objects nedded to not get a crash at the end!
    def __del__(self):
        plt.close()
        del self._hist
        del self._sig_hist
        del self._data_hist
        del self._hist_axis
        del self._fig


    ##------------------------------------------------------------------
    ## Public functions
    ##------------------------------------------------------------------
    ## Function to make the complete plot, after all definitions are set
    #
    # This function calls the different sub functions used to produce the
    # final plots and save it.
    # @param[in] out_name Name of the output file that should be produced
    def make_plot(self,out_name):
        self._Compiler()
        self._checker()
        self._Draw()
        self._SavePlot(out_name)

    ## Function to create the complete plot, after all definitions are set
    #
    # This function calls the different sub functions used to produce the
    # final plotsbut does not save it.
    # @param[out] _fig Created plot, to do your own custemization
    def create_plot(self):
        self._Compiler()
        self._checker()
        self._Draw()
        return self._fig

#    def Modify_annotations(self, lumi = self.lumi_val, cms = self.cms_val, AddText = self.additional_text, posx = self.cms_text_x, posy = self.cms_text_y, AddAlign = self.cms_text_alignment):
#        self.annotations_modified = True
#        self.lumi_val           = lumi
#        self.cms_val            = cms
#        self.additional_text    = AddText
#        self.cms_text_x         = posx
#        self.cms_text_y         = posy
#        self.cms_text_alignment = AddAlign

    ## Function to add the data histogram
    #
    # This function is used to add a data histogram and set the bool to
    # plot the data
    # @param[in] data_hist Data histogram that should be added
    # @param[in] doData Boolean if the data should be drawn (default = True)
    def Add_data(self, data_hist, doData = True):
        self._data = doData
        self._data_hist = data_hist

    ## Function to decide if you want to draw the data
    #
    # This function is used to set the bool to plot the data
    # @param[in] doData Boolean if the data should be drawn (default = True)
    def Draw_data(self, doData = True):
        self._data = doData

    ## Function to add a analysis plot to the figure
    #
    # This function is called to add an additional plot to the figure and
    # define its properties, like where it should be placed and how much
    # space of the figure should be taken by this plot.
    # At the moment 'Ratio', 'Diff', 'Signi', 'DiffRatio' and 'SoverSplusB' are available
    # as additional plots.
    # @param[in] plot String of the plot name that should be added (default = 'Ratio')
    # @param[in] pos Position where the plot should be added, 0 is on top of the main plot, 1 and 2 at the bottom (default = 0)
    # @param[in] height Height of the Plot from the whole plotting range in percent (default = 15)
    # @param[in] label Label of the y-axis for this additional plot (default = ''[Use the default of this specific analysis plot])
    def Add_plot(self, plot = 'Ratio', pos = 0, height = 15, label = ''):
        if self._add_plots[pos] == '':
            self._add_plots[pos] = plot
            self._add_plots_height[pos] = height
            if label == '':
                self._add_plots_labels[pos] = plot
            else:
                self._add_plots_labels[pos] = label
            self._Style_cont.InitStyle(addplots = self._add_plots, addheights = self._add_plots_height)
        else:
            print('\n\tfor pos %.0f is already %s planned, so that is not possible\t'%(pos,self.add_plots[pos]))
            sys.exit(42)

    ## Function to add a histogram to the background list
    #
    # This function is used to add an additional histogram to the list of
    # background histogram.
    # @param[in] histo Histogram that should be added
    def Add_histo(self, histo):
        self._hist.append(histo)

    ## Function to set the uncertainty histograms
    #
    # This function is used to add the systematic uncertainty histograms.
    # @param[in] histo List of histograms that contain as bin content the relativ systematic uncertainties.
    # @param[in] labels List of labels for the systematic uncertainties.
    # @param[in] band_center Parameter where the error band should be centered ('ref', at the reference line,
    #                        or 'val' around the e.g. ratio value) (default = 'ref')
    # @param[in] stacking String to identify how to stack different systematic uncertainties ('No' stacking,
    #                     'linear' stacking) (Default = 'No')
    def Add_error_hist(self, histo = [], labels = [], band_center = 'ref', stacking = 'No'):
        self._add_error_bands = True
        self._error_hist = histo
        if labels != []:
            self._Style_cont.Set_error_bands_labl(labels)
        self._Style_cont.Set_error_bands_center(band_center)
        self._Style_cont.Set_error_stacking(stacking)

    ## Function to set properties of the plotting axis
    #
    # This function sets axis properties like the y-range or
    # if any axis should be logarithmic.
    # @param[in] logx Boolean if the x-axis should be logarithmic (Default = False)
    # @param[in] logy Boolean if the y-axis should be logarithmic (Default = True)
    # @param[in] ymin Minimum plotting range for the y-axis (Default = -1 automatic values)
    # @param[in] ymax Maximum plotting range for the y-axis (Default = -1 automatic values)
    # @param[in] xmin Minimum plotting range for the x-axis (Default = -1 range from hist)
    # @param[in] xmax Maximum plotting range for the x-axis (Default = -1 range from hist)
    def Set_axis(self, logx = False, logy = True, ymin = -1, ymax = -1, xmin = -1, xmax = -1, grid = False):
        self._Style_cont.Set_axis(logx = logx, logy = logy, ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax, grid = grid)

    ## Function to show the complete plot in the matplotlib browser
    #
    # This function shows the plot in the matplotlib browser, so that the
    # user can modify it.
    def show_fig(self):
        self._fig.show()
        raw_input('bla')

    ## Function to save the complete plot
    #
    # This function saves the plot you which is stored in the object so create it first
    # @param[in] out_name name of the output file
    def SavePlot(self, out_name):
        self._SavePlot(out_name)

    def ChangeStyle(self,**kwargs):
        for key in kwargs:
            if hasattr(self._Style_cont,"_"+key):
                setattr(self._Style_cont,"_"+key,kwargs[key])
            else:
                print "\n\tNo attribute _%s in syle\n"%key
                sys.exit(42)

    ## Function to get the axis 0
    #
    # This function returns the axis for the top plot
    # @param[out] self._ax0 axis 0
    def Get_axis0(self):
        try:
            return self._ax0
        except(AttributeError):
            print("\n\tThe axis 0, you want to get does not exist\n")
            sys.exit(42)

    ## Function to get the axis 1
    #
    # This function returns the axis for the main plot
    # @param[out] self._ax1 axis 1
    def Get_axis1(self):
        try:
            return self._ax1
        except(AttributeError):
            print("\n\tThe axis 1, you want to get does not exist\n")
            sys.exit(42)

    ## Function to get the axis 2
    #
    # This function returns the axis for the plot below the main plot
    # @param[out] self._ax2 axis 2
    def Get_axis2(self):
        try:
            return self._ax2
        except(AttributeError):
            print("\n\tThe axis 2, you want to get does not exist\n")
            sys.exit(42)

    ## Function to get the axis 3
    #
    # This function returns the axis for the lowest plot
    # @param[out] self._ax3 axis 3
    def Get_axis3(self):
        try:
            return self._ax3
        except(AttributeError):
            print("\n\tThe axis 3, you want to get does not exist\n")
            sys.exit(42)

    ##------------------------------------------------------------------
    ## Private functions
    ##------------------------------------------------------------------
    def _Write_additional_text(self):
        if self._Style_cont.Get_add_lumi_text():
            self._Style_cont.Set_lumi_val(float(self._Style_cont.Get_lumi_val()))
            if self._Style_cont.Get_lumi_val() >= 1000:
                if len(self._hist_axis) > 0:
                    self._fig.text(0.915, 0.955, '$%.1f\,\mathrm{fb^{-1}} (%.0f\,\mathrm{TeV})$'%(self._Style_cont.Get_lumi_val()/1000,self._Style_cont.Get_cms_val()), va='bottom', ha='right', color=self._Style_cont.Get_annotation_text_color(), size=12)
                else:
                    self._fig.text(0.945, 0.955, '$%.1f\,\mathrm{fb^{-1}} (%.0f\,\mathrm{TeV})$'%(self._Style_cont.Get_lumi_val()/1000,self._Style_cont.Get_cms_val()), va='bottom', ha='right', color=self._Style_cont.Get_annotation_text_color(), size=12)
            else:
                if len(self._hist_axis) > 0:
                    self._fig.text(0.915, 0.955, '$%.0f\,\mathrm{pb^{-1}} (%.0f\,\mathrm{TeV})$'%(self._Style_cont.Get_lumi_val(),self._Style_cont.Get_cms_val()), va='bottom', ha='right', color=self._Style_cont.Get_annotation_text_color(), size=12)
                else:
                    self._fig.text(0.945, 0.955, '$%.0f\,\mathrm{pb^{-1}} (%.0f\,\mathrm{TeV})$'%(self._Style_cont.Get_lumi_val(),self._Style_cont.Get_cms_val()), va='bottom', ha='right', color=self._Style_cont.Get_annotation_text_color(), size=12)
        if self._Style_cont.Get_add_cms_text():
            if self._Style_cont.Get_cms_text_alignment() == 'row':
                self._fig.text(self._Style_cont.Get_cmsTextPosition().getX(), self._Style_cont.Get_cmsTextPosition().getY(), 'CMS', va='bottom', ha='left', color=self._Style_cont.Get_annotation_text_color(), size=14, weight='bold')
                self._fig.text(self._Style_cont.Get_cmsTextPosition().getX(), self._Style_cont.Get_cmsTextPosition().getY()-0.03, self._Style_cont.Get_additional_text(), va='bottom', ha='left', color=self._Style_cont.Get_annotation_text_color(), size=10, style = 'italic')
            elif self._Style_cont.Get_cms_text_alignment() == 'column':
                self._fig.text(self._Style_cont.Get_cmsTextPosition().getX(), self._Style_cont.Get_cmsTextPosition().getY(), 'CMS', va='bottom', ha='left', color=self._Style_cont.Get_annotation_text_color(), size=14, weight='bold')
                self._fig.text(self._Style_cont.Get_cmsTextPosition().getX() + 0.08, self._Style_cont.Get_cmsTextPosition().getY(), self._Style_cont.Get_additional_text(), va='bottom', ha='left', color=self._Style_cont.Get_annotation_text_color(), size=10, style = 'italic')
            else:
                print('At the moment only ''row'' and ''column'' are allowed alignment values')

    def _Add_legend(self):
        if self._Style_cont.Get_no_legend():
            return
        if self._add_plots[0] != '':
            self._Style_cont.Get_LegendPosition().addYspace(-(0.85 * self._add_plots_height[0] / 100.))
        if self._add_plots[1] != '':
            self._Style_cont.Get_LegendPosition().addYspace(  0.8 * self._add_plots_height[1] / 100.)
        if self._add_plots[2] != '':
            self._Style_cont.Get_LegendPosition().addYspace(  0.8 * self._add_plots_height[2] / 100.)
        if len(self._hist_axis) > 0:
            self._Style_cont.Get_LegendPosition().addXspace(  -0.04  )

        if self._Style_cont.Get_LegendPosition() == self._Style_cont.Get_cmsTextPosition():
            self._Style_cont.Get_LegendPosition().addYspace(self._Style_cont.Get_cmsTextPosition().getY()-self._Style_cont.Get_LegendPosition().getY()-0.02)
        handle_list = []
        label_list = []
        if self._Style_cont.Get_kind() == 'Standard':
            for item in self._hist:
                col_patch = mpatches.Patch(color = item.GetFillColor())
                handle_list.append(col_patch)
                label_list.append(item.GetTitle())
            for item in self._sig_hist:
                col_patch = mlines.Line2D([], [], color = item.GetLineColor(), markersize = 0)
                handle_list.append(col_patch)
                label_list.append(item.GetTitle())
            if self._add_error_bands:
                for i in range(0,len(self._error_hist)):
                    col_patch = mpatches.Patch(facecolor = self._Style_cont.Get_error_bands_fcol()[i], edgecolor = self._Style_cont.Get_error_bands_ecol()[i] , alpha = self._Style_cont.Get_error_bands_alph(), lw = 0.7)
                    handle_list.append(col_patch)
                    label_list.append(self._Style_cont.Get_error_bands_labl()[i])
                if self._Style_cont.Get_error_stacking() == 'No':
                    col_patch = mpatches.Patch(facecolor = 'grey', edgecolor = 'black' , alpha = 0.4 , lw = 0.7)
                    handle_list.append(col_patch)
                    label_list.append('syst. sum')
            for item in self._hist_axis:
                col_patch = mlines.Line2D([], [], color = item.GetLineColor(), markersize = 0)
                handle_list.append(col_patch)
                label_list.append(item.GetTitle())
            if self._data:
                dat_line=plt.errorbar([], [],xerr = False,yerr=True, markersize = self._Style_cont.Get_marker_size(),
                                  marker = self._Style_cont.Get_marker_style(),
                                  color = self._Style_cont.Get_marker_color(),
                                  capthick = self._Style_cont.Get_marker_error_cap_width())
                handle_list.append(dat_line)
                label_list.append(self._data_hist.GetTitle())
        elif self._Style_cont.Get_kind() == 'Lines':
            for item in self._hist:
                col_patch = mlines.Line2D([], [], color = item.GetLineColor(), markersize = 0)
                handle_list.append(col_patch)
                label_list.append(item.GetTitle())
            for item in self._sig_hist:
                col_patch = mlines.Line2D([], [], color = item.GetLineColor(), markersize = 0)
                handle_list.append(col_patch)
                label_list.append(item.GetTitle())
            for item in self._hist_axis:
                col_patch = mlines.Line2D([], [], color = item.GetLineColor(), markersize = 0)
                handle_list.append(col_patch)
                label_list.append(item.GetTitle())
        elif self._Style_cont.Get_kind() == 'Graphs':
            for item in self._hist:
                dat_line=plt.errorbar([], [],xerr = False,yerr=True, markersize = self._Style_cont.Get_marker_size(),
                                  marker = self._Style_cont.Get_marker_style(),
                                  color = item.GetLineColor(),
                                  capthick = self._Style_cont.Get_marker_error_cap_width())
                handle_list.append(dat_line)
                label_list.append(item.GetTitle())
            for item in self._sig_hist:
                dat_line=plt.errorbar([], [],xerr = False,yerr=True, markersize = self._Style_cont.Get_marker_size(),
                                  marker = self._Style_cont.Get_marker_style(),
                                  color = item.GetLineColor(),
                                  capthick = self._Style_cont.Get_marker_error_cap_width())
                handle_list.append(dat_line)
                label_list.append(item.GetTitle())
            for item in self._hist_axis:
                dat_line=plt.errorbar([], [],xerr = False,yerr=True, markersize = self._Style_cont.Get_marker_size(),
                                  marker = self._Style_cont.Get_marker_style(),
                                  color = item.GetLineColor(),
                                  capthick = self._Style_cont.Get_marker_error_cap_width())
                handle_list.append(dat_line)
                label_list.append(item.GetTitle())
            if self._add_error_bands:
                for i in range(0,len(self._error_hist)):
                    col_patch = mpatches.Patch(facecolor = self._Style_cont.Get_error_bands_fcol()[i], edgecolor = self._Style_cont.Get_error_bands_ecol()[i] , alpha = self._Style_cont.Get_error_bands_alph(), lw = 0.7)
                    handle_list.append(col_patch)
                    label_list.append(self._Style_cont.Get_error_bands_labl()[i])
                if self._Style_cont.Get_error_stacking() == 'No':
                    col_patch = mpatches.Patch(facecolor = 'grey', edgecolor = 'black' , alpha = 0.4 , lw = 0.7)
                    handle_list.append(col_patch)
                    label_list.append('syst. sum')
            if self._data:
                dat_line=plt.errorbar([], [],xerr = False,yerr=True, markersize = self._Style_cont.Get_marker_size(),
                                  marker = self._Style_cont.Get_marker_style(),
                                  color = self._Style_cont.Get_marker_color(),
                                  capthick = self._Style_cont.Get_marker_error_cap_width())
                handle_list.append(dat_line)
                label_list.append(self._data_hist.GetTitle())
        self.leg = plt.legend(handle_list, label_list,
                    loc = self._Style_cont.Get_LegendPosition().get_positiontext(),
                    bbox_to_anchor=(self._Style_cont.Get_LegendPosition().getX(),self._Style_cont.Get_LegendPosition().getY()),
                    bbox_transform=plt.gcf().transFigure,
                    numpoints = 1,
                    frameon = False,
                    fontsize = self._Style_cont.Get_legend_font_size())

        for text in self.leg.get_texts():
            text.set_color(self._Style_cont.Get_annotation_text_color())

    def _Compiler(self):
        if len(self._hist) == 0:
            self._add_plots[0] = ''
            self._add_plots[1] = ''
            self._add_plots[2] = ''
            self._add_error_bands = False
        if self._Style_cont.Get_kind() == 'Lines':
            self._add_plots[0] = ''
            self._add_plots[1] = ''
            self._add_plots[2] = ''
        if self._Style_cont.Get_kind() == 'Graphs':
            self._add_plots[0] = ''
            self._add_plots[1] = ''
            self._add_plots[2] = ''
            self._add_error_bands = False
        if self._add_plots[0] != '':
            self._hist_start = self._add_plots_height[0]
            self._hist_height -= self._add_plots_height[0]
        if self._add_plots[1] != '':
            self._hist_height -= self._add_plots_height[1]
        if self._add_plots[2] != '':
            self._hist_height -= self._add_plots_height[2]
        # sort syst hist by the integral
        self._error_hist = sorted(self._error_hist, key=methodcaller('Integral'), reverse=True)
        #matplotlib draws no errorbars in logy when the lower error = 0
        if self._Style_cont.Get_logy() and self._data:
            for ibin in self._data_hist.bins():
                if ibin.error==1:
                    ibin.error=1.-1e-12

    def cleanUnwantedBins(self,hist,toCleanHists):
        #remove=[]
        if toCleanHists is not None:
            for i in hist.bins():
                ignore=False
                for ihist in toCleanHists:
                    if ihist[i.idx].value==0:
                        ignore=True
                if ignore:
                    i.value=0
                    i.error=0

    def _Calc_additional_plot(self, plot, pos):
        if plot == 'Ratio':
            self._add_plots_labels[pos] = 'Data/MC'
            self._add_plots_ref_line[pos] = 1.
            return self._Calc_ratio()
        elif plot == 'Diff':
            self._add_plots_labels[pos] = 'Data - MC'
            self._add_plots_ref_line[pos] = 0.
            return self._Calc_diff()
        elif plot == 'Signi':
            self._add_plots_labels[pos] = 'Significance'
            self._add_plots_ref_line[pos] = 0.
            return self._Calc_signi()
        elif plot == 'DiffRatio':
            self._add_plots_labels[pos] = '$\mathdefault{\\frac{Data - MC}{MC}}$'
            self._add_plots_ref_line[pos] = 0.
            return self._Calc_diffratio()
        elif plot == 'SoverSplusB':
            self._add_plots_labels[pos] = '$\mathdefault{\\frac{Signal}{\sqrt{Signal + MC}}}$'
            self._add_plots_ref_line[pos] = 0.
            return self._Calc_SoverSpB()
        else:
            print('%s is not implemented yet as an additional plot, feel free to include this functionallity')

    def _Calc_ratio(self):
        sum_hist = self._hist[0].Clone('sum_hist')
        for i in range(1,len(self._hist)):
            sum_hist.Add(self._hist[i])
        ratio = self._data_hist.Clone('ratio')
        ratio.Divide(sum_hist)
        self.cleanUnwantedBins(ratio,[sum_hist,self._data_hist])
        x = []
        y = []
        err = []
        for j in range(0,len(self._error_hist)):
            x_i = []
            y_i = []
            err_i = []
            for i in range(sum_hist.GetNbinsX()+1):
                x_i.append(sum_hist.GetBinLowEdge(i))
                x_i.append(sum_hist.GetBinLowEdge(i) + sum_hist.GetBinWidth(i))
                if self._Style_cont.Get_error_bands_center() == 'ref':
                    y_i.append(1.)
                    y_i.append(1.)
                elif self._Style_cont.Get_error_bands_center() == 'val':
                    y_i.append(ratio.GetBinContent(i))
                    y_i.append(ratio.GetBinContent(i))
                err_i.append(ratio.GetBinContent(i) * self._error_hist[j].GetBinContent(i))
                err_i.append(ratio.GetBinContent(i) * self._error_hist[j].GetBinContent(i))
            x.append(np.array(x_i))
            y.append(np.array(y_i))
            err.append(np.array(err_i))
        return ratio, x, y, err

    def _Calc_diff(self):
        sum_hist = self._hist[0].Clone('sum_hist')
        for i in range(1,len(self._hist)):
            sum_hist.Add(self._hist[i])
        diff = self._data_hist.Clone('diff')
        diff.Add(sum_hist,-1)
        self.cleanUnwantedBins(diff,[sum_hist,self._data_hist])
        x = []
        y = []
        err = []
        for j in range(0,len(self._error_hist)):
            x_i = []
            y_i = []
            err_i = []
            for i in range(sum_hist.GetNbinsX()+1):
                x_i.append(sum_hist.GetBinLowEdge(i))
                x_i.append(sum_hist.GetBinLowEdge(i) + sum_hist.GetBinWidth(i))
                if self._Style_cont.Get_error_bands_center() == 'ref':
                    y_i.append(0.)
                    y_i.append(0.)
                elif self._Style_cont.Get_error_bands_center() == 'val':
                    y_i.append(diff.GetBinContent(i))
                    y_i.append(diff.GetBinContent(i))
                err_i.append(sum_hist.GetBinContent(i) * self._error_hist[j].GetBinContent(i))
                err_i.append(sum_hist.GetBinContent(i) * self._error_hist[j].GetBinContent(i))
            x.append(np.array(x_i))
            y.append(np.array(y_i))
            err.append(np.array(err_i))
        return diff, x, y, err

    def _Calc_diffratio(self):
        diff = self._data_hist.Clone('diffratio')

        sum_hist = sum(self._hist)
        diff.Add(sum_hist,-1)
        diff.Divide(sum_hist)
        self.cleanUnwantedBins(diff,[sum_hist,self._data_hist])

        x = []
        y = []
        err = []
        for j in range(0,len(self._error_hist)):
            x_i = []
            y_i = []
            err_i = []
            for i in range(sum_hist.GetNbinsX()+1):
                x_i.append(sum_hist.GetBinLowEdge(i))
                x_i.append(sum_hist.GetBinLowEdge(i) + sum_hist.GetBinWidth(i))
                if self._Style_cont.Get_error_bands_center() == 'ref':
                    y_i.append(0.)
                    y_i.append(0.)
                elif self._Style_cont.Get_error_bands_center() == 'val':
                    y_i.append(diff.GetBinContent(i))
                    y_i.append(diff.GetBinContent(i))
                if sum_hist.GetBinContent(i) > 0:
                    err_i.append(self._data_hist.GetBinContent(i) / sum_hist.GetBinContent(i) * self._error_hist[j].GetBinContent(i))
                    err_i.append(self._data_hist.GetBinContent(i) / sum_hist.GetBinContent(i) * self._error_hist[j].GetBinContent(i))
                else:
                    err_i.append(0)
                    err_i.append(0)

            x.append(np.array(x_i))
            y.append(np.array(y_i))
            err.append(np.array(err_i))
        return diff, x, y, err

    def _Calc_signi(self):
        sum_hist = self._hist[0].Clone('sum_hist')
        for i in range(1,len(self._hist)):
            sum_hist.Add(self._hist[i])
        signi = self._data_hist.Clone('signi')
        for i in range(signi.GetNbinsX()+1):
            if self._data_hist.GetBinContent(i)!=0 and sum_hist.GetBinContent(i)!=0:
                value = float(self._data_hist.GetBinContent(i) - sum_hist.GetBinContent(i))
                denominator = np.sqrt(float(pow(self._data_hist.GetBinError(i),2) + pow(sum_hist.GetBinError(i),2)))
                if denominator!=0:
                    value /= denominator
                    signi.SetBinContent(i,value)
                    signi.SetBinError(i,1)
        x = []
        y = []
        err = []
        for j in range(0,len(self._error_hist)):
            x_i = []
            x_i = []
            y_i = []
            err_i = []
            for i in range(signi.GetNbinsX()+1):
                x_i.append(sum_hist.GetBinLowEdge(i))
                x_i.append(sum_hist.GetBinLowEdge(i) + sum_hist.GetBinWidth(i))
                if self._Style_cont.Get_error_bands_center() == 'ref':
                    y_i.append(0.)
                    y_i.append(0.)
                elif self._Style_cont.Get_error_bands_center() == 'val':
                    y_i.append(signi.GetBinContent(i))
                    y_i.append(signi.GetBinContent(i))
                denominator = np.sqrt(float(pow(self._data_hist.GetBinError(i),2) + pow(sum_hist.GetBinError(i),2)))
                if denominator!=0:
                    err_i.append(sum_hist.GetBinContent(i) / denominator * self._error_hist[j].GetBinContent(i))
                    err_i.append(sum_hist.GetBinContent(i) / denominator * self._error_hist[j].GetBinContent(i))
                else:
                    err_i.append(0.)
                    err_i.append(0.)
            x.append(np.array(x_i))
            y.append(np.array(y_i))
            err.append(np.array(err_i))
        return signi, x, y, err

    def _Calc_SoverSpB(self):
        sum_hist = self._hist[0].Clone('sum_hist')
        for i in range(1,len(self._hist)):
            sum_hist.Add(self._hist[i])
        soverspb = self._sig_hist[0].Clone('soverspb')
        for i in range(soverspb.GetNbinsX()+1):
            if self._sig_hist[0].GetBinContent(i)!=0 and sum_hist.GetBinContent(i)!=0:
                value = float(self._sig_hist[0].GetBinContent(i))
                denominator = np.sqrt(float(self._sig_hist[0].GetBinContent(i) + sum_hist.GetBinContent(i)))
                if denominator!=0:
                    value /= denominator
                    soverspb.SetBinContent(i,value)
                    soverspb.SetBinError(i,0)
        x = []
        y = []
        err = []
        for j in range(0,len(self._error_hist)):
            x_i = []
            y_i = []
            err_i = []
            for i in range(soverspb.GetNbinsX()+1):
                x_i.append(sum_hist.GetBinLowEdge(i))
                x_i.append(sum_hist.GetBinLowEdge(i) + sum_hist.GetBinWidth(i))
                if self._Style_cont.Get_error_bands_center() == 'ref':
                    y_i.append(0.)
                    y_i.append(0.)
                elif self._Style_cont.Get_error_bands_center() == 'val':
                    y_i.append(soverspb.GetBinContent(i))
                    y_i.append(soverspb.GetBinContent(i))
                denominator = np.sqrt(float(self._sig_hist[0].GetBinContent(i) + sum_hist.GetBinContent(i)))
                if denominator!=0:
                    err_i.append(0.)
                    err_i.append(0.)
                else:
                    err_i.append(0.)
                    err_i.append(0.)
            x.append(np.array(x_i))
            y.append(np.array(y_i))
            err.append(np.array(err_i))
        return soverspb, x, y, err

    def _show_only_some(self, x, pos):
        s = str(int(x))
        if s[0] in ('4'):
            return s
        return ''

    def _Draw_Error_Bands(self, axis1):
        sum_hist = self._hist[0].Clone('sum_hist')
        for i in range(1,len(self._hist)):
            sum_hist.Add(self._hist[i])
        x = []
        y = []
        err = []
        for j in range(0,len(self._error_hist)):
            x_i = []
            y_i = []
            err_i = []
            for i in range(sum_hist.GetNbinsX()+1):
                x_i.append(sum_hist.GetBinLowEdge(i))
                y_i.append(sum_hist.GetBinContent(i))
                x_i.append(sum_hist.GetBinLowEdge(i) + sum_hist.GetBinWidth(i))
                y_i.append(sum_hist.GetBinContent(i))
                err_i.append(sum_hist.GetBinContent(i)*abs(self._error_hist[j].GetBinContent(i)))
                err_i.append(sum_hist.GetBinContent(i)*abs(self._error_hist[j].GetBinContent(i)))
            x.append(np.array(x_i))
            y.append(np.array(y_i))
            err.append(np.array(err_i))
        self._Draw_Any_uncertainty_band(axis1, x, y, err)

    def _Draw_Any_uncertainty_band(self, axis, x, y, err):
        x_vals = x[0]
        plt.fill_between(x_vals, y[0] - np.absolute(err[0]), y[0] + np.absolute(err[0]),
                         alpha = self._Style_cont.Get_error_bands_alph(),
                         edgecolor = self._Style_cont.Get_error_bands_ecol()[0],
                         facecolor = self._Style_cont.Get_error_bands_fcol()[0],
                         lw = 0.7, axes = axis, zorder = 2.1)
        dummy_y_p = np.copy(y[0])
        dummy_y_m = np.copy(y[0])
        dummy_err_sum = np.copy(np.square(err[0]))
        if self._Style_cont.Get_error_stacking() == 'linear':
            dummy_y_p = np.add(dummy_y_p, np.absolute(err[0]))
            dummy_y_m = np.subtract(dummy_y_m, np.absolute(err[0]))
        for i in range(1,len(self._error_hist)):
            plt.fill_between(x_vals, dummy_y_p, dummy_y_p + np.absolute(err[i]),
                             alpha = self._Style_cont.Get_error_bands_alph(),
                             edgecolor = self._Style_cont.Get_error_bands_ecol()[i],
                             facecolor = self._Style_cont.Get_error_bands_fcol()[i],
                             lw = 0.7, axes = axis, zorder = 2.1)
            plt.fill_between(x_vals, dummy_y_m - np.absolute(err[i]), dummy_y_m,
                             alpha = self._Style_cont.Get_error_bands_alph(),
                             edgecolor = self._Style_cont.Get_error_bands_ecol()[i],
                             facecolor = self._Style_cont.Get_error_bands_fcol()[i],
                             lw = 0.7, axes = axis, zorder = 2.1)
            if self._Style_cont.Get_error_stacking() == 'linear':
                dummy_y_p = np.add(dummy_y_p, np.absolute(err[i]))
                dummy_y_m = np.subtract(dummy_y_m, np.absolute(err[i]))
            elif self._Style_cont.Get_error_stacking() == 'No':
                dummy_err_sum = np.add(dummy_err_sum,np.square(err[i]))
        if self._Style_cont.Get_error_stacking() == 'No':
            dummy_err_sum = np.sqrt(dummy_err_sum)
            plt.fill_between(x_vals, dummy_y_p - dummy_err_sum, dummy_y_p + dummy_err_sum,
                             alpha = 0.4,
                             edgecolor = 'black',
                             facecolor = 'grey',
                             lw = 0.7, axes = axis, zorder = 2.2)

    def _Draw_0(self):
        ## Plot a derived distribution on top of the main distribution on axis 0
        if self._add_plots[0] != '':
            self._ax0 = plt.subplot2grid((100,1), (0,0), rowspan = self._hist_start, colspan=1, sharex = self._ax1, axisbg = self._Style_cont.Get_bg_color())
            self._ax0.spines['bottom'].set_color(self._Style_cont.Get_spine_color())
            self._ax0.spines['bottom'].set_linewidth(self._Style_cont.Get_spine_line_width())
            self._ax0.spines['top'].set_color(self._Style_cont.Get_spine_color())
            self._ax0.spines['top'].set_linewidth(self._Style_cont.Get_spine_line_width())
            self._ax0.spines['left'].set_color(self._Style_cont.Get_spine_color())
            self._ax0.spines['left'].set_linewidth(self._Style_cont.Get_spine_line_width())
            self._ax0.spines['right'].set_color(self._Style_cont.Get_spine_color())
            self._ax0.spines['right'].set_linewidth(self._Style_cont.Get_spine_line_width())
            self._ax0.tick_params(axis='y', colors = self._Style_cont.Get_tick_color())
            self._ax0.tick_params(axis='x', colors = self._Style_cont.Get_tick_color())
            add_hist, x, y, err = self._Calc_additional_plot(self._add_plots[0],0)
            duke_errorbar(add_hist, xerr = False, emptybins = False, axes=self._ax0,
                          markersize = self._Style_cont.Get_marker_size(),
                          label = self._add_plots_labels[0],
                          marker = self._Style_cont.Get_marker_style(),
                          ecolor = self._Style_cont.Get_marker_color(),
                          markerfacecolor = self._Style_cont.Get_marker_color(),
                          markeredgecolor = self._Style_cont.Get_marker_color(),
                          capthick = self._Style_cont.Get_marker_error_cap_width(),
                          ignore_binns=[self._data_hist,sum(self._hist)],
                          zorder = 2.2)
            if self._add_error_bands:
                self._Draw_Any_uncertainty_band(self._ax0, x, y, err)
            self._ax0.set_ylim(ymin = add_hist.min()*1.1, ymax = add_hist.max()*1.1)
            if self._Style_cont.Get_xmin() != -1 and self._Style_cont.Get_xmax() != -1:
                self._ax0.set_xlim(xmin = self._Style_cont.Get_xmin(), xmax = self._Style_cont.Get_xmax())
            self._ax0.axhline(self._add_plots_ref_line[0], color = self._Style_cont.Get_ref_line_color())
            self._ax0.set_ylabel(self._add_plots_labels[0], color = self._Style_cont.Get_label_text_color(), va='top', ha='left', size = self._Style_cont.Get_axis_title_font_size(), weight = 'medium')
            self._ax0.yaxis.set_label_coords(self._Style_cont.Get_y_label_offset(),1.)
            self._ax0.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5, prune='lower'))
            plt.setp(self._ax0.get_xticklabels(), visible=False)
            return None
        return None

    def _Draw_main(self):
        ## Create the figure for all subplots
        self._fig = plt.figure(figsize=(6, 6), dpi=100, facecolor=self._Style_cont.Get_bg_color())
        ## Create the subplot for the main distribution
        self._ax1 = plt.subplot2grid((100,1), (self._hist_start,0), rowspan = self._hist_height, colspan = 1, axisbg = self._Style_cont.Get_bg_color())
        ## If specified in the style container set logarithmic axis
        if self._Style_cont.Get_logy():
            self._ax1.set_yscale('log')
        if self._Style_cont.Get_logx():
            self._ax1.set_xscale('log')
        if self._Style_cont.Get_grid():
            self._ax1.grid(True)
            gridlines = self._ax1.get_xgridlines()
            gridlines.extend( self._ax1.get_ygridlines() )
            for line in gridlines:
                line.set_linestyle(self._Style_cont.Get_grid_style())
                line.set_linewidth(self._Style_cont.Get_grid_width())
                line.set_color(self._Style_cont.Get_grid_color())
        ## Crete the standard plots with histograms
        if self._Style_cont.Get_kind() == 'Standard' or self._Style_cont.Get_kind() == 'Lines':
            if len(self._hist) == 0:
                if not self._data and len(self._sig_hist) == 0:
                    print('\n\tyou have to add some histogram that should be plotted,')
                    print('\tthere are no background, signal or data histograms.\n')
                    sys.exit(42)
                if self._data:
                    data_handle = rplt.errorbar(self._data_hist, xerr = False, emptybins = False, axes = self._ax1,
                                  markersize = self._Style_cont.Get_marker_size(),
                                  marker = self._Style_cont.Get_marker_style(),
                                  ecolor = self._Style_cont.Get_marker_color(),
                                  markerfacecolor = self._Style_cont.Get_marker_color(),
                                  markeredgecolor = self._Style_cont.Get_marker_color(),
                                  capthick = self._Style_cont.Get_marker_error_cap_width())
                if len(self._sig_hist) > 0:
                    rplt.hist(self._sig_hist, stacked = False, axes = self._ax1)
            else:
                hist_handle = rplt.hist(self._hist, stacked = True, axes = self._ax1, zorder = 2)
                if self._data:
                    data_handle = rplt.errorbar(self._data_hist, xerr = False, emptybins = False, axes = self._ax1,
                                  markersize = self._Style_cont.Get_marker_size(),
                                  marker = self._Style_cont.Get_marker_style(),
                                  ecolor = self._Style_cont.Get_marker_color(),
                                  markerfacecolor = self._Style_cont.Get_marker_color(),
                                  markeredgecolor = self._Style_cont.Get_marker_color(),
                                  capthick = self._Style_cont.Get_marker_error_cap_width())
                if len(self._sig_hist) > 0:
                    rplt.hist(self._sig_hist, stacked = False, axes = self._ax1)
        ## Create the main plot with graphs
        elif self._Style_cont.Get_kind() == 'Graphs':
            if len(self._hist) == 0 and not self._data and len(self._sig_hist) == 0:
                print('\n\tyou have to add some histogram that should be plotted,')
                print('\tthere are no background, signal or data histograms.\n')
                sys.exit(42)
            else:
                for item in self._hist:
                    graph_handle = rplt.errorbar(item, xerr = False, emptybins = False, axes = self._ax1,
                                   markersize = self._Style_cont.Get_marker_size(),
                                   marker = self._Style_cont.Get_marker_style(),
                                   ecolor = item.GetLineColor(),
                                   markerfacecolor = item.GetLineColor(),
                                   markeredgecolor = item.GetLineColor(),
                                   capthick = self._Style_cont.Get_marker_error_cap_width())
                for item in self._sig_hist:
                    graph_handle = rplt.errorbar(item, xerr = False, emptybins = False, axes = self._ax1,
                                   markersize = self._Style_cont.Get_marker_size(),
                                   marker = self._Style_cont.Get_marker_style(),
                                   ecolor = item.GetLineColor(),
                                   markerfacecolor = item.GetLineColor(),
                                   markeredgecolor = item.GetLineColor(),
                                   capthick = self._Style_cont.Get_marker_error_cap_width())
                if self._data:
                    data_handle = rplt.errorbar(self._data_hist, xerr = False, emptybins = False, axes = self._ax1,
                                  markersize = self._Style_cont.Get_marker_size(),
                                  marker = self._Style_cont.Get_marker_style(),
                                  ecolor = self._Style_cont.Get_marker_color(),
                                  markerfacecolor = self._Style_cont.Get_marker_color(),
                                  markeredgecolor = self._Style_cont.Get_marker_color(),
                                  capthick = self._Style_cont.Get_marker_error_cap_width())
        ## If defined draw error bands
        if self._add_error_bands:
            self._Draw_Error_Bands(self._ax1)
        ## If specified change the axis ranges
        if self._Style_cont.Get_ymin() != -1 and self._Style_cont.Get_ymax() != -1:
            self._ax1.set_ylim(ymin = self._Style_cont.Get_ymin(), ymax = self._Style_cont.Get_ymax())
        if self._Style_cont.Get_xmin() != -1 and self._Style_cont.Get_xmax() != -1:
            self._ax1.set_xlim(xmin = self._Style_cont.Get_xmin(), xmax = self._Style_cont.Get_xmax())
        ## Set the y-axis title and its options
        self._ax1.set_ylabel(self._Style_cont.Get_yaxis_title(), color=self._Style_cont.Get_label_text_color(), va='top', ha='left', size = self._Style_cont.Get_axis_title_font_size(), weight = 'medium')
        self._ax1.yaxis.set_label_coords(self._Style_cont.Get_y_label_offset(),0.9)
        ## If no other additional plots, set the x-axis title
        if not (self._add_plots[1] != '' or self._add_plots[2] != ''):
            plt.xlabel(self._Style_cont.Get_xaxis_title(), color = self._Style_cont.Get_label_text_color(), position = (1., -0.1), va = 'top', ha = 'right', size = self._Style_cont.Get_axis_title_font_size(), weight = 'medium')
        ## If defined show the minor tick marks
        if self._Style_cont.Get_show_minor_tick_labels():
            self._ax1.yaxis.set_minor_formatter(plt.FormatStrFormatter('%d'))
            self._ax1.yaxis.set_minor_formatter(plt.FuncFormatter(self._show_only_some))
        ## Set the properties of the plot spine
        self._ax1.spines['bottom'].set_color(self._Style_cont.Get_spine_color())
        self._ax1.spines['bottom'].set_linewidth(self._Style_cont.Get_spine_line_width())
        self._ax1.spines['top'].set_color(self._Style_cont.Get_spine_color())
        self._ax1.spines['top'].set_linewidth(self._Style_cont.Get_spine_line_width())
        self._ax1.spines['left'].set_color(self._Style_cont.Get_spine_color())
        self._ax1.spines['left'].set_linewidth(self._Style_cont.Get_spine_line_width())
        self._ax1.spines['right'].set_color(self._Style_cont.Get_spine_color())
        self._ax1.spines['right'].set_linewidth(self._Style_cont.Get_spine_line_width())
        ## Set the properties of the tick marks
        self._ax1.tick_params(axis = 'y', colors = self._Style_cont.Get_tick_color())
        self._ax1.tick_params(axis = 'x', colors = self._Style_cont.Get_tick_color())
        ## Add the legend
        self._Add_legend()
        return None

    def _Draw_main_axis(self):
        ## Create the figure for all subplots
        self._fig = plt.figure(figsize=(6, 6), dpi=100, facecolor=self._Style_cont.Get_bg_color())
        ## Create the subplot for the main distribution
        self._ax1 = plt.subplot2grid((100,1), (self._hist_start,0), rowspan = self._hist_height, colspan = 1, axisbg = self._Style_cont.Get_bg_color())
        par1 = self._ax1.twinx()
        ## If specified in the style container set logarithmic axis
        if self._Style_cont.Get_logy():
            self._ax1.set_yscale('log')
        if self._Style_cont.Get_histaxis_logy():
            par1.set_yscale('log')
        if self._Style_cont.Get_logx():
            self._ax1.set_xscale('log')
        if self._Style_cont.Get_grid():
            self._ax1.grid(True)
            gridlines = self._ax1.get_xgridlines()
            gridlines.extend( self._ax1.get_ygridlines() )
            for line in gridlines:
                line.set_linestyle(self._Style_cont.Get_grid_style())
                line.set_linewidth(self._Style_cont.Get_grid_width())
                line.set_color(self._Style_cont.Get_grid_color())
        ## Crete the standard plots with histograms
        if self._Style_cont.Get_kind() == 'Standard' or self._Style_cont.Get_kind() == 'Lines':
            if len(self._hist) == 0:
                if not self._data and len(self._sig_hist) == 0:
                    print('\n\tyou have to add some histogram that should be plotted,')
                    print('\tthere are no background, signal or data histograms.\n')
                    sys.exit(42)
                if self._data:
                    data_handle = rplt.errorbar(self._data_hist, xerr = False, emptybins = False, axes = self._ax1,
                                  markersize = self._Style_cont.Get_marker_size(),
                                  marker = self._Style_cont.Get_marker_style(),
                                  ecolor = self._Style_cont.Get_marker_color(),
                                  markerfacecolor = self._Style_cont.Get_marker_color(),
                                  markeredgecolor = self._Style_cont.Get_marker_color(),
                                  capthick = self._Style_cont.Get_marker_error_cap_width())
                if len(self._sig_hist) > 0:
                    rplt.hist(self._sig_hist, stacked = False, axes = self._ax1)
                rplt.hist(self._hist_axis, stacked = False, axes = par1)
            else:
                hist_handle = rplt.hist(self._hist, stacked = True, axes = self._ax1, zorder = 2)
                if self._data:
                    data_handle = rplt.errorbar(self._data_hist, xerr = False, emptybins = False, axes = self._ax1,
                                  markersize = self._Style_cont.Get_marker_size(),
                                  marker = self._Style_cont.Get_marker_style(),
                                  ecolor = self._Style_cont.Get_marker_color(),
                                  markerfacecolor = self._Style_cont.Get_marker_color(),
                                  markeredgecolor = self._Style_cont.Get_marker_color(),
                                  capthick = self._Style_cont.Get_marker_error_cap_width())
                if len(self._sig_hist) > 0:
                    rplt.hist(self._sig_hist, stacked = False, axes = self._ax1)
                rplt.hist(self._hist_axis, stacked = False, axes = par1)
        ## Create the main plot with graphs
        elif self._Style_cont.Get_kind() == 'Graphs':
            if len(self._hist) == 0 and not self._data and len(self._sig_hist) == 0:
                print('\n\tyou have to add some histogram that should be plotted,')
                print('\tthere are no background, signal or data histograms.\n')
                sys.exit(42)
            else:
                for item in self._hist:
                    graph_handle = rplt.errorbar(item, xerr = False, emptybins = False, axes = self._ax1,
                                   markersize = self._Style_cont.Get_marker_size(),
                                   marker = self._Style_cont.Get_marker_style(),
                                   ecolor = item.GetLineColor(),
                                   markerfacecolor = item.GetLineColor(),
                                   markeredgecolor = item.GetLineColor(),
                                   capthick = self._Style_cont.Get_marker_error_cap_width())
                for item in self._sig_hist:
                    graph_handle = rplt.errorbar(item, xerr = False, emptybins = False, axes = self._ax1,
                                   markersize = self._Style_cont.Get_marker_size(),
                                   marker = self._Style_cont.Get_marker_style(),
                                   ecolor = item.GetLineColor(),
                                   markerfacecolor = item.GetLineColor(),
                                   markeredgecolor = item.GetLineColor(),
                                   capthick = self._Style_cont.Get_marker_error_cap_width())
                if self._data:
                    data_handle = rplt.errorbar(self._data_hist, xerr = False, emptybins = False, axes = self._ax1,
                                  markersize = self._Style_cont.Get_marker_size(),
                                  marker = self._Style_cont.Get_marker_style(),
                                  ecolor = self._Style_cont.Get_marker_color(),
                                  markerfacecolor = self._Style_cont.Get_marker_color(),
                                  markeredgecolor = self._Style_cont.Get_marker_color(),
                                  capthick = self._Style_cont.Get_marker_error_cap_width())
                for item in self._hist_axis:
                    axishist_handle = rplt.errorbar(item, xerr = False, emptybins = False, axes = par1,
                                   markersize = self._Style_cont.Get_marker_size(),
                                   marker = self._Style_cont.Get_marker_style(),
                                   ecolor = item.GetLineColor(),
                                   markerfacecolor = item.GetLineColor(),
                                   markeredgecolor = item.GetLineColor(),
                                   capthick = self._Style_cont.Get_marker_error_cap_width())
        ## If defined draw error bands
        if self._add_error_bands:
            self._Draw_Error_Bands(self._ax1)
        ## If specified change the axis ranges
        if self._Style_cont.Get_ymin() != -1 and self._Style_cont.Get_ymax() != -1:
            self._ax1.set_ylim(ymin = self._Style_cont.Get_ymin(), ymax = self._Style_cont.Get_ymax())
        if self._Style_cont.Get_histaxis_ymin() != -1 and self._Style_cont.Get_histaxis_ymax() != -1:
            par1.set_ylim(ymin = self._Style_cont.Get_histaxis_ymin(), ymax = self._Style_cont.Get_histaxis_ymax())
        if self._Style_cont.Get_xmin() != -1 and self._Style_cont.Get_xmax() != -1:
            self._ax1.set_xlim(xmin = self._Style_cont.Get_xmin(), xmax = self._Style_cont.Get_xmax())
        ## Set the y-axis title and its options
        self._ax1.set_ylabel(self._Style_cont.Get_yaxis_title(), color=self._Style_cont.Get_label_text_color(), va='top', ha='left', size = self._Style_cont.Get_axis_title_font_size(), weight = 'medium')
        self._ax1.yaxis.set_label_coords(self._Style_cont.Get_y_label_offset(),0.9)
        par1.set_ylabel(self._Style_cont.Get_histaxis_yaxis_title(), color=self._Style_cont.Get_histaxis_label_text_color(), va='top', ha='left', size = self._Style_cont.Get_axis_title_font_size(), weight = 'medium')
        par1.yaxis.set_label_coords(self._Style_cont.Get_histaxis_y_label_offset(),0.9)
        ## If no other additional plots, set the x-axis title
        if not (self._add_plots[1] != '' or self._add_plots[2] != ''):
            plt.xlabel(self._Style_cont.Get_xaxis_title(), color = self._Style_cont.Get_label_text_color(), position = (1., -0.1), va = 'top', ha = 'right', size = self._Style_cont.Get_axis_title_font_size(), weight = 'medium')
        ## If defined show the minor tick marks
        if self._Style_cont.Get_show_minor_tick_labels():
            self._ax1.yaxis.set_minor_formatter(plt.FormatStrFormatter('%d'))
            self._ax1.yaxis.set_minor_formatter(plt.FuncFormatter(self._show_only_some))
        ## Set the properties of the plot spine
        self._ax1.spines['bottom'].set_color(self._Style_cont.Get_spine_color())
        self._ax1.spines['bottom'].set_linewidth(self._Style_cont.Get_spine_line_width())
        self._ax1.spines['top'].set_color(self._Style_cont.Get_spine_color())
        self._ax1.spines['top'].set_linewidth(self._Style_cont.Get_spine_line_width())
        self._ax1.spines['left'].set_color(self._Style_cont.Get_spine_color())
        self._ax1.spines['left'].set_linewidth(self._Style_cont.Get_spine_line_width())
        self._ax1.spines['right'].set_color(self._Style_cont.Get_spine_color())
        self._ax1.spines['right'].set_linewidth(self._Style_cont.Get_spine_line_width())
        ## Set the properties of the tick marks
        self._ax1.tick_params(axis = 'y', colors = self._Style_cont.Get_tick_color())
        par1.tick_params(axis = 'y', colors = self._Style_cont.Get_histaxis_label_text_color())
        self._ax1.tick_params(axis = 'x', colors = self._Style_cont.Get_tick_color())
        ## Add the legend
        self._Add_legend()
        return None


    def _Draw_2(self):
        ## Plot a derived distribution below the main distribution on axis 2
        if self._add_plots[1] != '':
            self._ax2 = plt.subplot2grid((100,1), (self._hist_start + self._hist_height,0), rowspan = self._add_plots_height[1], colspan = 1, sharex = self._ax1, axisbg = self._Style_cont.Get_bg_color())
            add_hist, x, y, err = self._Calc_additional_plot(self._add_plots[1],1)
            duke_errorbar(add_hist, xerr = False, emptybins = False, axes = self._ax2,
                          markersize = self._Style_cont.Get_marker_size(),
                          label = self._add_plots_labels[1],
                          marker = self._Style_cont.Get_marker_style(),
                          ecolor = self._Style_cont.Get_marker_color(),
                          markerfacecolor = self._Style_cont.Get_marker_color(),
                          markeredgecolor = self._Style_cont.Get_marker_color(),
                          capthick = self._Style_cont.Get_marker_error_cap_width(),
                          ignore_binns=[self._data_hist,sum(self._hist)],
                          zorder = 2.2)
            if self._add_error_bands:
                self._Draw_Any_uncertainty_band(self._ax2, x, y, err)
            self._ax2.set_ylim(ymin = add_hist.min()*1.1, ymax = add_hist.max()*1.1)
            if self._Style_cont.Get_xmin() != -1 and self._Style_cont.Get_xmax() != -1:
                self._ax2.set_xlim(xmin = self._Style_cont.Get_xmin(), xmax = self._Style_cont.Get_xmax())
            self._ax2.axhline(self._add_plots_ref_line[1], color = self._Style_cont.Get_ref_line_color())
            self._ax2.set_ylabel(self._add_plots_labels[1], color = self._Style_cont.Get_label_text_color(), va='top', ha='left', size = self._Style_cont.Get_axis_title_font_size(), weight = 'medium')
            self._ax2.yaxis.set_label_coords(self._Style_cont.Get_y_label_offset(),1.)
            self._ax2.spines['bottom'].set_color(self._Style_cont.Get_spine_color())
            self._ax2.spines['bottom'].set_linewidth(self._Style_cont.Get_spine_line_width())
            self._ax2.spines['top'].set_color(self._Style_cont.Get_spine_color())
            self._ax2.spines['top'].set_linewidth(self._Style_cont.Get_spine_line_width())
            self._ax2.spines['left'].set_color(self._Style_cont.Get_spine_color())
            self._ax2.spines['left'].set_linewidth(self._Style_cont.Get_spine_line_width())
            self._ax2.spines['right'].set_color(self._Style_cont.Get_spine_color())
            self._ax2.spines['right'].set_linewidth(self._Style_cont.Get_spine_line_width())
            self._ax2.tick_params(axis = 'y', colors = self._Style_cont.Get_tick_color())
            self._ax2.tick_params(axis = 'x', colors = self._Style_cont.Get_tick_color())
            if self._add_plots[2] != '':
                plt.setp(self._ax2.get_xticklabels(), visible = False)
                self._ax2.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5, prune='both'))
                plt.xlabel(self._Style_cont.Get_xaxis_title(), color = self._Style_cont.Get_label_text_color(), position = (1., -0.1), va = 'top', ha = 'right', size = self._Style_cont.Get_axis_title_font_size(), weight = 'medium')
            else:
                #self._ax2.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5, prune='upper'))
                self._ax2.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5, prune='both'))
                #self._ax2.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5, prune='lower'))
                plt.xlabel(self._Style_cont.Get_xaxis_title(), color=self._Style_cont.Get_label_text_color(), position = (1., -0.1), va = 'top', ha = 'right', size = self._Style_cont.Get_axis_title_font_size(), weight = 'medium')
            plt.setp(self._ax1.get_xticklabels(), visible = False)
            return None
        return None

    def _Draw_3(self):
        ## Plot a derived distribution at the very bottom of the main distribution on axis 3
        if self._add_plots[2] != '':
            self._ax3 = plt.subplot2grid((100,1), (100 - self._add_plots_height[2],0), rowspan = self._add_plots_height[2], colspan = 1, sharex = self._ax1, axisbg = self._Style_cont.Get_bg_color())
            add_hist, x, y, err = self._Calc_additional_plot(self._add_plots[2],2)
            duke_errorbar(add_hist, xerr = False, emptybins = False, axes = self._ax3,
                          markersize = self._Style_cont.Get_marker_size(),
                          label = self._add_plots_labels[2],
                          marker = self._Style_cont.Get_marker_style(),
                          ecolor = self._Style_cont.Get_marker_color(),
                          markerfacecolor = self._Style_cont.Get_marker_color(),
                          markeredgecolor = self._Style_cont.Get_marker_color(),
                          capthick = self._Style_cont.Get_marker_error_cap_width(),
                          ignore_binns=[self._data_hist,sum(self._hist)],
                          zorder = 2.2)
            if self._add_error_bands:
                self._Draw_Any_uncertainty_band(self._ax3, x, y, err)
            self._ax3.set_ylim(ymin = add_hist.min()*1.1, ymax = add_hist.max()*1.1)
            if self._Style_cont.Get_xmin() != -1 and self._Style_cont.Get_xmax() != -1:
                self._ax3.set_xlim(xmin = self._Style_cont.Get_xmin(), xmax = self._Style_cont.Get_xmax())
            self._ax3.axhline(self._add_plots_ref_line[2], color = self._Style_cont.Get_ref_line_color())
            self._ax3.set_ylabel(self._add_plots_labels[2], color = self._Style_cont.Get_label_text_color(), va = 'top', ha = 'left', size = self._Style_cont.Get_axis_title_font_size(), weight = 'medium')
            self._ax3.yaxis.set_label_coords(self._Style_cont.Get_y_label_offset(),1.)
            #self._ax3.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5, prune='upper'))
            self._ax3.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5, prune='both'))
            self._ax3.spines['bottom'].set_color(self._Style_cont.Get_spine_color())
            self._ax3.spines['bottom'].set_linewidth(self._Style_cont.Get_spine_line_width())
            self._ax3.spines['top'].set_color(self._Style_cont.Get_spine_color())
            self._ax3.spines['top'].set_linewidth(self._Style_cont.Get_spine_line_width())
            self._ax3.spines['left'].set_color(self._Style_cont.Get_spine_color())
            self._ax3.spines['left'].set_linewidth(self._Style_cont.Get_spine_line_width())
            self._ax3.spines['right'].set_color(self._Style_cont.Get_spine_color())
            self._ax3.spines['right'].set_linewidth(self._Style_cont.Get_spine_line_width())
            self._ax3.tick_params(axis = 'y', colors = self._Style_cont.Get_tick_color())
            self._ax3.tick_params(axis = 'x', colors = self._Style_cont.Get_tick_color())
            plt.setp(self._ax1.get_xticklabels(), visible = False)
            plt.xlabel(self._Style_cont.Get_xaxis_title(), color = self._Style_cont.Get_label_text_color(), position = (1., -0.1), va = 'top', ha = 'right', size = self._Style_cont.Get_axis_title_font_size(), weight = 'medium')
            return None
        return None

    def _Draw(self):
        if self._useRoot:
            self.DrawRoot()
            return

        if len(self._hist_axis) > 0:
            self._Draw_main_axis()
        else:
            self._Draw_main()

        self._Draw_0()

        self._Draw_2()

        self._Draw_3()

        if len(self._hist_axis) > 0:
            plt.subplots_adjust(left = .10, bottom = .08, right =  .91, top = .95, wspace = .2, hspace = .0)
        else:
            plt.subplots_adjust(left = .10, bottom = .08, right =  .95, top = .95, wspace = .2, hspace = .0)
        self._Write_additional_text()

    def _SavePlot(self, out_name):
        if self._useRoot:
            self._fig.SaveAs(out_name)
            return
        if out_name[-3:] == 'pdf':
            try:
                plt.savefig(out_name, facecolor = self._fig.get_facecolor())
            except(AssertionError):
                plt.savefig(out_name[:-3] + 'eps', facecolor = self._fig.get_facecolor())
                command="epstopdf %s"%(out_name[:-3] + 'eps')
                command=command.split(" ")
                subprocess.call(command)
                command="rm %s"%(out_name[:-3] + 'eps')
                command=command.split(" ")
                subprocess.call(command)
        else:
            plt.savefig(out_name, facecolor = self._fig.get_facecolor())

    def _AddRootLegend(self):
        if self._Style_cont.Get_no_legend():
            return
        if self._Style_cont.Get_LegendPosition() == self._Style_cont.Get_cmsTextPosition():
            self._Style_cont.Get_LegendPosition().addYspace(self._Style_cont.Get_cmsTextPosition().getY()-self._Style_cont.Get_LegendPosition().getY()-0.02)

        numberOfEntries=len(self._hist)+len(self._sig_hist)
        if self._data:
            numberOfEntries+=1
        textSize=0.7*self._referenceHeight
        self.leg = Legend(numberOfEntries,rightmargin=1.-self._Style_cont.Get_LegendPosition().getX(),topmargin=1.-self._Style_cont.Get_LegendPosition().getY(),textfont=self._Style_cont.additionalTextFont,textsize=self._Style_cont.legendTextSize,entryheight=textSize,entrysep=textSize*0.1)
        self.leg.SetFillStyle(0)
        self.leg.SetBorderSize(0)
        #self.leg.SetFillColor(ROOT.kWhite)
        for h in self._hist:
            self.leg.AddEntry(h,h.GetTitle(), "f")
        for h in self._sig_hist:
            self.leg.AddEntry(h,h.GetTitle(), "l")
        if self._data:
            self.leg.AddEntry(self._data_hist,"data","ep")

    def _AddPlotBelow(self):
        ## setup the window and pads to draw a ratio

        nAdditionalPlots=0
        for i in self._add_plots:
            if i !='':
                nAdditionalPlots+=1

        self._canvas.cd()

        #make the canvas bigger
        addedPadheight=sum(self._add_plots_height)
        expansion_factor=1.+addedPadheight*0.01
        self._canvas.height=int(self._canvas.height*expansion_factor)

        #self._referenceHeight*=(1.+addedPadheight)

        ##set all the margins correct
        # this can be imporved perhaps
        if nAdditionalPlots!=0:
            self._mainPad = Pad(0,addedPadheight,1,1)
            self._mainPad.SetBottomMargin(0.0)
        else:
            self._mainPad = Pad(0,0,1,1)

        if nAdditionalPlots==1:
            self.additionalPad.append(Pad(0,0,1,self._add_plots_height[0]))
            self.additionalPad[0].SetBottomMargin(0.12/self._add_plots_height[0]*(1.-self._add_plots_height[0]))
        elif nAdditionalPlots==2:
            self.additionalPad.append(Pad(0,self._add_plots_height[1],1,addedPadheight))
            self.additionalPad.append(Pad(0,0,1,self._add_plots_height[1]))
            self.additionalPad[0].SetBottomMargin(0.)
            self.additionalPad[1].SetBottomMargin(0.12/self._add_plots_height[1]*(1.-self._add_plots_height[1]))

        #we dont want to have a scale here
        for pad in self.additionalPad:
            pad.SetTopMargin(0.0)
            pad.SetLeftMargin(0.1)
        self._mainPad.SetLeftMargin(0.1)

        self._mainPad.Draw()
        for pad in self.additionalPad:
            pad.Draw()
        ffactor=1
        for i in range(len(self._add_plots)):
            if self._add_plots[i] !='':
                add_hist, x, y, err = self._Calc_additional_plot(self._add_plots[i],i)
                self.additionalPad[i].cd()
                add_hist.Draw()
                add_hist.GetXaxis().SetTitle(self._Style_cont._xaxis_title.replace("$\\mathsf{","").replace("}$",""))
                add_hist.GetYaxis().SetTitle(self._add_plots_labels[i])


                add_hist.GetXaxis().SetTitleFont(self._Style_cont.additionalTextFont)
                add_hist.GetYaxis().SetTitleFont(self._Style_cont.additionalTextFont)
                add_hist.GetYaxis().SetTitleSize(self._Style_cont.axisTitleTextSize)
                add_hist.GetXaxis().SetTitleSize(self._Style_cont.axisTitleTextSize)

                add_hist.GetXaxis().SetLabelFont(self._Style_cont.additionalTextFont)
                add_hist.GetYaxis().SetLabelFont(self._Style_cont.additionalTextFont)
                add_hist.GetXaxis().SetLabelSize(self._Style_cont.axisLabelTextSize)
                add_hist.GetYaxis().SetLabelSize(self._Style_cont.axisLabelTextSize)


                add_hist.GetYaxis().SetTitleOffset(self._Style_cont.axisOffsetY)
                add_hist.GetYaxis().SetNdivisions(205)
                add_hist.GetXaxis().SetTitleOffset(self._Style_cont.axisOffsetX/addedPadheight/1.5)
        self._mainPad.cd()


    def DrawRoot(self):
        import rootplotlib as rooLib

        rooLib.init()

        self._canvas= Canvas()
        self._referenceHeight=0.05/self._canvas.GetAbsHNDC() *self._canvas.GetWw() / self._canvas.GetWh()
        #fix axis if not root style:
        for h in self._allHists:
            if h is None:
                continue
            if "\\" in h.GetXaxis().GetTitle():
                h.GetXaxis().SetTitle( h.GetXaxis().GetTitle().replace("$\\mathsf{","").replace("\\","#").replace("}$","")  )
            else:
                break


        self._AddPlotBelow()

        self._Draw_main_root()

        rnd=rounding(sigdigits=3)
        lumitext=""
        if self._Style_cont.Get_lumi_val() > 1000:
            lumitext='%s fb^{-1} (%.0f TeV)'%(rnd.latex(self._Style_cont.Get_lumi_val()/1000.),self._Style_cont.Get_cms_val())
        else:
            lumitext='%.1f pb^{-1} (%.0f TeV)'%(self._Style_cont.Get_lumi_val(),self._Style_cont.Get_cms_val())
        deco=rooLib.CmsDecoration(sc_obj=self._Style_cont ,extraText=self._Style_cont.Get_additional_text(), additionalText=None, lumiText=lumitext, position=self._Style_cont.Get_cmsTextPosition(), pad=ROOT.gPad)
        deco.Draw()
        self._canvas.Update()
        self._fig=self._canvas

    def _Draw_main_root(self):
        drawnObjects=[]
        same=""
        if len(self._hist)>0:
            self._hist[0].Draw("AXIS")
            drawnObjects.append(self._hist[0])
            same=" same"
            hs=HistStack(hists=self._hist)
            hs.Draw("hist"+same)
            drawnObjects.append(hs)
        for sg_hist in self._sig_hist:
            sg_hist.Draw("hist"+same)
            drawnObjects.append(sg_hist)
            same=" same"
        if self._data:
            self._data_hist.Draw("E"+same)
            drawnObjects.append(self._data_hist)
            same=" same"
        self._AddRootLegend()
        self._mainPad.SetLogy(self._Style_cont.Get_logy())
        self._mainPad.SetLogx(self._Style_cont.Get_logx())


        if self._Style_cont.Get_ymin() != -1 and self._Style_cont.Get_ymax() != -1:
            drawnObjects[0].GetYaxis().SetRangeUser(self._Style_cont.Get_ymin(),self._Style_cont.Get_ymax())
        else:
            maximum=None
            minimum=None
            for h in drawnObjects:
                if maximum is None:
                    maximum=h.max()
                    minimum=h.min()
                else:
                    if maximum<h.max():
                        maximum=h.max()
                    if minimum>h.min():
                        minimum=h.min()
            if minimum<=0:
                minimum=1
            drawnObjects[0].GetYaxis().SetRangeUser(minimum*0.02,maximum*100.)
        if self._Style_cont.Get_xmin() != -1 and self._Style_cont.Get_xmax() != -1:
            drawnObjects[0].GetXaxis().SetRangeUser(self._Style_cont.Get_xmin(),self._Style_cont.Get_xmax())

        drawnObjects[0].GetXaxis().SetLabelFont(self._Style_cont.additionalTextFont)
        drawnObjects[0].GetYaxis().SetLabelFont(self._Style_cont.additionalTextFont)
        drawnObjects[0].GetXaxis().SetLabelSize(self._Style_cont.axisLabelTextSize)
        drawnObjects[0].GetYaxis().SetLabelSize(self._Style_cont.axisLabelTextSize)

        drawnObjects[0].GetXaxis().SetTitle(drawnObjects[0].GetXaxis().GetTitle().replace("$\\mathsf{","").replace("}$",""))
        drawnObjects[0].GetXaxis().SetTitleFont(self._Style_cont.additionalTextFont)
        drawnObjects[0].GetYaxis().SetTitleFont(self._Style_cont.additionalTextFont)
        drawnObjects[0].GetYaxis().SetTitleSize(self._Style_cont.axisTitleTextSize)
        drawnObjects[0].GetXaxis().SetTitleSize(self._Style_cont.axisTitleTextSize)
        drawnObjects[0].GetYaxis().SetTitleOffset(self._Style_cont.axisOffsetY)
        drawnObjects[0].GetXaxis().SetTitleOffset(self._Style_cont.axisOffsetX)


        ROOT.gPad.RedrawAxis("g")
        self.leg.Draw()

        #for easy debug
        #self._canvas.Update()
        #raw_input()

    def _checker(self):
        pass
        #try:
            #print('histo with name:' + self._hist[0].GetName())
        #except AttributeError:
            #print('No histogram added')
        #print('with height: ' + str(self._hist_height) + ' and start: ' + str(self._hist_start))
