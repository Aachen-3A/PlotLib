#!/bin/env python

import ROOT
import numpy as np
from rootpy.plotting import Hist, HistStack, Legend, Canvas, Graph, Pad
import matplotlib
from matplotlib import rc
import matplotlib.ticker as mticker
import rootpy.plotting.root2matplotlib as rplt
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from plotlib import duke_errorbar
from operator import methodcaller
from rounding import rounding



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
    # @param[in] style String of which style should be used (default = 'Plain')
    # @param[in] hist List of background histograms (default = [])
    # @param[in] sig List of signal histograms (default = [])
    # @param[in] data_hist Data histogram (default = None)
    # @param[in] data Bool if data should be plotted (default = False)
    # @param[in] cms double to specify displayed center of mass energy (default = 13 (TeV))
    # @param[in] lumi double to specify displayed luminosity value (default = 42000 (pb-1))
    # @param[in] data Bool if data should be plotted (default = False)
    # @param[in] kwargs dict of key word arguments that will be passed to style
    def __init__(self, style = 'Plain', hist = [], sig = [], data_hist = None, data = False, useRoot=False, cms = 13, lumi = 42000, **kwargs):
        ## style variables
        self._style                = style
        ## BG histograms
        self._hist                 = hist
        self._hist_height          = 100
        self._hist_start           = 0
        ## SG histograms
        self._sig_hist             = sig
        ## Data histograms
        self._data                 = data
        self._data_hist            = data_hist
        ## Additional plots
        self._add_plots            = ['', '', '']
        self._add_plots_height     = [0, 0, 0]
        self._add_plots_labels     = ['', '', '']
        self._add_plots_ref_line   = [0, 0, 0]
        self._annotations_modified = False
        self._add_error_bands      = False
        self._error_hist           = []
        self._fig                  = None
        self._useRoot              = useRoot
        self._cms_val              = cms
        self._lumi_val             = lumi
        self._allHists=self._hist+self._sig_hist+[self._data_hist]
        self._Set_style(**kwargs)

    ## del function
    #
    # This deletes the main objects nedded to not get a crash at the end!
    def __del__(self):
        plt.close()
        del self._hist
        del self._sig_hist
        del self._data_hist
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
    # At the moment 'Ratio', 'Diff', 'Signi' and 'DiffRatio' are available
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
            self._Set_style(cmsPositon=self.cmsTextPosition.getText(), legendPosition= self.LegendPosition.getText())
        else:
            print('for pos %.0f is already %s planned, so that is not possible'%(pos,self.add_plots[pos]))

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
            self._error_bands_labl = labels
        self._error_bands_center = band_center
        self._error_stacking = stacking

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
    def Set_axis(self, logx = False, logy = True, ymin = -1, ymax = -1, xmin = -1, xmax = -1):
        self._logx = logx
        self._logy = logy
        self._ymin = ymin
        self._ymax = ymax
        self._xmin = xmin
        self._xmax = xmax

    ## Function to save the complete plot
    #
    # This function saves the plot you which is stored in the object so create it first
    # @param[in] out_name name of the output file
    def SavePlot(self, out_name):
        self._SavePlot(out_name)

    def ChangeStyle(self,**kwargs):

        for key in kwargs:
            if hasattr(self,"_"+key):
                setattr(self,"_"+key,kwargs[key])
            else:
                print "No attribute _%s in plotter"%key


    ##------------------------------------------------------------------
    ## Private functions
    ##------------------------------------------------------------------
    def _Set_style(self,cmsPositon="upper right",legendPosition="upper right"):
        if self._useRoot:
            self.cmsTextFont          = 61   # Fonts
            self.lumiTextFont         = 42
            self.extraTextFont        = 52
            self.additionalTextFont   = 42
            self.cmsTextSize          = 0.9  #Text sizes
            self.lumiTextSize         = 0.6
            self.extraTextSize        = 0.76*self.cmsTextSize
            self.additionalTextSize   = 1.0*self.extraTextSize
            self.legendTextSize       = self.extraTextSize*0.8
            self.lumiTextOffset       = 0.2
            self.extraTextOffset      = 2.5  # only used in outOfFrame version
            self.axisTextSize         = 0.9
            self.axisOffset           = 1.3
            self._ratio_pad           ={}
            self.rootMemory           =[]
        else:
            matplotlib.rcParams.update({'font.size': 10})
            matplotlib.rcParams.update({'lines.linewidth' : 1})
        #rc('text', usetex=True)
        self._xaxis_title      = self._allHists[0].xaxis.GetTitle()
        self._yaxis_title      = self._allHists[0].yaxis.GetTitle()
        self._additional_text  = 'Preliminary'
        self._y_label_offset   = -0.11
        self._error_bands_ecol = ['darkmagenta','darkcyan']
        self._error_bands_fcol = ['m','cyan']
        self._error_bands_alph = 0.7
        self._error_bands_labl = ['Sys. uncert. 1','Sys. uncert. 2']
        self._error_bands_center = 'ref'
        self._error_stacking = 'No'
        self._spine_line_width = 0.5
        self._logx = False
        self._logy = True
        self._ymin = -1
        self._ymax = -1
        self._xmin = -1
        self._xmax = -1
        self.cmsTextPosition=position(cmsPositon,isText=True)
        self.LegendPosition=position(legendPosition)
        if self._style == 'CMS':
            self._add_cms_text           = True
            self._add_lumi_text          = True
            self._label_text_color       = 'black'
            self._annotation_text_color  = 'black'
            self._bg_color               = 'w'
            self._ref_line_color         = 'blue'
            self._spine_color            = 'black'
            self._tick_color             = 'black'
            self._marker_style           = 'o'
            self._marker_size            = 3
            self._marker_color           = 'black'
            self._marker_error_cap_width = 0
            self._cms_text_alignment     = 'row'
            self._show_minor_tick_labels = False
            self._legend_font_size       = 9
            if self._add_plots[0] != '':
                self.cmsTextPosition.addYspace(  -0.9 * self._add_plots_height[0] / 100.)
            if self._add_plots[1] != '':
                self.cmsTextPosition.addYspace(  0.9 * self._add_plots_height[1] / 100.)
            if self._add_plots[2] != '':
                self.cmsTextPosition.addYspace(  0.9 * self._add_plots_height[2] / 100.)
        elif self._style == 'Plain':
            self._add_cms_text           = False
            self._add_lumi_text          = False
            self._label_text_color       = 'black'
            self._annotation_text_color  = 'black'
            self._bg_color               = 'w'
            self._ref_line_color         = 'blue'
            self._spine_color            = 'black'
            self._tick_color             = 'black'
            self._marker_style           = 'o'
            self._marker_size            = 4
            self._marker_color           = 'black'
            self._marker_error_cap_width = 1
            self._cms_text_alignment     = 'row'
            self._show_minor_tick_labels = True
            self._legend_font_size       = 10
            if self._add_plots[0] != '':
                self.cmsTextPosition.addYspace(  -0.8 * self._add_plots_height[0] / 100.)
            if self._add_plots[1] != '':
                self.cmsTextPosition.addYspace(  0.8 * self._add_plots_height[1] / 100.)
            if self._add_plots[2] != '':
                self.cmsTextPosition.addYspace(  0.8 * self._add_plots_height[2] / 100.)
        elif self._style == 'Cool':
            self._add_cms_text           = True
            self._add_lumi_text          = True
            self._label_text_color       = 'white'
            self._annotation_text_color  = 'white'
            self._bg_color               = '#07000d'
            self._ref_line_color         = 'y'
            self._spine_color            = '#5998ff'
            self._tick_color             = 'w'
            self._marker_style           = 'o'
            self._marker_size            = 3
            self._marker_color           = 'lightgray'
            self._marker_error_cap_width = 0
            self._cms_text_alignment     = 'column'
            self._show_minor_tick_labels = False
            self._legend_font_size       = 9

    def _Write_additional_text(self):
        if self._add_lumi_text:
            self._lumi_val=float(self._lumi_val)
            if self._lumi_val > 1000:
                self._fig.text(0.945, 0.955, '$%.1f\,\mathrm{fb^{-1}} (%.0f\,\mathrm{TeV})$'%(self._lumi_val/1000,self._cms_val), va='bottom', ha='right', color=self._annotation_text_color, size=12)
            else:
                self._fig.text(0.945, 0.955, '$%.0f\,\mathrm{pb^{-1}} (%.0f\,\mathrm{TeV})$'%(self._lumi_val,self._cms_val), va='bottom', ha='right', color=self._annotation_text_color, size=12)
        if self._add_cms_text:
            if self._cms_text_alignment == 'row':
                self._fig.text(self.cmsTextPosition.getX(), self.cmsTextPosition.getY(), 'CMS', va='bottom', ha='left', color=self._annotation_text_color, size=14, weight='bold')
                self._fig.text(self.cmsTextPosition.getX(), self.cmsTextPosition.getY()-0.03, self._additional_text, va='bottom', ha='left', color=self._annotation_text_color, size=10, style = 'italic')
            elif self._cms_text_alignment == 'column':
                self._fig.text(self.cmsTextPosition.getX(), self.cmsTextPosition.getY(), 'CMS', va='bottom', ha='left', color=self._annotation_text_color, size=14, weight='bold')
                self._fig.text(self.cmsTextPosition.getX() + 0.08, self.cmsTextPosition.getY(), self._additional_text, va='bottom', ha='left', color=self._annotation_text_color, size=10, style = 'italic')
            else:
                print('At the moment only ''row'' and ''column'' are allowed alignment values')

    def _Add_legend(self):
        if self._add_plots[0] != '':
            self.LegendPosition.addYspace(-(0.85 * self._add_plots_height[0] / 100.))
        if self._add_plots[1] != '':
            self.LegendPosition.addYspace(  0.8 * self._add_plots_height[1] / 100.)
        if self._add_plots[2] != '':
            self.LegendPosition.addYspace(  0.8 * self._add_plots_height[2] / 100.)

        if self.LegendPosition==self.cmsTextPosition:
            self.LegendPosition.addYspace(self.cmsTextPosition.getY()-self.LegendPosition.getY()-0.02)
        handle_list = []
        label_list = []
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
                col_patch = mpatches.Patch(facecolor = self._error_bands_fcol[i], edgecolor = self._error_bands_ecol[i] , alpha = self._error_bands_alph , lw = 0.7)
                handle_list.append(col_patch)
                label_list.append(self._error_bands_labl[i])
            if self._error_stacking == 'No':
                col_patch = mpatches.Patch(facecolor = 'grey', edgecolor = 'black' , alpha = 0.4 , lw = 0.7)
                handle_list.append(col_patch)
                label_list.append('syst. sum')
        if self._data:
            dat_line=plt.errorbar([], [],xerr = False,yerr=True, markersize = self._marker_size,
                              marker = self._marker_style,
                              color = self._marker_color,
                              capthick = self._marker_error_cap_width)
            handle_list.append(dat_line)
            label_list.append(self._data_hist.GetTitle())

        self.leg = plt.legend(handle_list, label_list,
                    loc = 'upper right',
                    bbox_to_anchor=(self.LegendPosition.getX(),self.LegendPosition.getY()),
                    bbox_transform=plt.gcf().transFigure,
                    numpoints = 1,
                    frameon = False,
                    fontsize = self._legend_font_size)

        for text in self.leg.get_texts():
            text.set_color(self._annotation_text_color)

    def _Compiler(self):
        if len(self._hist) == 0:
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
        if self._logy and self._data:
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
                if self._error_bands_center == 'ref':
                    y_i.append(1.)
                    y_i.append(1.)
                elif self._error_bands_center == 'val':
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
                if self._error_bands_center == 'ref':
                    y_i.append(0.)
                    y_i.append(0.)
                elif self._error_bands_center == 'val':
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
                if self._error_bands_center == 'ref':
                    y_i.append(0.)
                    y_i.append(0.)
                elif self._error_bands_center == 'val':
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
            y_i = []
            err_i = []
            for i in range(signi.GetNbinsX()+1):
                x_i.append(sum_hist.GetBinLowEdge(i))
                x_i.append(sum_hist.GetBinLowEdge(i) + sum_hist.GetBinWidth(i))
                if self._error_bands_center == 'ref':
                    y_i.append(0.)
                    y_i.append(0.)
                elif self._error_bands_center == 'val':
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
                         alpha = self._error_bands_alph,
                         edgecolor = self._error_bands_ecol[0],
                         facecolor = self._error_bands_fcol[0],
                         lw = 0.7, axes = axis, zorder = 2.1)
        dummy_y_p = np.copy(y[0])
        dummy_y_m = np.copy(y[0])
        dummy_err_sum = np.copy(np.square(err[0]))
        if self._error_stacking == 'linear':
            dummy_y_p = np.add(dummy_y_p, np.absolute(err[0]))
            dummy_y_m = np.subtract(dummy_y_m, np.absolute(err[0]))
        for i in range(1,len(self._error_hist)):
            plt.fill_between(x_vals, dummy_y_p, dummy_y_p + np.absolute(err[i]),
                             alpha = self._error_bands_alph,
                             edgecolor = self._error_bands_ecol[i],
                             facecolor = self._error_bands_fcol[i],
                             lw = 0.7, axes = axis, zorder = 2.1)
            plt.fill_between(x_vals, dummy_y_m - np.absolute(err[i]), dummy_y_m,
                             alpha = self._error_bands_alph,
                             edgecolor = self._error_bands_ecol[i],
                             facecolor = self._error_bands_fcol[i],
                             lw = 0.7, axes = axis, zorder = 2.1)
            if self._error_stacking == 'linear':
                dummy_y_p = np.add(dummy_y_p, np.absolute(err[i]))
                dummy_y_m = np.subtract(dummy_y_m, np.absolute(err[i]))
            elif self._error_stacking == 'No':
                dummy_err_sum = np.add(dummy_err_sum,np.square(err[i]))
        if self._error_stacking == 'No':
            dummy_err_sum = np.sqrt(dummy_err_sum)
            plt.fill_between(x_vals, dummy_y_p - dummy_err_sum, dummy_y_p + dummy_err_sum,
                             alpha = 0.4,
                             edgecolor = 'black',
                             facecolor = 'grey',
                             lw = 0.7, axes = axis, zorder = 2.2)

    def _Draw_0(self, axis1):
        ## Plot a derived distribution on top of the main distribution on axis 0
        if self._add_plots[0] != '':
            ax0 = plt.subplot2grid((100,1), (0,0), rowspan = self._hist_start, colspan=1, sharex = axis1, axisbg = self._bg_color)
            ax0.spines['bottom'].set_color(self._spine_color)
            ax0.spines['bottom'].set_linewidth(self._spine_line_width)
            ax0.spines['top'].set_color(self._spine_color)
            ax0.spines['top'].set_linewidth(self._spine_line_width)
            ax0.spines['left'].set_color(self._spine_color)
            ax0.spines['left'].set_linewidth(self._spine_line_width)
            ax0.spines['right'].set_color(self._spine_color)
            ax0.spines['right'].set_linewidth(self._spine_line_width)
            ax0.tick_params(axis='y', colors = self._tick_color)
            ax0.tick_params(axis='x', colors = self._tick_color)
            add_hist, x, y, err = self._Calc_additional_plot(self._add_plots[0],0)
            duke_errorbar(add_hist, xerr = False, emptybins = False, axes=ax0,
                          markersize = self._marker_size,
                          label = self._add_plots_labels[0],
                          marker = self._marker_style,
                          ecolor = self._marker_color,
                          markerfacecolor = self._marker_color,
                          markeredgecolor = self._marker_color,
                          capthick = self._marker_error_cap_width,
                          ignore_binns=[self._data_hist,sum(self._hist)],
                          zorder = 2.2)
            if self._add_error_bands:
                self._Draw_Any_uncertainty_band(ax0, x, y, err)
            ax0.set_ylim(ymin = add_hist.min()*1.1, ymax = add_hist.max()*1.1)
            if self._xmin != -1 and self._xmax != -1:
                ax0.set_xlim(xmin = self._xmin, xmax = self._xmax)
            ax0.axhline(self._add_plots_ref_line[0], color = self._ref_line_color)
            ax0.set_ylabel(self._add_plots_labels[0], color = self._label_text_color, va='top', ha='left')
            ax0.yaxis.set_label_coords(self._y_label_offset,1.)
            ax0.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5, prune='lower'))
            plt.setp(ax0.get_xticklabels(), visible=False)
            return ax0
        return None

    def _Draw_main(self):
        self._fig = plt.figure(figsize=(6, 6), dpi=100, facecolor=self._bg_color)
        ## Plot the main distribution on axis 1
        ax1 = plt.subplot2grid((100,1), (self._hist_start,0), rowspan = self._hist_height, colspan = 1, axisbg = self._bg_color)
        if self._logy:
            ax1.set_yscale('log')
        if self._logx:
            ax1.set_xscale('log')
        if len(self._hist) == 0:
            if not self._data and len(self._sig_hist) == 0:
                print('you have to add some histogram that should be plotted, there are no background, signal or data histograms.')
            if self._data:
                data_handle = rplt.errorbar(self._data_hist, xerr = False, emptybins = False, axes = ax1,
                              markersize = self._marker_size,
                              marker = self._marker_style,
                              ecolor = self._marker_color,
                              markerfacecolor = self._marker_color,
                              markeredgecolor = self._marker_color,
                              capthick = self._marker_error_cap_width)
            if len(self._sig_hist) > 0:
                rplt.hist(self._sig_hist, stacked = False, axes = ax1)
        elif len(self._hist) == 1:
            hist_handle = rplt.hist(self._hist[0], stacked = False, axes = ax1, zorder = 2)
            if self._data:
                data_handle = rplt.errorbar(self._data_hist, xerr = False, emptybins = False, axes = ax1,
                              markersize = self._marker_size,
                              marker = self._marker_style,
                              ecolor = self._marker_color,
                              markerfacecolor = self._marker_color,
                              markeredgecolor = self._marker_color,
                              capthick = self._marker_error_cap_width)
            if len(self._sig_hist) > 0:
                rplt.hist(self._sig_hist, stacked = False, axes = ax1)
        else:
            hist_handle = rplt.hist(self._hist, stacked = True, axes = ax1, zorder = 2)
            if self._data:
                data_handle = rplt.errorbar(self._data_hist, xerr = False, emptybins = False, axes = ax1,
                              markersize = self._marker_size,
                              marker = self._marker_style,
                              ecolor = self._marker_color,
                              markerfacecolor = self._marker_color,
                              markeredgecolor = self._marker_color,
                              capthick = self._marker_error_cap_width)
            if len(self._sig_hist) > 0:
                rplt.hist(self._sig_hist, stacked = False, axes = ax1)
        if self._add_error_bands:
            self._Draw_Error_Bands(ax1)
        if self._ymin != -1 and self._ymax != -1:
            ax1.set_ylim(ymin = self._ymin, ymax = self._ymax)
        if self._xmin != -1 and self._xmax != -1:
            ax1.set_xlim(xmin = self._xmin, xmax = self._xmax)
        ax1.set_ylabel(self._yaxis_title, color=self._label_text_color, va='top', ha='left')
        ax1.yaxis.set_label_coords(self._y_label_offset,0.9)
        if not (self._add_plots[1] != '' or self._add_plots[2] != ''):
            plt.xlabel(self._xaxis_title, color = self._label_text_color, position = (1., -0.1), va = 'top', ha = 'right')
        if self._show_minor_tick_labels:
            ax1.yaxis.set_minor_formatter(plt.FormatStrFormatter('%d'))
            ax1.yaxis.set_minor_formatter(plt.FuncFormatter(self._show_only_some))
        #ax1.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5, prune='lower'))
        ax1.spines['bottom'].set_color(self._spine_color)
        ax1.spines['bottom'].set_linewidth(self._spine_line_width)
        ax1.spines['top'].set_color(self._spine_color)
        ax1.spines['top'].set_linewidth(self._spine_line_width)
        ax1.spines['left'].set_color(self._spine_color)
        ax1.spines['left'].set_linewidth(self._spine_line_width)
        ax1.spines['right'].set_color(self._spine_color)
        ax1.spines['right'].set_linewidth(self._spine_line_width)
        ax1.tick_params(axis = 'y', colors = self._tick_color)
        ax1.tick_params(axis = 'x', colors = self._tick_color)
        self._Add_legend()
        return ax1

    def _Draw_2(self, axis1):
        ## Plot a derived distribution below the main distribution on axis 2
        if self._add_plots[1] != '':
            ax2 = plt.subplot2grid((100,1), (self._hist_start + self._hist_height,0), rowspan = self._add_plots_height[1], colspan = 1, sharex = axis1, axisbg = self._bg_color)
            add_hist, x, y, err = self._Calc_additional_plot(self._add_plots[1],1)
            duke_errorbar(add_hist, xerr = False, emptybins = False, axes = ax2,
                          markersize = self._marker_size,
                          label = self._add_plots_labels[1],
                          marker = self._marker_style,
                          ecolor = self._marker_color,
                          markerfacecolor = self._marker_color,
                          markeredgecolor = self._marker_color,
                          capthick = self._marker_error_cap_width,
                          ignore_binns=[self._data_hist,sum(self._hist)],
                          zorder = 2.2)
            if self._add_error_bands:
                self._Draw_Any_uncertainty_band(ax2, x, y, err)
            ax2.set_ylim(ymin = add_hist.min()*1.1, ymax = add_hist.max()*1.1)
            if self._xmin != -1 and self._xmax != -1:
                ax2.set_xlim(xmin = self._xmin, xmax = self._xmax)
            ax2.axhline(self._add_plots_ref_line[1], color = self._ref_line_color)
            ax2.set_ylabel(self._add_plots_labels[1], color = self._label_text_color, va='top', ha='left')
            ax2.yaxis.set_label_coords(self._y_label_offset,1.)
            ax2.spines['bottom'].set_color(self._spine_color)
            ax2.spines['bottom'].set_linewidth(self._spine_line_width)
            ax2.spines['top'].set_color(self._spine_color)
            ax2.spines['top'].set_linewidth(self._spine_line_width)
            ax2.spines['left'].set_color(self._spine_color)
            ax2.spines['left'].set_linewidth(self._spine_line_width)
            ax2.spines['right'].set_color(self._spine_color)
            ax2.spines['right'].set_linewidth(self._spine_line_width)
            ax2.tick_params(axis = 'y', colors = self._tick_color)
            ax2.tick_params(axis = 'x', colors = self._tick_color)
            if self._add_plots[2] != '':
                plt.setp(ax2.get_xticklabels(), visible = False)
                ax2.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5, prune='both'))
                plt.xlabel(self._xaxis_title, color = self._label_text_color, position = (1., -0.1), va = 'top', ha = 'right')
            else:
                #ax2.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5, prune='upper'))
                ax2.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5, prune='both'))
                #ax2.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5, prune='lower'))
                plt.xlabel(self._xaxis_title, color=self._label_text_color, position = (1., -0.1), va = 'top', ha = 'right')
            plt.setp(axis1.get_xticklabels(), visible = False)
            return ax2
        return None

    def _Draw_3(self, axis1):
        ## Plot a derived distribution at the very bottom of the main distribution on axis 3
        if self._add_plots[2] != '':
            ax3 = plt.subplot2grid((100,1), (100 - self._add_plots_height[2],0), rowspan = self._add_plots_height[2], colspan = 1, sharex = axis1, axisbg = self._bg_color)
            add_hist, x, y, err = self._Calc_additional_plot(self._add_plots[2],2)
            duke_errorbar(add_hist, xerr = False, emptybins = False, axes = ax3,
                          markersize = self._marker_size,
                          label = self._add_plots_labels[2],
                          marker = self._marker_style,
                          ecolor = self._marker_color,
                          markerfacecolor = self._marker_color,
                          markeredgecolor = self._marker_color,
                          capthick = self._marker_error_cap_width,
                          ignore_binns=[self._data_hist,sum(self._hist)],
                          zorder = 2.2)
            if self._add_error_bands:
                self._Draw_Any_uncertainty_band(ax3, x, y, err)
            ax3.set_ylim(ymin = add_hist.min()*1.1, ymax = add_hist.max()*1.1)
            if self._xmin != -1 and self._xmax != -1:
                ax3.set_xlim(xmin = self._xmin, xmax = self._xmax)
            ax3.axhline(self._add_plots_ref_line[2], color = self._ref_line_color)
            ax3.set_ylabel(self._add_plots_labels[2], color = self._label_text_color, va = 'top', ha = 'left')
            ax3.yaxis.set_label_coords(self._y_label_offset,1.)
            #ax3.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5, prune='upper'))
            ax3.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5, prune='both'))
            ax3.spines['bottom'].set_color(self._spine_color)
            ax3.spines['bottom'].set_linewidth(self._spine_line_width)
            ax3.spines['top'].set_color(self._spine_color)
            ax3.spines['top'].set_linewidth(self._spine_line_width)
            ax3.spines['left'].set_color(self._spine_color)
            ax3.spines['left'].set_linewidth(self._spine_line_width)
            ax3.spines['right'].set_color(self._spine_color)
            ax3.spines['right'].set_linewidth(self._spine_line_width)
            ax3.tick_params(axis = 'y', colors = self._tick_color)
            ax3.tick_params(axis = 'x', colors = self._tick_color)
            plt.setp(axis1.get_xticklabels(), visible = False)
            plt.xlabel(self._xaxis_title, color = self._label_text_color, position = (1., -0.1), va = 'top', ha = 'right')
            return ax3
        return None

    def _Draw(self):
        if self._useRoot:
            self.DrawRoot()
            return
        ax1 = self._Draw_main()

        ax0 = self._Draw_0(ax1)

        ax2 = self._Draw_2(ax1)

        ax3 = self._Draw_3(ax1)

        plt.subplots_adjust(left = .10, bottom = .08, right =  .95, top = .95, wspace = .2, hspace = .0)
        self._Write_additional_text()

    def _SavePlot(self, out_name):
        if self._useRoot:
            self._fig.SaveAs(out_name)
            return
        plt.savefig(out_name, facecolor = self._fig.get_facecolor())

    def _AddRootLegend(self):
        if self.LegendPosition==self.cmsTextPosition:
            self.LegendPosition.addYspace(self.cmsTextPosition.getY()-self.LegendPosition.getY()-0.02)

        numberOfEntries=len(self._hist)+len(self._sig_hist)
        if self._data:
            numberOfEntries+=1
        textSize=self.legendTextSize*self._referenceHeight
        self.leg = Legend(numberOfEntries,rightmargin=1.-self.LegendPosition.getX(),topmargin=1.-self.LegendPosition.getY(),textfont=42,textsize=textSize,entryheight=textSize,entrysep=textSize*0.1)
        self.leg.SetFillStyle(0)
        self.leg.SetBorderSize(0)
        self.leg.SetFillColor(ROOT.kWhite)
        for h in self._hist:
            self.leg.AddEntry(h,h.GetTitle(), "f")
        for h in self._sig_hist:
            self.leg.AddEntry(h,h.GetTitle(), "l")
        if self._data:
            self.leg.AddEntry(self._data_hist,"data","ep")


    def _AddPlotBelow(self, pos=2):
        ## setup the window and pads to draw a ratio
        if self._add_plots[pos] != '':
            self._canvas.cd()

            expansion_factor=1.+self._add_plots_height[pos]*0.01
            ## expand canvas
            #self._canvas.SetWindowSize(self._canvas.width,self._canvas.height*expansion_factor)
            #height(self._canvas.height()*expansion_factor )
            self._canvas.height=int(self._canvas.height*expansion_factor)

            # resize drawing pad
            # base length - ( base length / expansion factor )
            y_ndc = 0.97 - (0.97 / expansion_factor)
            ROOT.gPad.SetPad(0.01, y_ndc, 0.98, 0.98)
            #update_pad()

            # draw new pad for ratio on canvas
            self._canvas.canvas.cd()
            # base length - ( drawing pad bottom margin * base length / expansion factor )
            y_ndc = 0.97 - ((1 - ROOT.gPad.GetBottomMargin()) * 0.97 / expansion_factor)
            self._ratio_pad[pos] = Pad(0.01, 0.01, 0.98, y_ndc)

            # adjust settings, due to different scale
            self._ratio_pad[pos].SetTopMargin(0.0)
            self._ratio_pad[pos].SetBottomMargin(0.30)
            self._ratio_pad[pos].Draw()

            self._ratio_pad[pos].cd()

            add_hist, x, y, err = self._Calc_additional_plot(self._add_plots[pos],pos)
            self.rootMemory.append(add_hist)
            add_hist.Draw()



    #def update_pad(self):
        #"""Updates the pad and redraws the axis"""

        #if not ROOT.gPad:
            #print "No active pad."
            #return

        #ROOT.gPad.Modified()
        #ROOT.gPad.Update()
        #ROOT.gPad.RedrawAxis()


    def DrawRoot(self):
        import rootplotlib as rooLib

        rooLib.init()

        self._canvas= Canvas()
        self._referenceHeight=0.05/self._canvas.GetAbsHNDC() *self._canvas.GetWw() / self._canvas.GetWh()
        self._AddPlotBelow()

        self._Draw_main_root()

        rnd=rounding(sigdigits=3)
        lumitext=""
        if self._lumi_val > 1000:
            lumitext='%s fb^{-1} (%.0f TeV)'%(rnd.latex(self._lumi_val/1000.),self._cms_val)
        else:
            lumitext='%.1f pb^{-1} (%.0f TeV)'%(self._lumi_val,self._cms_val)
        deco=rooLib.CmsDecoration(extraText=self._additional_text, additionalText=None, lumiText=lumitext, align="left", valign="top", pad=ROOT.gPad)
        deco.Draw()
        self._canvas.Update()
        self._fig=self._canvas

    def _Draw_main_root(self):
        print self._canvas.find_all_primitives()
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
        if self._logy:
            self._canvas.SetLogy(True)
        if self._ymin != -1 and self._ymax != -1:
            drawnObjects[0].GetYaxis().SetRangeUser(self._ymin,self._ymax)
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
        if self._xmin != -1 and self._xmax != -1:
            drawnObjects[0].GetXaxis().SetRangeUser(self._xmin,self._xmax)
        drawnObjects[0].GetXaxis().SetTitle(drawnObjects[0].GetXaxis().GetTitle().replace("$\\mathsf{","").replace("}$",""))
        print drawnObjects[0].GetXaxis().GetTitle()
        drawnObjects[0].GetYaxis().SetTitleSize(self.axisTextSize*self._referenceHeight)
        drawnObjects[0].GetXaxis().SetTitleSize(self.axisTextSize*self._referenceHeight)
        drawnObjects[0].GetYaxis().SetTitleOffset(self.axisOffset)



        ROOT.gPad.RedrawAxis("g")
        self.leg.Draw()


    def _checker(self):
        pass
        #try:
            #print('histo with name:' + self._hist[0].GetName())
        #except AttributeError:
            #print('No histogram added')
        #print('with height: ' + str(self._hist_height) + ' and start: ' + str(self._hist_start))

class position():
    def __init__(self,positiontext="upper right", refference="", isText=False):

        self._positiontext=positiontext
        if not isinstance(positiontext,str):
            self._definedCoorinates=True
            self._valign="left"
            self._align="left"
        else:
            self._definedCoorinates=False
            self._valign=self._positiontext.split(" ")[0]
            self._align=self._positiontext.split(" ")[1]
        self.addY=0
        self.addX=0
        self._isText=isText
        self._correctcms={"left":0.,
                    "middle":0.,
                    "right":-0.15,
                    "upper":-0.04,
                    "center":0.,
                    "lower":0.,
        }
        if self._definedCoorinates:
            self._x=self._positiontext[0]
            self._y=self._positiontext[1]

    def __eq__(self,other):
        return (self._positiontext==other._positiontext)

    def addYspace(self,y):
        if self._valign=="upper" and y<0.:
            self.addY+=y
        elif self._valign!="upper" and self.getY()+y>0.1:
            self.addY+=y

    def addXspace(self,x):
        self.addX+=x



    def setPosition(self,positiontext):
        self._positiontext=positiontext
        self.valign=self._positiontext.split(" ")[0]
        self.align=self._positiontext.split(" ")[1]

    def getText(self):
        return self._positiontext

    def getX(self):
        if self._definedCoorinates:
            return self._x
        alignDict={
                    "left":0.12,
                    "middle":0.5,
                    "right":0.95,
        }
        if self._isText:
            return self.addX+alignDict[self._align]+self._correctcms[self._align]
        return self.addX+alignDict[self._align]

    def getY(self):
        if self._definedCoorinates:
            return self._y
        alignDict={
                    "upper":0.95,
                    "center":0.5,
                    "lower":0.12,
        }
        if self._isText:
            return self.addY+alignDict[self._valign]+self._correctcms[self._valign]
        return self.addY+alignDict[self._valign]












