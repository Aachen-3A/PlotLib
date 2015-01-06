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
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

##@class plotter
# Class to collect matplotlib functions
#
# To use the various matplotlib functions to produce the standard
# plots (with data[optional], signal[optional], background and
# uncertainties[optional]). Also different analysis distributions
# like ratio or siginficance can be added.
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
    def __init__(self, style = 'Plain', hist = [], sig = [], data_hist = None, data = False):
        ## style variables
        self.style                = style
        ## BG histograms
        self.hist                 = hist
        self.hist_height          = 100
        self.hist_start           = 0
        ## SG histograms
        self.sig_hist             = sig
        ## Data histograms
        self.data                 = data
        self.data_hist            = data_hist
        ## Additional plots
        self.add_plots            = ['', '', '']
        self.add_plots_height     = [0, 0, 0]
        self.add_plots_labels     = ['', '', '']
        self.add_plots_ref_line   = [0, 0, 0]
        self.annotations_modified = False
        self.add_error_bands      = False
        self._Set_style()

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
        self.data = doData
        self.data_hist = data_hist

    ## Function to decide if you want to draw the data
    #
    # This function is used to set the bool to plot the data
    # @param[in] doData Boolean if the data should be drawn (default = True)
    def Draw_data(self, doData = True):
        self.data = doData

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
        if self.add_plots[pos] == '':
            self.add_plots[pos] = plot
            self.add_plots_height[pos] = height
            if label == '':
                self.add_plots_labels[pos] = plot
            else:
                self.add_plots_labels[pos] = label
            self._Set_style()
        else:
            print('for pos %.0f is already %s planned, so that is not possible'%(pos,self.add_plots[pos]))

    ## Function to add a histogram to the background list
    #
    # This function is used to add an additional histogram to the list of
    # background histogram.
    # @param[in] histo Histogram that should be added
    def Add_histo(self, histo):
        self.hist.append(histo)

    ## Function to set the uncertainty histogram
    #
    # This function is used to add the systematic uncertainty histogram.
    # @param[in] histo Histogram that contains as bin content the relativ systematic uncertainty.
    def Add_error_hist(self, histo):
        self.add_error_bands = True
        self.error_hist = histo

    ##------------------------------------------------------------------
    ## Private functions
    ##------------------------------------------------------------------
    def _Set_style(self):
        matplotlib.rcParams.update({'font.size': 10})
        #rc('text', usetex=True)
        self.xaxis_title      = self.hist[0].xaxis.GetTitle()
        self.yaxis_title      = self.hist[0].yaxis.GetTitle()
        self.lumi_val         = 42000
        self.cms_val          = 13
        self.additional_text  = '$Preliminary$'
        self.y_label_offset   = -0.11
        self.error_bands_ecol = 'gray'
        self.error_bands_fcol = 'lightgray'
        self.error_bands_alph = 0.7
        self.error_bands_labl = 'Sys. uncert.'
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
            if self.add_plots[0] == '':
                self.cms_text_x         = 0.8
                self.cms_text_y         = 0.9
            else:
                self.cms_text_x         = 0.8
                self.cms_text_y         = 0.9 - (0.8 * self.add_plots_height[0] / 100.)
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
            if self.add_plots[0] == '':
                self.cms_text_x         = 0.8
                self.cms_text_y         = 0.9
            else:
                self.cms_text_x         = 0.8
                self.cms_text_y         = 0.9 - (0.8 * self.add_plots_height[0] / 100.)
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

    def _Add_legend(self, axis):
        self.legend_x = 0.95
        if self.style == 'Cool':
            if self.add_plots[0] == '':
                self.legend_y = 0.9
            else:
                self.legend_y = 0.9 - (0.8 * self.add_plots_height[0] / 100.)
        else:
            self.legend_y = self.cms_text_y - 0.02
        handle_list = []
        label_list = []
        for item in self.hist:
            col_patch = mpatches.Patch(color=item.GetFillColor())
            handle_list.append(col_patch)
            label_list.append(item.GetTitle())
        for item in self.sig_hist:
            col_patch = mlines.Line2D([], [], color=item.GetLineColor(), markersize=0)
            handle_list.append(col_patch)
            label_list.append(item.GetTitle())
        if self.add_error_bands:
            col_patch = mpatches.Patch(facecolor=self.error_bands_fcol,edgecolor=self.error_bands_ecol,alpha=self.error_bands_alph,lw=0.7)
            handle_list.append(col_patch)
            label_list.append(self.error_bands_labl)
        if self.data:
            dat_line = mlines.Line2D([], [], color=self.marker_color, marker=self.marker_style, markersize=self.marker_size)
            handle_list.append(dat_line)
            label_list.append(self.data_hist.GetTitle())
        leg = axis.legend(handle_list, label_list,
                    loc = 'upper right',
                    bbox_to_anchor=(self.legend_x, self.legend_y),
                    bbox_transform=plt.gcf().transFigure,
                    numpoints = 1,
                    frameon = False)
        for text in leg.get_texts():
            text.set_color(self.annotation_text_color)

    def _Compiler(self):
        if self.add_plots[0] != '':
            self.hist_start = self.add_plots_height[0]
            self.hist_height -= self.add_plots_height[0]
        if self.add_plots[1] != '':
            self.hist_height -= self.add_plots_height[1]
        if self.add_plots[2] != '':
            self.hist_height -= self.add_plots_height[2]

    def _Calc_additional_plot(self, plot, pos):
        if plot == 'Ratio':
            self.add_plots_labels[pos] = 'Data/MC'
            self.add_plots_ref_line[pos] = 1.
            return self._Calc_ratio()
        elif plot == 'Diff':
            self.add_plots_labels[pos] = 'Data - MC'
            self.add_plots_ref_line[pos] = 0.
            return self._Calc_diff()
        elif plot == 'Signi':
            self.add_plots_labels[pos] = 'Significance'
            self.add_plots_ref_line[pos] = 0.
            return self._Calc_signi()
        elif plot == 'DiffRatio':
            self.add_plots_labels[pos] = '(Data - MC)/MC'
            self.add_plots_ref_line[pos] = 0.
            return self._Calc_diffratio()
        else:
            print('%s is not implemented yet as an additional plot, feel free to include this functionallity')

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

    def _Calc_diffratio(self):
        diff = self.data_hist.Clone('diffratio')
        for item in self.hist:
            diff.Add(item,-1)
        sum_hist = self.hist[0].Clone('sum_hist')
        for i in range(1,len(self.hist)):
            sum_hist.Add(self.hist[i])
        diff.Divide(sum_hist)
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

    def _Draw_Error_Bands(self, axis1):
        sum_hist = self.hist[0].Clone('sum_hist')
        for i in range(1,len(self.hist)):
            sum_hist.Add(self.hist[i])
        x = []
        y = []
        error = []
        for i in range(sum_hist.GetNbinsX()+1):
            x.append(sum_hist.GetBinLowEdge(i))
            y.append(sum_hist.GetBinContent(i))
            x.append(sum_hist.GetBinLowEdge(i) + sum_hist.GetBinWidth(i))
            y.append(sum_hist.GetBinContent(i))
            error.append(sum_hist.GetBinContent(i)*abs(self.error_hist.GetBinContent(i)))
            error.append(sum_hist.GetBinContent(i)*abs(self.error_hist.GetBinContent(i)))
        plt.fill_between(np.array(x), np.array(y)-np.array(error), np.array(y)+np.array(error), alpha = self.error_bands_alph, edgecolor=self.error_bands_ecol, facecolor=self.error_bands_fcol, lw = 0.7, axes = axis1, zorder=2.1)

    def _Draw_0(self, axis1):
        ## Plot a derived distribution on top of the main distribution on axis 0
        if self.add_plots[0] != '':
            ax0 = plt.subplot2grid((100,1), (0,0), rowspan=self.hist_start, colspan=1, sharex = axis1, axisbg = self.bg_color)
            ax0.spines['bottom'].set_color(self.spine_color)
            ax0.spines['top'].set_color(self.spine_color)
            ax0.spines['left'].set_color(self.spine_color)
            ax0.spines['right'].set_color(self.spine_color)
            ax0.tick_params(axis='y', colors=self.tick_color)
            ax0.tick_params(axis='x', colors=self.tick_color)
            add_hist = self._Calc_additional_plot(self.add_plots[0],0)
            rplt.errorbar(add_hist, xerr=False, emptybins=False, axes=ax0,
                          markersize=self.marker_size,
                          label=self.add_plots_labels[0],
                          marker = self.marker_style,
                          ecolor = self.marker_color,
                          markerfacecolor = self.marker_color,
                          markeredgecolor = self.marker_color,
                          capthick = self.marker_error_cap_width)
            ax0.axhline(self.add_plots_ref_line[0], color=self.ref_line_color)
            ax0.set_ylabel(self.add_plots_labels[0], color=self.label_text_color, va='top', ha='left')
            ax0.yaxis.set_label_coords(self.y_label_offset,1.)
            ax0.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5, prune='lower'))
            plt.setp(ax0.get_xticklabels(), visible=False)
            return ax0
        return None

    def _Draw_main(self):
        self.fig = plt.figure(figsize=(6, 6), dpi=100, facecolor=self.bg_color)
        ## Plot the main distribution on axis 1
        ax1 = plt.subplot2grid((100,1), (self.hist_start,0), rowspan=self.hist_height, colspan=1, axisbg = self.bg_color)
        if len(self.hist) == 1:
            hist_handle = rplt.hist(self.hist[0], stacked=False, axes=ax1, zorder=2)
            if self.data:
                data_handle = rplt.errorbar(self.data_hist, xerr=False, emptybins=False, axes=ax1, 
                              markersize=self.marker_size,
                              marker = self.marker_style,
                              ecolor = self.marker_color,
                              markerfacecolor = self.marker_color,
                              markeredgecolor = self.marker_color,
                              capthick = self.marker_error_cap_width)
        else:
            hist_handle = rplt.hist(self.hist, stacked=True, axes=ax1, zorder=2)
            if self.data:
                data_handle = rplt.errorbar(self.data_hist, xerr=False, emptybins=False, axes=ax1,
                              markersize=self.marker_size,
                              marker = self.marker_style,
                              ecolor = self.marker_color,
                              markerfacecolor = self.marker_color,
                              markeredgecolor = self.marker_color,
                              capthick = self.marker_error_cap_width)
        if self.add_error_bands:
            self._Draw_Error_Bands(ax1)
        if len(self.sig_hist) > 0:
            rplt.hist(self.sig_hist, stacked = False, axes = ax1)
        ax1.set_ylabel(self.yaxis_title, color=self.label_text_color, va='top', ha='left')
        ax1.yaxis.set_label_coords(self.y_label_offset,0.9)
        self._Add_legend(ax1)
        if not (self.add_plots[1] != '' or self.add_plots[2] != ''):
            plt.xlabel(self.xaxis_title, color=self.label_text_color, position=(1., -0.1), va='top', ha='right')
        ax1.yaxis.set_major_locator(mticker.MaxNLocator(prune='lower'))
        ax1.spines['bottom'].set_color(self.spine_color)
        ax1.spines['top'].set_color(self.spine_color)
        ax1.spines['left'].set_color(self.spine_color)
        ax1.spines['right'].set_color(self.spine_color)
        ax1.tick_params(axis='y', colors=self.tick_color)
        ax1.tick_params(axis='x', colors=self.tick_color)
        return ax1

    def _Draw_2(self, axis1):
        ## Plot a derived distribution below the main distribution on axis 2
        if self.add_plots[1] != '':
            ax2 = plt.subplot2grid((100,1), (self.hist_start+self.hist_height,0), rowspan=self.add_plots_height[1], colspan=1, sharex = axis1, axisbg = self.bg_color)
            add_hist = self._Calc_additional_plot(self.add_plots[1],1)
            rplt.errorbar(add_hist, xerr=False, emptybins=False, axes=ax2,
                          markersize=self.marker_size,
                          label=self.add_plots_labels[1],
                          marker = self.marker_style,
                          ecolor = self.marker_color,
                          markerfacecolor = self.marker_color,
                          markeredgecolor = self.marker_color,
                          capthick = self.marker_error_cap_width)
            ax2.axhline(self.add_plots_ref_line[1], color=self.ref_line_color)
            ax2.set_ylabel(self.add_plots_labels[1], color=self.label_text_color, va='top', ha='left')
            ax2.yaxis.set_label_coords(self.y_label_offset,1.)
            ax2.spines['bottom'].set_color(self.spine_color)
            ax2.spines['top'].set_color(self.spine_color)
            ax2.spines['left'].set_color(self.spine_color)
            ax2.spines['right'].set_color(self.spine_color)
            ax2.tick_params(axis='y', colors=self.tick_color)
            ax2.tick_params(axis='x', colors=self.tick_color)
            if self.add_plots[2] != '':
                plt.setp(ax2.get_xticklabels(), visible=False)
                ax2.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5, prune='both'))
                plt.xlabel(self.xaxis_title, color=self.label_text_color, position=(1., -0.1), va='top', ha='right')
            else:
                ax2.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5, prune='upper'))
                plt.xlabel(self.xaxis_title, color=self.label_text_color, position=(1., -0.1), va='top', ha='right')
            plt.setp(axis1.get_xticklabels(), visible=False)
            return ax2
        return None

    def _Draw_3(self, axis1):
        ## Plot a derived distribution at the very bottom of the main distribution on axis 3
        if self.add_plots[2] != '':
            ax3 = plt.subplot2grid((100,1), (100-self.add_plots_height[2],0), rowspan=self.add_plots_height[2], colspan=1, sharex = axis1, axisbg = self.bg_color)
            add_hist = self._Calc_additional_plot(self.add_plots[2],2)
            rplt.errorbar(add_hist, xerr=False, emptybins=False, axes=ax3,
                          markersize=self.marker_size,
                          label=self.add_plots_labels[2],
                          marker = self.marker_style,
                          ecolor = self.marker_color,
                          markerfacecolor = self.marker_color,
                          markeredgecolor = self.marker_color,
                          capthick = self.marker_error_cap_width)
            ax3.axhline(self.add_plots_ref_line[2], color=self.ref_line_color)
            ax3.set_ylabel(self.add_plots_labels[2], color=self.label_text_color, va='top', ha='left')
            ax3.yaxis.set_label_coords(self.y_label_offset,1.)
            ax3.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5, prune='upper'))
            ax3.spines['bottom'].set_color(self.spine_color)
            ax3.spines['top'].set_color(self.spine_color)
            ax3.spines['left'].set_color(self.spine_color)
            ax3.spines['right'].set_color(self.spine_color)
            ax3.tick_params(axis='y', colors=self.tick_color)
            ax3.tick_params(axis='x', colors=self.tick_color)
            plt.setp(axis1.get_xticklabels(), visible=False)
            plt.xlabel(self.xaxis_title, color=self.label_text_color, position=(1., -0.1), va='top', ha='right')
            return ax3
        return None

    def _Draw(self):

        ax1 = self._Draw_main()

        ax0 = self._Draw_0(ax1)

        ax2 = self._Draw_2(ax1)

        ax3 = self._Draw_3(ax1)

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
