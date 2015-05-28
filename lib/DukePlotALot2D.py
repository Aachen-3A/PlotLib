#!/usr/bin/python

# contourplot.py
# version 1.0
#
# Creates 2D plot from xlist, ylist and zlist with optional contours.
#
# Viktor Kutzner (kutzner@physik.rwth-aachen.de)

from __future__ import division
from numpy import linspace, meshgrid
from matplotlib.mlab import griddata
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter

import rootpy.plotting.root2matplotlib as rplt

from rootpy.plotting import Hist2D, Hist

from DukePlotALot import *

import style_class as sc

class plotter2D():
    def __init__(self, hist = None, style = sc.style_container()):
        self._hist = hist

        self._x_projection_size = -1
        self._y_projection_size = -1
        self._x_starting_point = 0
        self._y_starting_point = 0

        self._Style_cont = style
        self._Style_cont.InitStyle()

    def __del__(self):
        pass

    def make_plot(self, out_name, individual = False):
        self._hist = self._Get_rootpy_hist2d()
        self._Draw_main()
        self._SavePlot(out_name)
        if individual:
            test = plotter2D(hist = self._hist, style = self._Style_cont)
            test.Set_axis(ymin = self._Style_cont.Get_ymin(), ymax = self._Style_cont.Get_ymax())
            test.make_plot(out_name.replace('.','_2D.'))
            if self._x_projection_size > 0:
                test = plotter(hist = [self._Get_rootpy_hist1d(self._x_projection)],style=self._Style_cont)
                test.Set_axis(ymin = self._Style_cont.Get_zmin(), ymax = self._Style_cont.Get_zmax(), xmin = self._Style_cont.Get_xmin(), xmax = self._Style_cont.Get_xmax())
                test.make_plot(out_name.replace('.','_x.'))
            if self._y_projection_size > 0:
                test = plotter(hist = [self._Get_rootpy_hist1d(self._y_projection)],style=self._Style_cont)
                test.Set_axis(ymin = self._Style_cont.Get_zmin(), ymax = self._Style_cont.Get_zmax(), xmin = self._Style_cont.Get_ymin(), xmax = self._Style_cont.Get_ymax())
                test.make_plot(out_name.replace('.','_y.'))

    def create_plot(self):
        self._hist = self._Get_rootpy_hist2d()
        self._Draw_main()

    def show_fig(self):
        self._fig.show()
        raw_input('bla')

    def save_plot(self, out_name):
        self._SavePlot(out_name)

    def Get_2D_axis(self):
        return self._ax1

    def Get_x_projection_axis(self):
        return self._ax2

    def Get_y_projection_axis(self):
        return self._ax3

    def Get_z_axis(self):
        return self._z_axis

    def Get_x_projection_hist(self):
        ret_hist = self._x_projection
        return ret_hist

    def Get_y_projection_hist(self):
        return self._y_projection

    def Add_x_projection(self, size = 15):
        self._x_projection_size = size
        self._x_starting_point = size

    def Add_y_projection(self, size = 15):
        self._y_projection_size = size
        self._y_starting_point = size

    def Set_axis(self, logx = False, logy = False, ymin = -1, ymax = -1, xmin = -1, xmax = -1, zmin = -1, zmax = -1, grid = False):
        self._Style_cont.Set_axis(logx = logx, logy = logy, ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax, grid = grid)

    def _Write_additional_text(self):
        if self._Style_cont.Get_add_lumi_text():
            self._Style_cont.Set_lumi_val(float(self._Style_cont.Get_lumi_val()))
            if self._Style_cont.Get_lumi_val() == 0:
                self._fig.text(0.945, 0.955, '$%.0f\,\mathrm{TeV}$'%(self._Style_cont.Get_cms_val()), va='bottom', ha='right', color=self._Style_cont.Get_annotation_text_color(), size=12)
            elif self._Style_cont.Get_lumi_val() >= 1000:
                self._fig.text(0.945, 0.955, '$%.1f\,\mathrm{fb^{-1}} (%.0f\,\mathrm{TeV})$'%(self._Style_cont.Get_lumi_val()/1000,self._Style_cont.Get_cms_val()), va='bottom', ha='right', color=self._Style_cont.Get_annotation_text_color(), size=12)
            else:
                self._fig.text(0.945, 0.955, '$%.0f\,\mathrm{pb^{-1}} (%.0f\,\mathrm{TeV})$'%(self._Style_cont.Get_lumi_val(),self._Style_cont.Get_cms_val()), va='bottom', ha='right', color=self._Style_cont.Get_annotation_text_color(), size=12)
        if self._Style_cont.Get_add_cms_text():
            self._fig.text(self._Style_cont.Get_cmsTextPosition().getX(), self._Style_cont.Get_cmsTextPosition().getY(), 'CMS', va='bottom', ha='left', color=self._Style_cont.Get_annotation_text_color(), size=14, weight='bold')
            self._fig.text(self._Style_cont.Get_cmsTextPosition().getX() + 0.08, self._Style_cont.Get_cmsTextPosition().getY(), self._Style_cont.Get_additional_text(), va='bottom', ha='left', color=self._Style_cont.Get_annotation_text_color(), size=10, style = 'italic')

    def _SavePlot(self, out_name):
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

    def _Draw_main(self):
        self._fig = plt.figure(figsize=(7, 6), dpi=100, facecolor=self._Style_cont.Get_bg_color())
        self._ax1 = plt.subplot2grid((100,100), (0, self._y_starting_point), rowspan = 100 - self._x_projection_size, colspan = 100 - self._y_projection_size, axisbg = self._Style_cont.Get_bg_color())
        self._ax1.spines['bottom'].set_color(self._Style_cont.Get_spine_color())
        self._ax1.spines['bottom'].set_linewidth(self._Style_cont.Get_spine_line_width())
        self._ax1.spines['top'].set_color(self._Style_cont.Get_spine_color())
        self._ax1.spines['top'].set_linewidth(self._Style_cont.Get_spine_line_width())
        self._ax1.spines['left'].set_color(self._Style_cont.Get_spine_color())
        self._ax1.spines['left'].set_linewidth(self._Style_cont.Get_spine_line_width())
        self._ax1.spines['right'].set_color(self._Style_cont.Get_spine_color())
        self._ax1.spines['right'].set_linewidth(self._Style_cont.Get_spine_line_width())
        self._plotted_hist = rplt.imshow(self._hist, axes = self._ax1, cmap=mpl.cm.YlGnBu)
        # self._plotted_hist.set_interpolation('nearest')

        self._ax1.tick_params(axis = 'y', colors = self._Style_cont.Get_tick_color())
        self._ax1.tick_params(axis = 'x', colors = self._Style_cont.Get_tick_color())

        if self._x_projection_size > 0:
            self._Draw_x_projection()
            plt.setp(self._ax1.get_xticklabels(), visible = False)
        else:
            self._ax1.set_xlabel(self._hist.xaxis.GetTitle(), color = self._Style_cont.Get_label_text_color(), position = (1., -0.1), va = 'top', ha = 'right', size = self._Style_cont.Get_axis_title_font_size(), weight = 'medium')

        if self._y_projection_size > 0:
            self._Draw_y_projection()
            plt.setp(self._ax1.get_yticklabels(), visible = False)
        else:
            self._ax1.yaxis.set_label_coords(self._Style_cont.Get_y_label_offset(),0.9)
            self._ax1.set_ylabel(self._hist.yaxis.GetTitle(), color=self._Style_cont.Get_label_text_color(), va='top', ha='left', size = self._Style_cont.Get_axis_title_font_size(), weight = 'medium')

        if self._Style_cont.Get_ymin() != -1 and self._Style_cont.Get_ymax() != -1:
            self._ax1.set_ylim(ymin = self._Style_cont.Get_ymin(), ymax = self._Style_cont.Get_ymax())
        if self._Style_cont.Get_xmin() != -1 and self._Style_cont.Get_xmax() != -1:
            self._ax1.set_xlim(xmin = self._Style_cont.Get_xmin(), xmax = self._Style_cont.Get_xmax())

        plt.subplots_adjust(left = .10, bottom = .08, right =  .85, top = .95, wspace = .2, hspace = .0)

        self._Write_additional_text()

        self._z_axis = self._fig.add_axes([0.87, 0.08, 0.05, 0.87])
        self._fig.colorbar(self._plotted_hist, cax = self._z_axis)

        self._z_axis.set_ylabel(self._hist.zaxis.GetTitle(), color = self._Style_cont.Get_label_text_color(), position = (-0.1, 0.9), va = 'top', ha = 'left', size = self._Style_cont.Get_axis_title_font_size(), weight = 'medium')

        self._z_axis.tick_params(axis = 'y', colors = self._Style_cont.Get_tick_color())

    def _Calc_x_projection(self):
        if self._Style_cont.Get_content() == 'Efficiencies':
            pass
        else:
            self._x_projection = self._hist.ProjectionX('test_name_x', 0, -1, 'e')

    def _Draw_x_projection(self):
        self._Calc_x_projection()
        self._ax2 = plt.subplot2grid((100,100), (100 - self._x_projection_size, self._y_starting_point), rowspan = self._x_projection_size, colspan = 100 - self._y_projection_size, axisbg = self._Style_cont.Get_bg_color(), sharex = self._ax1)
        rplt.hist(self._Get_rootpy_hist1d(self._x_projection), axes = self._ax2)
        
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

        if self._y_projection_size > 0:
            self._ax2.yaxis.set_major_locator(mticker.MaxNLocator(nbins=3, prune='upper'))

        self._ax2.set_xlabel(self._hist.xaxis.GetTitle(), color = self._Style_cont.Get_label_text_color(), position = (1., -0.1), va = 'top', ha = 'right', size = self._Style_cont.Get_axis_title_font_size(), weight = 'medium')

    def _Calc_y_projection(self):
        if self._Style_cont.Get_content() == 'Efficiencies':
            pass
        else:
            self._y_projection = self._hist.ProjectionY('test_name_y', 0, -1, 'e')

    def _Draw_y_projection(self):
        self._Calc_y_projection()
        self._ax3 = plt.subplot2grid((100,100), (0, 0), rowspan = 100 - self._x_projection_size, colspan = self._y_projection_size, axisbg = self._Style_cont.Get_bg_color(), sharey = self._ax1)
        x_i = []
        x_ii = []
        for i in range(0, self._y_projection.GetNbinsX()):
            x_i.append(self._y_projection.GetBinContent(i))
            x_ii.append(i*self._y_projection.GetBinWidth(i))
        plt.hist(x_ii, weights = x_i, histtype = 'step', orientation='horizontal', axes = self._ax3)
        for label in self._ax3.xaxis.get_ticklabels():
            label.set_rotation(270)

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

        if self._x_projection_size > 0:
            self._ax3.xaxis.set_major_locator(mticker.MaxNLocator(nbins=3, prune='upper'))

        self._ax3.yaxis.set_label_coords(self._Style_cont.Get_y_label_offset() - 0.3, 0.9)
        self._ax3.set_ylabel(self._hist.yaxis.GetTitle(), color=self._Style_cont.Get_label_text_color(), va='top', ha='left', size = self._Style_cont.Get_axis_title_font_size(), weight = 'medium')

    def _Get_rootpy_hist2d(self):
        dummy_hist = Hist2D(self._hist.GetNbinsX(), self._hist.GetXaxis().GetXmin(), self._hist.GetXaxis().GetXmax(), self._hist.GetNbinsY(), self._hist.GetYaxis().GetXmin(), self._hist.GetYaxis().GetXmax())
        for i in range(0, self._hist.GetNbinsX()):
            for j in range(0, self._hist.GetNbinsY()):
                dummy_hist.SetBinContent(i,j,self._hist.GetBinContent(i,j))
                dummy_hist.SetBinError(i,j,self._hist.GetBinError(i,j))
        dummy_hist.xaxis.SetTitle('$' + self._hist.GetXaxis().GetTitle().replace('#','\\') + '$')
        dummy_hist.xaxis.SetTitle(dummy_hist.xaxis.GetTitle().replace('$$','$'))
        dummy_hist.yaxis.SetTitle('$' + self._hist.GetYaxis().GetTitle().replace('#','\\') + '$')
        dummy_hist.yaxis.SetTitle(dummy_hist.yaxis.GetTitle().replace('$$','$'))
        if self._hist.GetZaxis().GetTitle() != '' or self._hist.GetZaxis().GetTitle() == '$$':
            dummy_hist.zaxis.SetTitle('$' + self._hist.GetZaxis().GetTitle().replace('#','\\') + '$')
        else:
            dummy_hist.zaxis.SetTitle('$\epsilon$')
        if dummy_hist.zaxis.GetTitle() == '$$':
            dummy_hist.zaxis.SetTitle('$\epsilon$')
        else:
            dummy_hist.zaxis.SetTitle(dummy_hist.zaxis.GetTitle().replace('$$','$'))
        return dummy_hist

    def _Get_rootpy_hist1d(self, temp_hist):
        dummy_hist = Hist(temp_hist.GetNbinsX(), temp_hist.GetXaxis().GetXmin(), temp_hist.GetXaxis().GetXmax())
        dummy_hist.SetTitle(temp_hist.GetName())
        for i in range(0, temp_hist.GetNbinsX()):
            dummy_hist.SetBinContent(i, temp_hist.GetBinContent(i))
            dummy_hist.SetBinError(i, temp_hist.GetBinError(i))
        dummy_hist.xaxis.SetTitle('$' + temp_hist.GetXaxis().GetTitle().replace('#','\\') + '$')
        dummy_hist.xaxis.SetTitle(dummy_hist.xaxis.GetTitle().replace('$$','$'))
        dummy_hist.yaxis.SetTitle('$' + temp_hist.GetYaxis().GetTitle().replace('#','\\') + '$')
        if temp_hist.GetYaxis().GetTitle() == '$$':
            dummy_hist.yaxis.SetTitle('$\epsilon$')
        else:
            dummy_hist.yaxis.SetTitle(temp_hist.GetYaxis().GetTitle().replace('$$','$'))
        return dummy_hist

def grid(x, y, z, resX=100, resY=100):
    # convert three column data to grid
    xi = linspace(min(x), max(x), resX)
    yi = linspace(min(y), max(y), resY)
    Z = griddata(x, y, z, xi, yi)
    X, Y = meshgrid(xi, yi)
    return X, Y, Z

class FixedOrderFormatter(ScalarFormatter):
    """Formats axis ticks using scientific notation with a constant order of 
    magnitude"""
    def __init__(self, order_of_mag=0, useOffset=True, useMathText=False):
        self._order_of_mag = order_of_mag
        ScalarFormatter.__init__(self, useOffset=useOffset, 
                                 useMathText=useMathText)
    def _set_orderOfMagnitude(self, range):
        """Over-riding this to avoid having orderOfMagnitude reset elsewhere"""
        self.orderOfMagnitude = self._order_of_mag


def contourplot(xlist, ylist, zlist, levels=[], labels=[], overlay=True, resolution=100, interactive=True, showPoints=True, logz=False, filename="", contoursAlpha=0.5, contourLineWidth=4, CMSLabel=True, CMSLabelCoords=[5,87,5,84], xlabel="x", ylabel="y", zlabel="z", xyrange=[0,100,0,100], LegendLocation="upper right", LegendTitle="contourplot",LegendColumns=1):

    X,Y,Z = grid(xlist, ylist, zlist, resX=resolution, resY=resolution)
    fig = plt.figure(figsize=(16,10))

    if showPoints:
        datapoints = plt.plot(xlist, ylist, "x", alpha=contoursAlpha)

    if logz:
        TickLocator=mpl.ticker.LogLocator()
    else:
        TickLocator=mpl.ticker.LinearLocator()

    contours = []
    colors = ['r','g',"b","c","m","y","k"]
    for i in range(len(levels)):
        colors += colors
        contours.append( plt.contour(X,Y,Z, [levels[i]], colors=colors.pop(0), linewidth=contourLineWidth, locator=TickLocator) )

    if overlay:
        plt.contourf(X,Y,Z, alpha=contoursAlpha, locator=TickLocator, cmap=mpl.cm.Blues)
        cbar = plt.colorbar()
        cbar.set_label(zlabel)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xlim(xyrange[0],xyrange[1])
    plt.ylim(xyrange[2],xyrange[3])

    # Do the legend:
    lines = []
    for i in range(len(levels)):
        lines.append(contours[i].collections[0])
    if len(levels) != 0:
        plt.legend(lines, labels, loc=LegendLocation, title=LegendTitle, ncol=LegendColumns)

    # CMS Annotations:
    if CMSLabel:
        plt.text(CMSLabelCoords[0],CMSLabelCoords[1], "CMS", fontweight='semibold')
        plt.text(CMSLabelCoords[2],CMSLabelCoords[3], "Simulation", style="italic")
    
    # Show and/or save:
    if filename != "":  plt.savefig(filename)
    if interactive:     plt.show()
    return contours


def main():

    font = {'size': 15}
    mpl.rc('font', **font)
    mpl.rcParams['lines.linewidth'] = 2

    # do a sample plot:   
    xlist = [0,   0,    100,  100]
    ylist = [0,   100,  0,    100]
    zlist = [500, 1e-5, 1e-5, 1e-8]
    
    levels = [10,  20,  30]
    labels = ["a", "b", "c"]

    contourplot(xlist, ylist, zlist, levels, labels, logz=True, LegendLocation="lower left", filename='test.pdf')


if __name__ == "__main__":
    main()
