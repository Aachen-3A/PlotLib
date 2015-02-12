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

import rootpy.plotting.root2matplotlib as rplt

from rootpy.plotting import Hist2D, Hist

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

    def make_plot(self,out_name):
        self._hist = self._Get_rootpy_hist2d()
        self._Draw_main()
        self._SavePlot(out_name)

    def Add_x_projection(self, size = 15):
        self._x_projection_size = size
        self._x_starting_point = size

    def Add_y_projection(self, size = 15):
        self._y_projection_size = size
        self._y_starting_point = size

    def Set_axis(self, logx = False, logy = False, ymin = -1, ymax = -1, xmin = -1, xmax = -1, grid = False):
        self._Style_cont.Set_axis(logx = logx, logy = logy, ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax, grid = grid)

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
        ax1 = plt.subplot2grid((100,100), (0, self._y_starting_point), rowspan = 100 - self._x_projection_size, colspan = 100 - self._y_projection_size, axisbg = self._Style_cont.Get_bg_color())
        # ax1.set_zscale('log')
        self._plotted_hist = rplt.imshow(self._hist, axes = ax1, cmap=mpl.cm.brg)

        if self._x_projection_size > 0:
            self._Draw_x_projection(ax1)
            plt.setp(ax1.get_xticklabels(), visible = False)

        if self._y_projection_size > 0:
            self._Draw_y_projection(ax1)
            plt.setp(ax1.get_yticklabels(), visible = False)

        if self._Style_cont.Get_ymin() != -1 and self._Style_cont.Get_ymax() != -1:
            ax1.set_ylim(ymin = self._Style_cont.Get_ymin(), ymax = self._Style_cont.Get_ymax())
        if self._Style_cont.Get_xmin() != -1 and self._Style_cont.Get_xmax() != -1:
            ax1.set_xlim(xmin = self._Style_cont.Get_xmin(), xmax = self._Style_cont.Get_xmax())

        plt.subplots_adjust(left = .10, bottom = .08, right =  .85, top = .95, wspace = .2, hspace = .0)

        cbar_ax = self._fig.add_axes([0.87, 0.08, 0.05, 0.87])
        self._fig.colorbar(self._plotted_hist, cax = cbar_ax)

        plt.show()

    def _Draw_x_projection(self, axis1):
        ax2 = plt.subplot2grid((100,100), (100 - self._x_projection_size, self._y_starting_point), rowspan = self._x_projection_size, colspan = 100 - self._y_projection_size, axisbg = self._Style_cont.Get_bg_color(), sharex = axis1)
        self._x_projection = self._hist.ProjectionX('test_name_x', 0, -1, 'e')
        rplt.hist(self._Get_rootpy_hist1d(self._x_projection), axes = ax2)

    def _Draw_y_projection(self, axis1):
        ax3 = plt.subplot2grid((100,100), (0, 0), rowspan = 100 - self._x_projection_size, colspan = self._y_projection_size, axisbg = self._Style_cont.Get_bg_color(), sharey = axis1)
        self._y_projection = self._hist.ProjectionY('test_name_y', 0, -1, 'e')
        x_i = []
        x_ii = []
        for i in range(0, self._y_projection.GetNbinsX()):
            x_i.append(self._y_projection.GetBinContent(i))
            x_ii.append(i*self._y_projection.GetBinWidth(i))
        plt.hist(x_ii, weights = x_i, histtype = 'step', orientation='horizontal', axes = ax3)
        for label in ax3.xaxis.get_ticklabels():
            label.set_rotation(270)

    def _Get_rootpy_hist2d(self):
        dummy_hist = Hist2D(self._hist.GetNbinsX(), self._hist.GetXaxis().GetXmin(), self._hist.GetXaxis().GetXmax(), self._hist.GetNbinsY(), self._hist.GetYaxis().GetXmin(), self._hist.GetYaxis().GetXmax())
        for i in range(0, self._hist.GetNbinsX()):
            for j in range(0, self._hist.GetNbinsY()):
                dummy_hist.SetBinContent(i,j,self._hist.GetBinContent(i,j))
                dummy_hist.SetBinError(i,j,self._hist.GetBinError(i,j))
        return dummy_hist

    def _Get_rootpy_hist1d(self, temp_hist):
        dummy_hist = Hist(temp_hist.GetNbinsX(), temp_hist.GetXaxis().GetXmin(), temp_hist.GetXaxis().GetXmax())
        for i in range(0, temp_hist.GetNbinsX()):
            dummy_hist.SetBinContent(i, temp_hist.GetBinContent(i))
            dummy_hist.SetBinError(i, temp_hist.GetBinError(i))
        return dummy_hist

def grid(x, y, z, resX=100, resY=100):
    # convert three column data to grid
    xi = linspace(min(x), max(x), resX)
    yi = linspace(min(y), max(y), resY)
    Z = griddata(x, y, z, xi, yi)
    X, Y = meshgrid(xi, yi)
    return X, Y, Z


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
