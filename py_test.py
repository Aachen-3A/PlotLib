#!/bin/env python

import ROOT
import numpy as np
from rootpy.plotting import Hist, Hist2D, Hist3D, HistStack, Legend, Canvas, Graph
from rootpy.plotting.style import get_style, set_style
from rootpy.interactive import wait
import random
import rootpy.plotting.root2matplotlib as rplt
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

# create a simple 1D histogram with 10 constant-width bins between 0 and 1
h_simple = Hist(10, 0, 1)
print h_simple.name

# If the name is not specified, a UUID is used so that ROOT never complains
# about two histograms having the same name.
# Alternatively you can specify the name (and the title or any other style
# attributes) in the constructor:
h_simple = Hist(10, -4, 12, name='my hist', title='Some Data',
                drawstyle='hist',
                legendstyle='F',
                fillstyle='/')

# fill the histogram
for i in xrange(1000):
    # all ROOT CamelCase methods are aliased by equivalent snake_case methods
    # so you can call fill() instead of Fill()
    h_simple.Fill(random.gauss(4, 3))


# easily set visual attributes
h_simple.linecolor = 'blue'
h_simple.fillcolor = 'green'
h_simple.fillstyle = '/'

# attributes may be accessed in the same way
print h_simple.name
print h_simple.title
print h_simple.markersize

# plot
canvas = Canvas(width=700, height=500)
canvas.SetLeftMargin(0.15)
canvas.SetBottomMargin(0.15)
canvas.SetTopMargin(0.10)
canvas.SetRightMargin(0.05)
h_simple.Draw()

# create the legend
legend = Legend([h_simple], pad=canvas,
                header='Header',
                leftmargin=0.05,
                rightmargin=0.5)
legend.Draw()

# 2D and 3D histograms are handled in the same way
# the constructor arguments are repetitions of #bins, left bound, right bound.
h2d = Hist2D(10, 0, 1, 50, -40, 10, name='2d hist')
h3d = Hist3D(3, -1, 4, 10, -1000, -200, 2, 0, 1, name='3d hist')

# variable-width bins may be created by passing the bin edges directly:
h1d_variable = Hist([1, 4, 10, 100])
h2d_variable = Hist2D([2, 4, 7, 100, 200], [-100, -50, 0, 10, 20])
h3d_variable = Hist3D([1, 3, 10], [20, 50, 100], [-10, -5, 10, 20])

# variable-width and constant-width bins can be mixed:
h2d_mixed = Hist2D([2, 10, 30], 10, 1, 5)

# wait for you to close all open canvases before exiting
# wait() will have no effect if ROOT is in batch mode:
# ROOT.gROOT.SetBatch(True)
wait()

# set the random seed
ROOT.gRandom.SetSeed(42)
np.random.seed(42)

# points
x = np.sort(np.random.random(10)) * 3500
y = np.random.random(10)

# set style for ROOT
set_style('CMSTDR')

# create graph
graph = Graph(x.shape[0])
for i, (xx, yy) in enumerate(zip(x, y)):
    graph.SetPoint(i, xx, yy)

# set visual attributes
graph.linecolor = 'blue'
graph.markercolor = 'blue'
graph.xaxis.SetTitle("E_{T} [GeV]")
graph.yaxis.SetTitle("d#sigma_{jet}/dE_{T,jet} [fb/GeV]")
graph.xaxis.SetRangeUser(0, 3500)
graph.yaxis.SetRangeUser(0, 1)

# plot with ROOT
canvas = Canvas()
graph.Draw("APL")

label = ROOT.TText(0.4, 0.8, "ROOT")
label.SetTextFont(43)
label.SetTextSize(25)
label.SetNDC()
label.Draw()
canvas.Modified()
canvas.Update()

# plot with matplotlib

def plot_with_matplotlib():
    fig, axes = plt.subplots()

    axes.plot(x, y, 'o-', markeredgewidth=0)
    axes.set_xlabel(r"$E_T$ [GeV]",
                    horizontalalignment="right", x=1, labelpad=20)
    axes.set_ylabel(r"$d\sigma_{jet}/dE_{T,jet}$ [fb/GeV]",
                    horizontalalignment="right", y=1, labelpad=32)
    axes.set_xlim(0, 3500)
    axes.set_ylim(0, 1)

    return fig, axes

# plot without style
fig1, axes1 = plot_with_matplotlib()
axes1.text(0.4, 0.8, 'matplotlib (no style)',
           verticalalignment='center', horizontalalignment='center',
           transform=axes1.transAxes, fontsize=20)

# plot with ATLAS style
set_style('ATLAS', mpl=True)
fig2, axes2 = plot_with_matplotlib()
axes2.text(0.4, 0.8, 'matplotlib',
           verticalalignment='center', horizontalalignment='center',
           transform=axes2.transAxes, fontsize=20)
axes2.xaxis.set_minor_locator(AutoMinorLocator())
axes2.yaxis.set_minor_locator(AutoMinorLocator())

if not ROOT.gROOT.IsBatch():
    plt.show()

# wait for you to close the canvas before exiting
wait(True)
