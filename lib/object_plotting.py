#!/bin/env python

import sys
import matplotlib
import numpy as np
import rootpy
import ROOT as r

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

import rootpy.plotting.root2matplotlib as rplt

def plot_upper_cut(axis1, value, cut_arrow_position = 0.5):
    print('will now plot an upper cut marker at %f'%value)

    cut_color = 'y'
    cut_level = 4

    ymin, ymax = axis1.get_ylim()
    height = ymax - ymin

    xmin, xmax = axis1.get_xlim()
    width = xmax - xmin

    axis1.axvline(x = value, color = cut_color, zorder = cut_level, linewidth = 2)
    axis1.arrow((value - xmin) / width, cut_arrow_position, -0.02, 0, head_width = 0.1, head_length = 0.02, fc = cut_color, ec = cut_color, zorder = cut_level, linewidth = 2, transform = axis1.transAxes)

def plot_lower_cut(axis1, value, cut_arrow_position = 0.5):
    print('will now plot an lower cut marker at %f'%value)

    cut_color = 'y'
    cut_level = 4

    ymin, ymax = axis1.get_ylim()
    height = ymax - ymin

    xmin, xmax = axis1.get_xlim()
    width = xmax - xmin

    axis1.axvline(x = value, color = cut_color, zorder = cut_level, linewidth = 2)
    axis1.arrow((value - xmin) / width, cut_arrow_position, 0.02, 0, head_width = 0.1, head_length = 0.02, fc = cut_color, ec = cut_color, zorder = cut_level, linewidth = 2, transform = axis1.transAxes)

def plot_both_cuts(axis1, value1, value2, include = True, cut_arrow_position = 0.5):
    if include:
        if value1 > value2:
            plot_lower_cut(axis1, value2, cut_arrow_position = cut_arrow_position)
            plot_upper_cut(axis1, value1, cut_arrow_position = cut_arrow_position)
        else:
            plot_lower_cut(axis1, value1, cut_arrow_position = cut_arrow_position)
            plot_upper_cut(axis1, value2, cut_arrow_position = cut_arrow_position)
    else:
        if value1 > value2:
            plot_lower_cut(axis1, value1, cut_arrow_position = cut_arrow_position)
            plot_upper_cut(axis1, value2, cut_arrow_position = cut_arrow_position)
        else:
            plot_lower_cut(axis1, value2, cut_arrow_position = cut_arrow_position)
            plot_upper_cut(axis1, value1, cut_arrow_position = cut_arrow_position)

def gaussian(x, norm, mu, sig):
    return norm * np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def plot_gauss_fit(axis1, histo, xmin = 0, xmax = 0, color = 'black', write_results = False):
    print('Now fitting a gaus to %s and plotting it'%histo.GetTitle())
    if xmin == 0 and xmax == 0:
        xmin, xmax = axis1.get_xlim()

    try:
        fit_res = histo.Fit('gaus', 'N0S', '', xmin, xmax)
        dummy = fit_res.Parameter(0)
    except rootpy.ROOTError as bla:
        print('\n\terror in the Chi2 fitting: ')
        print('\t'+str(bla))
        print('\ttrying again with binned likelihood method\n')
        try:
            fit_res = histo.Fit('gaus', 'WLN0S', '', xmin, xmax)
            dummy = fit_res.Parameter(0)
        except rootpy.ROOTError as bla:
            print('\n\terror in the binned likelihood fitting: ')
            print('\t' + str(bla))
            print('\tThis doesn\'t seem to work either, I think you have to fix something!\n')
            sys.exit(42)
        except:
            print('An unexpected error occured in the binned likelihood fitting:')
            print(sys.exc_info()[0], sys.exc_info()[1])
            sys.exit(42)
    except:
        print('An unexpected error occured in the Chi2 fitting:')
        print(sys.exc_info()[0], sys.exc_info()[1])
        sys.exit(42)

    x = np.linspace(xmin, xmax, num = 500)
    y = gaussian(x, fit_res.Parameter(0), fit_res.Parameter(1), fit_res.Parameter(2))
    line, = axis1.plot(x, y, '-', linewidth=1, zorder = 0, color = color)

    if write_results:
        text = r'''Fit of $A \cdot \exp\left(\frac{(x - B)^2}{2 \cdot C^2}\right)$
    
$ \chi^2/Ndf = %.2f/%.0f = %.2f$

$ A = %.2f \pm %.2f$
$ B = %.5f \pm %.5f$
$ C = %.5f \pm %.5f$'''%(fit_res.Chi2(),
                                 fit_res.Ndf(),
                                 fit_res.Chi2()/fit_res.Ndf(),
                                 fit_res.Parameter(0),fit_res.ParError(0),
                                 fit_res.Parameter(1),fit_res.ParError(1),
                                 fit_res.Parameter(2),fit_res.ParError(2))
    
        plt.text(0.3, 400, text)

def eff_function(x, par):
    return par[0] + par[1] / (x[0] + par[2]) + par[3] * x[0]

def plot_efficiency_fit(axis1, histo, xmin = 0, xmax = 0, startvals = [], plottrange = [], color = 'black'):
    print('Now fitting a efficiency curve to %s and plotting it'%histo.GetTitle())
    if xmin == 0 and xmax == 0:
        xmin, xmax = axis1.get_xlim()

    func = r.TF1('func',eff_function,xmin,xmax,4)

    if startvals != []:
        func.SetParameter(0, startvals[0])
        func.SetParameter(1, startvals[1])
        func.SetParameter(2, startvals[2])
        func.SetParameter(3, startvals[3])

    try:
        fit_res = histo.Fit(func, 'N0S', '', xmin, xmax)
        dummy = fit_res.Parameter(0)
    except rootpy.ROOTError as bla:
        print('\n\terror in the Chi2 fitting: ')
        print('\t'+str(bla))
        print('\ttrying again with binned likelihood method\n')
        try:
            fit_res = histo.Fit(func, 'WLN0S', '', xmin, xmax)
            dummy = fit_res.Parameter(0)
        except rootpy.ROOTError as bla:
            print('\n\terror in the binned likelihood fitting: ')
            print('\t' + str(bla))
            print('\tThis doesn\'t seem to work either, I think you have to fix something!\n')
            sys.exit(42)
        except:
            print('An unexpected error occured in the binned likelihood fitting:')
            print(sys.exc_info()[0], sys.exc_info()[1])
            sys.exit(42)
    except:
        print('An unexpected error occured in the Chi2 fitting:')
        print(sys.exc_info()[0], sys.exc_info()[1])
        sys.exit(42)

    print('Chi2/Ndf: %.2f/%.0f = %.2f'%(fit_res.Chi2(),fit_res.Ndf(),fit_res.Chi2()/fit_res.Ndf()))

    if plottrange != []:
        x = np.linspace(plottrange[0], plottrange[1])
    else:
        x = np.linspace(xmin, xmax)
    y = eff_function([x], [fit_res.Parameter(0), fit_res.Parameter(1), fit_res.Parameter(2), fit_res.Parameter(3)])
    line, = axis1.plot(x, y, '-', linewidth=2, zorder = 0, color = color)

    text = r'''Fit of $A + \frac{B}{C + M} + D \cdot M$
    
$ \chi^2/Ndf = %.2f/%.0f = %.2f$

$ A = %.4f \pm %.4f$
$ B = %.1f \pm %.1f$
$ C = %f \pm %f$
$ D = %f \pm %f$'''%(fit_res.Chi2(),
                             fit_res.Ndf(),
                             fit_res.Chi2()/fit_res.Ndf(),
                             fit_res.Parameter(0),fit_res.ParError(0),
                             fit_res.Parameter(1),fit_res.ParError(1),
                             fit_res.Parameter(2),fit_res.ParError(2),
                             fit_res.Parameter(3),fit_res.ParError(3))

    plt.text(2000, 0.25, text)

def res_function(x, par):
    return par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0]

def plot_resolution_fit(axis1, histo, xmin = 0, xmax = 0, startvals = [], plottrange = [], color = 'black'):
    print('Now fitting a resolution curve to %s and plotting it'%histo.GetTitle())
    if xmin == 0 and xmax == 0:
        xmin, xmax = axis1.get_xlim()

    func = r.TF1('func',res_function,xmin,xmax,4)

    if startvals != []:
        func.SetParameter(0, startvals[0])
        func.SetParameter(1, startvals[1])
        func.SetParameter(2, startvals[2])
        func.SetParameter(3, startvals[3])

    try:
        fit_res = histo.Fit(func, 'N0S', '', xmin, xmax)
        dummy = fit_res.Parameter(0)
    except rootpy.ROOTError as bla:
        print('\n\terror in the Chi2 fitting: ')
        print('\t'+str(bla))
        print('\ttrying again with binned likelihood method\n')
        try:
            fit_res = histo.Fit(func, 'WLN0S', '', xmin, xmax)
            dummy = fit_res.Parameter(0)
        except rootpy.ROOTError as bla:
            print('\n\terror in the binned likelihood fitting: ')
            print('\t' + str(bla))
            print('\tThis doesn\'t seem to work either, I think you have to fix something!\n')
            sys.exit(42)
        except:
            print('An unexpected error occured in the binned likelihood fitting:')
            print(sys.exc_info()[0], sys.exc_info()[1])
            sys.exit(42)
    except:
        print('An unexpected error occured in the Chi2 fitting:')
        print(sys.exc_info()[0], sys.exc_info()[1])
        sys.exit(42)

    print('Chi2/Ndf: %.2f/%.0f = %.2f'%(fit_res.Chi2(),fit_res.Ndf(),fit_res.Chi2()/fit_res.Ndf()))

    if plottrange != []:
        x = np.linspace(plottrange[0], plottrange[1])
    else:
        x = np.linspace(xmin, xmax)
    y = res_function([x], [fit_res.Parameter(0), fit_res.Parameter(1), fit_res.Parameter(2), fit_res.Parameter(3)])
    line, = axis1.plot(x, y, '-', linewidth=2, zorder = 0, color = color)

    text = r'''Fit of $A + B \cdot M + C \cdot M^2 + D \cdot M^3$
    
$ \chi^2/Ndf = %.2f/%.0f = %.2f$

$ A = %f \pm %f$
$ B = %f \pm %f$
$ C = %f \pm %f$
$ D = %f \pm %f$'''%(fit_res.Chi2(),
                             fit_res.Ndf(),
                             fit_res.Chi2()/fit_res.Ndf(),
                             fit_res.Parameter(0),fit_res.ParError(0),
                             fit_res.Parameter(1),fit_res.ParError(1),
                             fit_res.Parameter(2),fit_res.ParError(2),
                             fit_res.Parameter(3),fit_res.ParError(3))

    plt.text(2500, 0.015, text)
    return func

def draw_lines(axis):
    for idx_st in range(1,5):
        nSectors = 12
        if idx_st == 4: nSectors = 14
        for idx_wh in range(-1,3):
            xline = (idx_st - 1)*60 + (idx_wh + 2)*nSectors
            axis.axvline(x = xline, color = 'slategray', linewidth = 1, linestyle = ':')

    for idx in range(1,4):
        xline = idx*60
        axis.axvline(x = xline, color = 'dimgray', linewidth = 2, linestyle = '--')

    for idx in range(1,5):
        xlabel = (idx - 1)*60 + 20
        a,b = axis.get_ylim()
        ylabel = a + 0.63*(b -a)

        strSt = "MB%d" % idx
        axis.text(xlabel,ylabel, strSt, va='bottom', ha='left', size=12)

    stations = (1,2,3,4)
    wheels = (-2,-1,0,1,2)
    bin_labels = [''] * 250
    for st in stations:
        for wh in wheels:
            nSectors = 12
            if st == 4: nSectors = 14 
            for sec in range(1,nSectors+1):
                if sec == 1:
                    label = "Wheel %d" % wh
                    if wh == -2: label += " MB%d" % st
                    binHistoNew = (st - 1)*60 + (wh + 2)*nSectors + sec
                    bin_labels[binHistoNew+2] = label

    axis.xaxis.set_major_locator(mticker.MaxNLocator(nbins=len(bin_labels)))
    xtickNames = plt.setp(axis, xticklabels=bin_labels)
    plt.setp(xtickNames, rotation=-45, fontsize=8, va='top', ha='left')
    axis.tick_params('x', length=0, width=0, which='major')
