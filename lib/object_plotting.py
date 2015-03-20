#!/bin/env python

import sys
import matplotlib
import numpy as np
import rootpy

import matplotlib.pyplot as plt

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

def plot_gauss_fit(axis1, histo, xmin = 0, xmax = 0):
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

    x = np.linspace(xmin, xmax)
    y = gaussian(x, fit_res.Parameter(0), fit_res.Parameter(1), fit_res.Parameter(2))
    line, = axis1.plot(x, y, '-', linewidth=1, zorder = 0)

