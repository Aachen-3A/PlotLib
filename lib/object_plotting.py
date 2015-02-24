#!/bin/env python

import sys
import matplotlib

import matplotlib.pyplot as plt

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
