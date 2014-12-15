# Creates 2D plots of L and std(L) versus n and N, where N is the
# number of bins, and n the number of objects per bin.

import sys
import os
#sys.path.append(sys.argv[1])
import numpy as np
import pylab as pl
import argparse
from data_analysers.BiasAnalyser import *
from bias_methods import *    

parser = argparse.ArgumentParser(description='Calculate labeling bias of a data-set.')

parser.add_argument("table_file",
                    help = "table with variables and labels.")
parser.add_argument("label_name", help = "name of the column containing the labels.")
parser.add_argument("--number_objects", metavar = "N", 
                    default = [10, 20, 50], type = np.mat,
                    help = "number of objects per bin to be used for calculating the bias.")
parser.add_argument("--int_pars", default = ["petroRad_r_kpc","absPetroMag_r"], 
                    nargs = 2,
                    help = "Name of the columns used as intrinsic parameters.")
parser.add_argument("--obs_pars", default = ["z", "corrMag_r", "petroRad_r_psf"],
                    help = "Name of the columns used as observable parameters.")
parser.add_argument("--N_iter", default = 5, type = int,
                    help = "Number of calculations of L to calculate means and standard deviations.")
parser.add_argument("--bins_int", default = [[2, 2], [5, 5], [10, 10], 
                                         [15, 15], [20, 20]],
                    type = np.mat,
                    help = "Bins in intrinsic and observable parameters.")
parser.add_argument("--bins_obs", default = [2, 5, 10, 20],
                    type = np.mat,
                    help = "Bins in intrinsic and observable parameters.")
parser.add_argument("--out_path", default = "./",
                    help = "Path where to leave results.")
parser.add_argument("--plot_N_objs",
                    help = "N_objs for plot L versus bins")

args = parser.parse_args()

field_int = args.int_pars
fields_obs = args.obs_pars
classf = args.label_name
N_iter = args.N_iter
N_objs = np.array(args.number_objects, dtype = int)[0]
bins_int = np.array(args.bins_int, dtype = int)
bins_obs = np.array(args.bins_obs, dtype = int)
out_path = args.out_path
bins = np.zeros((bins_int.shape[0] * len(bins_obs), 3), dtype = int)
for i in range (bins_int.shape[0]):
    for j in range (len(bins_obs)):
        #bins[i*len(bins_obs) + j][:] = [bins_int[i][0], bins_int[i][1], bins_obs[j]]
        bins[i + j * bins_int.shape[0]][:] = [bins_int[i][0], bins_int[i][1], bins_obs[j]]

print bins

dataAn = BiasAnalyser (args.table_file)

fields_cut = []
cuts = []

# Calculate L versus N and save in files
Ls, Ls_std = saveLs (dataAn, field_int[0], field_int[1], fields_obs, classf, 
                     fields_cut, cuts, bins, N_objs, N_iter, out_path)

print Ls

y_ticks = []
for i in range (bins.shape[0]):
    y_ticks.append(str(bins[i][0]) + "x" + str(bins[i][1]) + "x" + str(bins[i][2])
                   + " = " + str (bins[i].prod()))
def plotImage (Ls, bins, N_objs, filename):
    # determine ticks for bins
    
    pl.clf()
    pl.imshow(Ls.transpose(), interpolation = "nearest")
    pl.gca().invert_yaxis()
    pl.colorbar()
    pl.xticks (np.arange(len(N_objs)), N_objs, rotation = "90")
    pl.yticks (np.arange(bins.shape[0]), y_ticks)
    pl.savefig(filename)

plotImage (Ls, bins, N_objs, "L_vs_N")
plotImage (Ls_std, bins, N_objs, "L_vs_N_std")

pl.clf()
for i in range (bins.shape[0]):
    leg = (str(bins[i][0]) + "x" + str(bins[i][1]) + "x" + str(bins[i][2])
           + " = " + str (bins[i].prod()))

    pl.errorbar(N_objs[Ls[:, i]!=0], Ls[:, i][Ls[:, i]!=0], 
                yerr = Ls_std[:, i][Ls[:, i]!=0], label = leg)
pl.xscale("log")
pl.legend()
pl.savefig ("L_vs_N_objs")

# plotting L versus bins for "plot_N_objs"
if args.plot_N_objs:
    pl.clf ()
    pl.errorbar(np.arange(Ls.shape[1]), Ls[N_objs == args.plot_N_objs], 
                yerr = Ls_std[N_objs == args.plot_N_objs])
    pl.xticks (np.arange(bins.shape[0]), y_ticks, rotation = "60", size = "small")
    pl.subplots_adjust(bottom=0.25)
    pl.margins(0.05)
    pl.savefig ("L_vs_bins_" + str(args.plot_N_objs))
