# Creates 2D plots of L and std(L) versus n and N, where N is the
# number of bins, and n the number of objects per bin.

import sys
import os
#sys.path.append(sys.argv[1])
import numpy as np
import pylab as pl
import argparse
from data_analysers.BiasAnalyser import *
import bias_methods as bm   

parser = argparse.ArgumentParser(description='Calculate labeling bias of a data-set.')

parser.add_argument("table_file",
                    help = "table with variables and labels.")
parser.add_argument("label_name", nargs = "+",
                    help = "name of the column containing the labels.")
parser.add_argument("--number_objects", metavar = "N", 
                    default = [10, 20, 50], type = np.array,
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
                    help = "Bins in intrinsic and observable parameters 'b11 b12;b21 b22.")
parser.add_argument("--bins_obs", default = [2, 5, 10, 20],
                    type = np.mat,
                    help = "Bins in intrinsic and observable parameters.")
parser.add_argument("--out_path", default = "./",
                    help = "Path where to leave results.")
parser.add_argument("--plot_N_objs", type = int,
                    help = "N_objs for plot L versus bins")
parser.add_argument("--plot_N_bins_int", type = np.mat,
                    help = "number of bins to use in the L versus N_objs plot")
parser.add_argument("--pbb_thresholds", nargs = "+",
                    help = "Threshold on the probabilities")
parser.add_argument("--no_zeros", action='store_const', const = True,
                    help = "Do not consider labels that don't match the pbb. thresholds")

args = parser.parse_args()

field_int = args.int_pars
fields_obs = args.obs_pars
classf = args.label_name
if args.pbb_thresholds:
    pbb_thresholds = np.array (args.pbb_thresholds).astype(np.float)
else:
    pbb_thresholds = args.pbb_thresholds
N_iter = args.N_iter
N_objs = np.array(args.number_objects, dtype = int)
bins_int = np.array(args.bins_int, dtype = int)
bins_obs = np.array(args.bins_obs, dtype = int)
out_path = args.out_path
bins = np.zeros((bins_int.shape[0] * len(bins_obs), 3), dtype = int)
for i in range (bins_int.shape[0]):
    for j in range (len(bins_obs)):
        #bins[i*len(bins_obs) + j][:] = [bins_int[i][0], bins_int[i][1], bins_obs[j]]
        bins[i + j * bins_int.shape[0]][:] = [bins_int[i][0], bins_int[i][1], bins_obs[j]]

# Intrinsic parameters. Currently accept exactly two.
field_r = args.int_pars[0]
field_c = args.int_pars[1]

print bins

dataAn = BiasAnalyser ()

tbdata = pf.open(args.table_file)[1].data # Open fits table file.
N_tot = len(tbdata)

# i_s is used to randomize the data.
i_s = np.arange (N_tot)
np.random.shuffle(i_s)
tbdata = tbdata[i_s]

y = bm.createLabels (tbdata, classf, pbb_thresholds)

crit_zeros = np.ones (len(y), dtype = bool)
if args.no_zeros:
    crit_zeros = (y != 0)
y = y[crit_zeros]

# Get intrinsic paraeters.
intrinsic = np.array([tbdata.field (field_r)[crit_zeros],
                      tbdata.field (field_c)[crit_zeros]]).transpose()

# Get observable parameters.
observables = []
for obs in fields_obs:
    observables.append (tbdata.field (obs)[crit_zeros])
observables = np.array(observables).transpose()

labels = np.unique(y)

# Calculate L versus N and save in files
# Ls, Ls_std = bm.saveLs (dataAn, field_int[0], field_int[1], fields_obs, classf, 
#                         fields_cut, cuts, bins, N_objs, N_iter, out_path)
print "N_objs = ", N_objs
Ls, Ls_std = bm.saveLs (dataAn, intrinsic, observables, y, labels, bins,
                        N_objs, N_iter, increasing_bias = [True, True, False],
                        out_path= out_path)

print Ls

y_ticks = []
for i in range (bins.shape[0]):
    y_ticks.append(str(bins[i][0]) + "x" + str(bins[i][1]) + "x" + str(bins[i][2])
                   + " = " + str (bins[i].prod()))
def plotImage (Ls, bins, N_objs, filename):
    # determine ticks for bins
    
    pl.clf()
    pl.subplots_adjust(bottom=0.15)
    pl.imshow(Ls.transpose(), interpolation = "nearest")
    pl.gca().invert_yaxis()
    pl.colorbar()
    pl.xticks (np.arange(len(N_objs)), N_objs, rotation = "60")
    pl.yticks (np.arange(bins.shape[0]), y_ticks)
    pl.xlabel ("number of objects per bin")
    pl.ylabel ("number of bins in intrinsic and observable parameters:\n" + 
               r"$N_R \times N_M \times N_{\mathcal{A}_{j,q}}$", 
               multialignment='center')
    pl.savefig(filename)

plotImage (Ls, bins, N_objs, "L_vs_N")
plotImage (Ls_std, bins, N_objs, "L_vs_N_std")

pl.clf()

if not np.isscalar(args.plot_N_bins_int):
    sel_bins = np.array(args.plot_N_bins_int, dtype = int)[0]
    crit = ((bins[:, 0:2] == sel_bins).prod(axis = 1))
else:
    crit = np.ones(bins.shape[0])

for i in range (bins.shape[0]):
    leg = (str(bins[i][0]) + "x" + str(bins[i][1]) + "x" + str(bins[i][2])
           + " = " + str (bins[i].prod()))

    #     pl.errorbar(N_objs[Ls[:, i]!=0], Ls[:, i][Ls[:, i]!=0], 
    #                 yerr = Ls_std[:, i][Ls[:, i]!=0], label = leg)
    # else:
    if crit[i]:
        pl.errorbar(N_objs[Ls[:, i]!=0], Ls[:, i][Ls[:, i]!=0], 
                    yerr = Ls_std[:, i][Ls[:, i]!=0], label = leg)
pl.xscale("log")
pl.xlabel ("number of objects per bin")
pl.ylabel (r"$L$")
pl.legend()
pl.savefig ("L_vs_N_objs")

# plotting L versus bins for "plot_N_objs"
if args.plot_N_objs:
    Ls_plot = Ls[N_objs == args.plot_N_objs][0]
    Ls_s_plot = Ls_std[N_objs == args.plot_N_objs][0]
    pl.clf ()
    pl.errorbar(np.arange(Ls.shape[1])[Ls_plot != 0], Ls_plot[Ls_plot != 0], 
                yerr = Ls_s_plot[Ls_plot != 0], fmt = "x")
    pl.xticks (np.arange(bins.shape[0]), y_ticks, rotation = "60", size = "small")
    pl.subplots_adjust(bottom=0.3)
    pl.margins(0.05)
    pl.xlabel ("number of bins in intrinsic and observable parameters:\n" + 
               r"$N_R \times N_M \times N_{\mathcal{A}_{j,q}}$", 
               multialignment='center')
    pl.ylabel (r"$L$")
    pl.savefig ("L_vs_bins_" + str(args.plot_N_objs))
