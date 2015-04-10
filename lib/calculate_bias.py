import sys
import os
#sys.path.append(sys.argv[1])
import numpy as np
import argparse
import bias_methods as bm

from data_analysers.BiasAnalyser import *
from data_analysers.FitsBiasAnalyser import *

parser = argparse.ArgumentParser(description='Calculate labeling bias of a data-set.')

parser.add_argument("table_file",
                    help = "table with variables and labels.")
parser.add_argument("label_name", nargs = "+",
                    help = "name of the column containing the labels.")
parser.add_argument("--number_objects", metavar = "N", type = int, default = -1, 
                    help = "number of objects per bin to be used for calculating the bias.")
parser.add_argument("--int_pars", default = ["petroRad_r_kpc","absPetroMag_r", "z"], 
                    help = "Name of the columns used as intrinsic parameters.")
parser.add_argument("--obs_pars", default = ["corrMag_r", "petroRad_r_psf"],
                    help = "Name of the columns used as observable parameters.")
parser.add_argument("--pbb_thresholds", nargs = "+",
                    help = "Threshold on the probabilities")
parser.add_argument("--no_zeros", action='store_const', const = True,
                    help = "Do not consider labels that don't match the pbb. thresholds")
parser.add_argument("--bins_obs", type = int, default = 20,
                    help = "log2 bins in intrinsic parameters.")
parser.add_argument("--log2_bins_int", type = int, default = 7,
                    help = "Bins in observable parameters.")
parser.add_argument("--N_iter", default = 5, type = int,
                    help = "Number of calculations of L to calculate means and standard deviations.")
parser.add_argument("--labels", nargs = "+",
                    help = "Labels to be used.")

args = parser.parse_args()

# Number of objects for calculating bias.
N_objs = args.number_objects

fields_int = args.int_pars
fields_beta = args.obs_pars   # Observable parameters. As many as
                              # needed.
classf = args.label_name      # Label or probability column name.
if args.pbb_thresholds:
    pbb_thresholds = np.array (args.pbb_thresholds).astype(np.float)
else:
    pbb_thresholds = args.pbb_thresholds

N_calc = args.N_iter

dataAn = BiasAnalyser ()

l2_bins_int = args.log2_bins_int # Number of bins for intrinsic pars.
bins_obs = args.bins_obs         # Number of bins for observable pars.
N_iter = args.N_iter             # Number of iterations for mean and std.
minElementsBin = 5               # Min. number of objs. per bin.
np.random.seed (0)

tbdata = pf.open(args.table_file)[1].data # Open fits table file.
N_tot = len(tbdata)

# i_s is used to randomize the data.
i_s = np.arange (N_tot)
np.random.shuffle(i_s)
tbdata = tbdata[i_s]

if args.pbb_thresholds:
    y = bm.createLabels (tbdata, classf, pbb_thresholds)
else:
    y = np.array(tbdata.field(classf[0]), dtype = int)

data_aux = tbdata.tolist()
data_aux = np.asarray(tbdata.tolist())
crit_zeros = np.ones (len(y), dtype = bool)
# if args.no_zeros:
#     crit_zeros = (y != 0)
# else:
#     crit_zeros = np.array(np.ones (len(y)), dtype = bool)
if args.no_zeros:
    crit_zeros = (y != 0) & ~np.isnan(data_aux).any(axis = 1)
else:
    crit_zeros = ~np.isnan(data_aux).any(axis = 1)
y = y[crit_zeros]

# Get intrinsic parameters.
intrinsic = []
for intr in args.int_pars:
    intrinsic.append (tbdata.field (intr)[crit_zeros])
intrinsic = np.array(intrinsic).transpose()

# Get observable parameters.
observables = []
for obs in fields_beta:
    observables.append (tbdata.field (obs)[crit_zeros])
observables = np.array(observables).transpose()

if args.labels:
    labels = np.array(args.labels, dtype = int)
else:
    labels = np.unique(y)
print "labels = ", labels
tot_lab = 0
for i in range (len(labels)):
    print "class ", labels[i], " = ", (y == labels[i]).sum()
    tot_lab += (y == labels[i]).sum()
print "Labeled total = ", tot_lab
dataAn = BiasAnalyser ()

# Ls, Ns = dataAn.getRandomL (intrinsic, observables, y, labels, N_iter, bins_int,
#                             bins_obs, minElementsBin, N_objs, 
#                             [True, True, False])
increasing_bias = [True, False]
print N_objs, l2_bins_int, bins_obs
print intrinsic.shape, observables.shape, y.shape, labels, increasing_bias, N_calc, l2_bins_int, bins_obs, N_objs, l2_bins_int, bins_obs

Ls, Ns = dataAn.getRandomL (intrinsic, observables, y, labels, 
                            increasing_bias, N_calc, l2_bins_int,
                            bins_obs, minElementsBin = N_objs, 
                            N_objs = N_objs*2**l2_bins_int*bins_obs,
                            bootstrap = True, 
                            kd_tree = "highest_fraction_diference")

print "-------"
print "L = ", Ls.mean(), " +- ", Ls.std()
