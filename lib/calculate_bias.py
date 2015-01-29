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
parser.add_argument("--number_objects", metavar = "N", type = int, default = 50000, 
                    help = "number of objects to be used for calculating the bias.")
parser.add_argument("--int_pars", default = ["petroRad_r_kpc","absPetroMag_r"], 
                    nargs = 2,
                    help = "Name of the columns used as intrinsic parameters.")
parser.add_argument("--obs_pars", default = ["z", "corrMag_r", "petroRad_r_psf"],
                    help = "Name of the columns used as observable parameters.")
parser.add_argument("--pbb_thresholds", nargs = "+",
                    help = "Threshold on the probabilities")
parser.add_argument("--no_zeros", action='store_const', const = True,
                    help = "Do not consider labels that don't match the pbb. thresholds")

args = parser.parse_args()

# Number of objects for calculating bias.
N_objs = args.number_objects

# Intrinsic parameters. Currently accept exactly two.
field_r = args.int_pars[0]
field_c = args.int_pars[1]

fields_beta = args.obs_pars   # Observable parameters. As many as
                              # needed.
classf = args.label_name      # Label or probability column name.
if args.pbb_thresholds:
    pbb_thresholds = np.array (args.pbb_thresholds).astype(np.float)
else:
    pbb_thresholds = args.pbb_thresholds

dataAn = BiasAnalyser ()

bins_int = np.array([10, 10]) # Number of bins for intrinsic pars.
bins_obs = 10                 # Number of bins for observable pars.
N_iter = 10                   # Number of iterations for mean and std.
minElementsBin = 5            # Min. number of objs. per bin.
np.random.seed (0)

tbdata = pf.open(args.table_file)[1].data # Open fits table file.
N_tot = len(tbdata)

# i_s is used to randomize the data.
i_s = np.arange (N_tot)
np.random.shuffle(i_s)
tbdata = tbdata[i_s]

y = bm.createLabels (tbdata, classf, pbb_thresholds)

crit_zeros = np.ones (len(y))
if args.no_zeros:
    crit_zeros = (y != 0)
y = y[crit_zeros]

# Get intrinsic paraeters.
intrinsic = np.array([tbdata.field (field_r)[crit_zeros],
                      tbdata.field (field_c)[crit_zeros]]).transpose()

# Get observable parameters.
observables = []
for obs in fields_beta:
    observables.append (tbdata.field (obs)[crit_zeros])
observables = np.array(observables).transpose()

labels = np.unique(y)
print "labels = ", labels
print (y == labels[0]).sum()
print (y == labels[1]).sum()

dataAn = BiasAnalyser ()

Ls, Ns = dataAn.getRandomL (intrinsic, observables, y, labels, N_iter, bins_int,
                            bins_obs, minElementsBin, N_objs, 
                            [True, True, False])

print "-------"
print "L = ", Ls.mean(), " +- ", Ls.std()
