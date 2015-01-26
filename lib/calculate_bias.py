import sys
import os
#sys.path.append(sys.argv[1])
import numpy as np
import argparse

def saveLs (dataAn, field_r, field_c, fields_beta, classf, fields_cut, cuts, bins, Ns,
            minElementsBin, N_calc, out_path, txt = ""):
    for i in np.arange (len(bins)):
        for j in np.arange (len(Ns)):
            print Ns[j]
            filenameLN = out_path + "L_N" + txt + str(bins[i]) + "_" + str(Ns[j]) + ".dat"
            if os.path.exists(filenameLN):
                continue
            L1, N1 = dataAn.getCutLFractionsMean (field_r, field_c, fields_beta, classf, 
                                                  fields_cut, cuts, 
                                                  bins = bins[i][:2], nx = bins[i][2], 
                                                  minElementsBin = minElementsBin, 
                                                  N_bins = Ns[j], N_calc = N_calc)
            print filenameLN
            file = open (filenameLN, "w")
            file.write ("L L_std N\n" + str(L1.mean()) + " " + str(L1.std()) + 
                        " " + str(N1.mean()) + "\n")
            file.close ()

from data_analysers.BiasAnalyser import *
from data_analysers.FitsBiasAnalyser import *

parser = argparse.ArgumentParser(description='Calculate labeling bias of a data-set.')

parser.add_argument("table_file",
                    help = "table with variables and labels.")
parser.add_argument("label_name", help = "name of the column containing the labels.")
parser.add_argument("--number_objects", metavar = "N", type = int, default = 50000, 
                    help = "number of objects to be used for calculating the bias.")
parser.add_argument("--int_pars", default = ["petroRad_r_kpc","absPetroMag_r"], 
                    nargs = 2,
                    help = "Name of the columns used as intrinsic parameters.")
parser.add_argument("--obs_pars", default = ["z", "corrMag_r", "petroRad_r_psf"],
                    help = "Name of the columns used as observable parameters.")
parser.add_argument("--pbb_threshold", 
                    help = "Threshold on the probabilities")

args = parser.parse_args()

# Number of objects for calculating bias.
N = args.number_objects

# Intrinsic parameters. Currently accept exactly two.
field_r = args.int_pars[0]
field_c = args.int_pars[1]

fields_beta = args.obs_pars   # Observable parameters. As many as
                              # needed.
classf = args.label_name      # Label or probability column name.

dataAn = BiasAnalyser ()

bins_int = np.array([10, 10]) # Number of bins for intrinsic pars.
bins_obs = 10                 # Number of bins for observable pars.
N_iter = 10                   # Number of iterations for mean and std.
minElementsBin = 5            # Min. number of objs. per bin.
np.random.seed (0)

tbdata = pf.open(args.table_file)[1].data # Open fits table file.

# i_s is used to randomize the data.
i_s = np.arange (len(tbdata))
np.random.shuffle(i_s)
tbdata = tbdata[i_s[:N]]

# Get intrinsic paraeters.
intrinsic = np.array([tbdata.field (field_r),
                      tbdata.field (field_c)]).transpose()

# Get observable parameters.
observables = []
for obs in fields_beta:
    observables.append (tbdata.field (obs))
observables = np.array(observables).transpose()

# Get labels / probabilities column.
y = tbdata.field (classf)
labels = np.unique(y)

dataAn = BiasAnalyser ()

Ls = np.zeros (N_iter)
for i in range (N_iter):
    # Calculate L for iteration i.
    L, N = dataAn.L (intrinsic, observables, y, labels, bins_int,
                     bins_obs, minElementsBin, N, [True, True, False])
    print i, "------"
    print "L = ", L
    print "N per bin = ", 1. *N / (bins_obs * bins_int.prod())
    print "N = ", N
    Ls[i] = L
    exit()

print "-------"
print "L = ", Ls.mean(), " +- ", Ls.std()
