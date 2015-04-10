import sys
import os
#sys.path.append(sys.argv[1])
import numpy as np
import argparse
from data_analysers.BiasAnalyser import *
import bias_methods as bm   

def calculateDL (dataAn, bins_obs, intrinsic, observables, 
                 y, labels, increasing_bias, log2_bins_int, N_calc):
    N_bins = 2**log2_bins_int*bins_obs
    L1, N1 = dataAn.getRandomL (intrinsic, observables, y, labels,
                                increasing_bias, N_calc, 
                                log2_bins_int, bins_obs,
                                minElementsBin = 0, 
                                N_objs = N_labeled, bootstrap = True)
    L2, N2 = dataAn.getRandomL (intrinsic, observables, y, labels, 
                                increasing_bias, N_calc, 
                                log2_bins_int, bins_obs,
                                minElementsBin = 0, 
                                N_objs = N_labeled - 10*N_bins, bootstrap = True)
    if L1.prod() == 0 or L2.prod() == 0:
        return np.inf, 0
    dL = (L2 - L1)
    dLm, dLs = L2.mean() - L1.mean(), dL.std()
    return dLm, dLs


def searchBestObsBins (dataAn, log_bins_ini, log_bins_end, 
                       dL_ini, dL_end, intrinsic, 
                       observables, y, labels, increasing_bias,
                       N_calc, tol):
    if (log_bins_end - log_bins_ini) <= 1:
        print "returning"
        return [], [], []
    current_bins = int (np.trunc((log_bins_end + log_bins_ini)*0.5))
    
    N_labeled = len(np.where(np.in1d (y, labels))[0])
    N_bins = current_bins**2
    obj_per_bin = int(np.trunc(1.* intrinsic.shape[0]/N_bins))
    if obj_per_bin < 10:
        dLm = np.inf
        dLs = 0
    else:
        # L1, N1 = dataAn.getRandomL (intrinsic, observables, y, labels,
        #                             increasing_bias, N_calc, 
        #                             log2_bins_int, current_bins,
        #                             minElementsBin = 0, 
        #                             N_objs = N_labeled, bootstrap = True)
        # L2, N2 = dataAn.getRandomL (intrinsic, observables, y, labels, 
        #                             increasing_bias, N_calc, 
        #                             log2_bins_int, current_bins,
        #                             minElementsBin = 0, 
        #                             N_objs = N_labeled - 10*N_bins, bootstrap = True)
        # dL = (L2 - L1)
        # dLm, dLs = L2.mean() - L1.mean(), dL.std()

        dLm, dLs = calculateDL (dataAn, 2**current_bins, intrinsic, observables, 
                                y, labels, increasing_bias, current_bins, N_calc)

        print "dL[", current_bins, "] = ", dLm, dLs, log_bins_ini, log_bins_end, obj_per_bin
        #print L1

    if dLm >= tol:
        ret_dLm, ret_dLs, ret_bins = searchBestObsBins (dataAn, log_bins_ini, 
                                                        current_bins, 
                                                        dL_ini, dLm, intrinsic, 
                                                        observables, 
                                                        y, labels, 
                                                        increasing_bias, N_calc, tol)
    else:
        ret_dLm, ret_dLs, ret_bins = searchBestObsBins (dataAn, current_bins, 
                                                        log_bins_end,
                                                        dLm, dL_end, 
                                                        intrinsic, observables, 
                                                        y, labels, 
                                                        increasing_bias, N_calc, tol)
    return [dLm] + ret_dLm, [dLs] + ret_dLs, [current_bins] + ret_bins
    
parser = argparse.ArgumentParser(description='Obtain least biased obs bin.')

parser.add_argument("table_file",
                    help = "table with variables and labels.")
parser.add_argument("label_name", nargs = "+",
                    help = "name of the column containing the labels.")
parser.add_argument("--labels", nargs = "+",
                    help = "Labels to be used.")
parser.add_argument("--int_pars", 
                    default = ["petroRad_r_kpc","absPetroMag_r", "z"], 
                    nargs = "+",
                    help = "Name of the columns used as intrinsic parameters.")
parser.add_argument("--obs_pars", default = ["corrMag_r", "petroRad_r_psf"],
                    help = "Name of the columns used as observable parameters.")
parser.add_argument("--N_iter", default = 5, type = int,
                    help = "Number of calculations of L to calculate means and standard deviations.")
parser.add_argument("--out_path", default = "./",
                    help = "Path where to leave results.")
parser.add_argument("--pbb_thresholds", nargs = "+",
                    help = "Threshold on the probabilities")
parser.add_argument("--no_zeros", action='store_const', const = True,
                    help = "Do not consider labels that don't match the pbb. thresholds")
parser.add_argument("--tol", default = 1E-2, type = float,
                    help = "Tolerance for accepting a bins")

args = parser.parse_args()

fields_int = args.int_pars
fields_obs = args.obs_pars
classf = args.label_name
if args.pbb_thresholds:
    pbb_thresholds = np.array (args.pbb_thresholds).astype(np.float)
else:
    pbb_thresholds = args.pbb_thresholds
N_iter = args.N_iter
#bins_obs_range = np.array(args.bins_obs_range, dtype = int)
out_path = args.out_path
tol = args.tol

dataAn = BiasAnalyser ()

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

crit_zeros = np.ones (len(y), dtype = bool)
if args.no_zeros:
    crit_zeros = (y != 0)
y = y[crit_zeros]

print "Number of objects = ", y.shape

# Get intrinsic parameters.
# intrinsic = np.array([tbdata.field (field_r)[crit_zeros],
#                       tbdata.field (field_c)[crit_zeros]]).transpose()
intrinsic = []
for intr in fields_int:
    intrinsic.append (tbdata.field (intr)[crit_zeros])
intrinsic = np.array(intrinsic).transpose()

# Get observable parameters.
observables = []
for obs in fields_obs:
    observables.append (tbdata.field (obs)[crit_zeros])
observables = np.array(observables).transpose()

if args.labels:
    labels = np.array(args.labels, dtype = int)
else:
    labels = np.unique(y)

N_labeled = len(np.where(np.in1d (y, labels))[0])
increasing_bias = [True, False]

# # Calculate initial dL (for min and max). max is the maximum number of
# # bins we can use to get at least two objects per bin.
# bins_min = 2
# obj_per_bin = int(np.trunc(intrinsic.shape[0]/(2.**log2_bins_int*bins_min)))
# Lmin1, N1 = dataAn.L (intrinsic, observables, y, labels,
#                       increasing_bias, log2_bins_int, bins_min,
#                       minElementsBin = obj_per_bin, N = N_labeled)
# Lmin2, N2 = dataAn.L (intrinsic, observables, y, labels,
#                       increasing_bias, log2_bins_int, bins_min,
#                       minElementsBin = obj_per_bin, 
#                       N = N_labeled - 1)
# print "L_min = ", Lmin1, Lmin2
# dL_min = (Lmin1 - Lmin2)
# 
# bins_max = int (intrinsic.shape[0]/(2.**log2_bins_int)/2)
# obj_per_bin = 0#int(np.trunc(intrinsic.shape[0]/(2.**log2_bins_int*bins_max)))
# print bins_max, obj_per_bin
# Lmax1, N1 = dataAn.L (intrinsic, observables, y, labels,
#                       increasing_bias, log2_bins_int, bins_max,
#                       minElementsBin = obj_per_bin,
#                       N_objs = N_labeled, bootstrap = True)
# Lmax2, N2 = dataAn.L (intrinsic, observables, y, labels,
#                       increasing_bias, log2_bins_int, bins_max,
#                       minElementsBin = obj_per_bin, 
#                       N_objs = N_labeled - 10*N_bins, bootstrap = True)
# print "L_max = ", Lmax1, Lmax2
# dL_max = (Lmax1 - Lmax2)

log_bins_min = 1
dL_min, dLs_min = calculateDL (dataAn, 2**log_bins_min, intrinsic, observables, 
                               y, labels, increasing_bias, log_bins_min, N_iter)
N_obj_per_bin = 15
bins_max = np.sqrt(1.*intrinsic.shape[0]/N_obj_per_bin)
log_bins_max = int ( np.trunc(np.log2 (bins_max)))
dL_max, dLs_max = calculateDL (dataAn, 2**log_bins_max, intrinsic, observables, 
                               y, labels, increasing_bias, log_bins_max, N_iter)

dLm, dLs, bins_obs = searchBestObsBins (dataAn, log_bins_min, log_bins_max, 
                                        dL_min, dL_max, 
                                        intrinsic, observables, y, labels, 
                                        increasing_bias, N_iter, tol)

dLm = [dL_min] + dLm + [dL_max]
bins_obs = [log_bins_min] + bins_obs + [log_bins_max]
dLm = np.array(dLm)
bins_obs = np.array(bins_obs)
print dLm
print bins_obs

print "best = ", np.max(bins_obs[(dLm < tol) & (dLm != 0)])
print bins_obs[(dLm < tol)], (dLm < tol)

pl.clf()
#pl.errorbar (bins_obs, dLm, yerr = dLs)
pl.axhline (tol)
pl.plot(bins_obs, dLm, "x")
pl.xscale("log")
pl.xlim([np.min(bins_obs)*0.5, np.max(bins_obs)*1.5, ])
pl.savefig ("leastBiasedObsBins.eps")


