import sys
import pylab as pl
import pyfits as pf
import numpy as np
import matplotlib as mpl
from scipy.interpolate import *
from matplotlib import cm, colors

class BiasAnalyser ():

    def __init__ (self, N = -1):
        self.N = N

    def getFractions (self, x, classes, labels, nx, ini_crit = True,
                      get_is = False, get_N = False, equalBins = False):

        # Create criteria for selecting sources that correspond to
        # desired labels.
        if np.isscalar (ini_crit):
            crit = np.where(np.in1d (classes, labels))
        else:
            crit = np.where(np.in1d (classes, labels))[ini_crit]

        x = x[crit]
        classes = classes[crit]

        N = len(x)

        # If no sources present return empy array.
        if N == 0:
            print "BiasAnalyzer.getFractions N = 0"
            ret = np.array([]), np.array([]), np.array([]), np.array([])
            if get_is:
                ret = ret + (np.array([]),)
            if get_N:
                ret += (N,)
            return ret

        # Create array containing equal number of sources per bin.
        ns = np.array(np.round(np.linspace (0, N, nx + 1)), dtype = int)
        i_s = x.argsort()
        x_means = []

        # Create array for fractions per bin.
        fractions = np.zeros((len(ns)-1, len(labels)))
        i_ret = []
        dx = (x[-1]-x[0]) / nx

        for i in range(len(ns) - 2):
            if equalBins:
                print i, len(equalBins)
                x1 = equalBins[i]
                x2 = equalBins[i+1]
                i1 = i_s [(x[i_s] >= x1) & (x[i_s] < x2)]
                Ni = len(i1)
            else:
                # Select objects in bin i
                i1 = i_s [ns[i]:ns[i+1]]
                Ni = ns[i+1] - ns[i]
            x_means.append(x[i1].mean())

            # For each label calculate fractions
            for k in range(len(labels)):
                fractions[i, k] = 1.*(classes[i1] == labels[k]).sum()/Ni
            i_ret.append(i1)
        i1 = i_s [ns[-2]:]
        Ni = N - ns[-2]
        x_means.append(x[i1].mean())
        for k in range(len(labels)):
            fractions[len(ns) - 2, k] = 1.*(classes[i1] == labels[k]).sum()/Ni
        i_ret.append(i1)
        ret = np.array(x_means), fractions
        if get_is:
            ret = ret + (np.array(i_ret),)
        if get_N:
            ret += (N,)
        return ret

    # splits x into nx sorted vectors
    def getEqualNumberSplit (self, x, nx):
        N = len(x)
        dn = 1. * N / nx
        ns = np.array(np.round(np.linspace(0, N, nx+1)), dtype = "int")
        i_s = x.argsort()
        ret = []
        for i in range(len(ns) - 1):
            i1 = i_s [ns[i]:ns[i+1]]
            ret.append (i1)
        ret = np.array(ret)
        return ret

    def KDTree (self, x, pow_n, i_array, i_var = 0, i_ret = 0, 
                crit = "iterative", y = None, labels = None):
        if pow_n == 0:
            return np.zeros(len(i_array)) + i_ret
            
        if crit == "iterative":
            this_i_var = i_var
            sep = np.median(x[:, i_var])
            # print sep
        elif crit == "highest_fraction_difference":
            if y is None or labels is None:
                print "ERROR BiasAnalyser.KDTree: y and labels required for ", crit
                exit()
            f_d = np.zeros(x.shape[1])  # Fraction differences
            seps = np.zeros(x.shape[1]) # Thresholds for each
                                        # intrinsic variable
            
            # calculate splitting fractions
            for i in range(x.shape[1]):
                seps[i] = np.median(x[:, i])
                if (x[:, i] < seps[i]).sum() == 0 or (x[:, i] >= seps[i]).sum() == 0:
                    seps[i] = np.mean(x[:, i])
                f_left = 0
                f_right = 0
                for k in labels:
                    if len(x) != len(y):
                        print len(x[:, i] < seps[i]), len(y), len(x)
                        print (y[x[:, i] < seps[i]] == k).sum()
                        print (x[:, i] < seps[i]).sum()
                        exit()
                    f_left  = (1. * (y[x[:, i] < seps[i]] == k).sum() / 
                               (x[:, i] < seps[i]).sum())
                    f_right = (1. * (y[x[:, i] >= seps[i]] == k).sum() / 
                               (x[:, i] >= seps[i]).sum())
                    f_d[i] += np.abs (f_left - f_right)
                f_d [i] = f_d [i] / len (labels)
            this_i_var = np.argmax (f_d)
            sep = seps [this_i_var]
            # print " " , i_ret, " KDTree f_d = ", f_d, this_i_var, x.shape, f_d.prod()
            # if f_d.prod() == 0.:
            #     i = np.argmin (f_d)
            #     print seps[i]
            #     print i, (y[x[:, i] < seps[i]]).sum(), (y[x[:, i] < seps[i]]).sum()
            #     exit()

        else:
            print "ERROR BiasAnalyser.KDTree: ", crit, " criteria not implemented."
            exit()

        crit_left = np.arange(len(x))[x[:, this_i_var] < sep]
        crit_right = np.arange(len(x))[x[:, this_i_var] >= sep]
        i_var_new = (this_i_var + 1)%x.shape[1]
        
        left = self.KDTree (x[crit_left], pow_n - 1, crit_left, 
                             i_var_new, i_ret, crit, y[crit_left], labels)
        right = self.KDTree (x[crit_right], pow_n - 1, crit_right, 
                             i_var_new, i_ret + 2**(pow_n-1), crit, y[crit_right], labels)
        ret = np.zeros(x.shape[0], dtype = "int")
        ret[crit_left] = left
        ret[crit_right] = right
        return ret
        
    def calculateSigma2 (self, intrinsic, observable, y, labels, 
                         log2_bins_int, bins_obs, increasing_bias, 
                         kd_tree = "iterative"):
        kd_tree = self.KDTree (intrinsic, log2_bins_int, 
                               np.arange(intrinsic.shape[0]), crit = kd_tree,
                               y = y, labels = labels)
        kd_keys = np.unique(kd_tree)
        sigma2 = np.zeros((observable.shape[1], len(labels), 2**log2_bins_int))
        
        for j in range (observable.shape[1]):
            for q in kd_keys:
                i_bin_int = (kd_tree == q)
                # rs.shape = (bins_obs, len(labels))
                fs, rs = self.getFractions (observable[:, j][i_bin_int], 
                                            y[i_bin_int], labels, bins_obs)
                for k in range (len(labels)):
                    if increasing_bias[j]:
                        sigma2 [j, k, q] = ((rs[:, k] - rs[0, k])**2).sum()/bins_obs
                    else:
                        sigma2 [j, k, q] = ((rs[:, k] - rs[-1, k])**2).sum()/bins_obs
        return sigma2
    
    def getFractionsPerObject (self, intrinsic, observable, y, labels, 
                               log2_bins_int, bins_obs, increasing_bias, 
                               kd_tree = "iterative"):
        kd_tree = self.KDTree (intrinsic, log2_bins_int, 
                               np.arange(intrinsic.shape[0]), crit = kd_tree,
                               y = y, labels = labels)
        kd_keys = np.unique(kd_tree)
        sigma2 = np.zeros((observable.shape[1], len(labels), 2**log2_bins_int))
        output_size = np.concatenate((observable.shape, [len(labels)]))
        print "output_size = ", output_size
        int_frac = np.zeros(output_size)
        obs_frac = np.zeros(output_size)

        for j in range (observable.shape[1]):
            for q in kd_keys:
                i_bin_int = (kd_tree == q)
                is_bin = np.arange(observable.shape[0])[i_bin_int]
                # rs.shape = (bins_obs, len(labels))
                fs, rs, i_s = self.getFractions (observable[:, j][i_bin_int], 
                                            y[i_bin_int], labels, bins_obs,
                                            get_is = True)
                for k in range (len(labels)):
                    for i in range(len(i_s)):
                        print is_bin[i_s[i]].min(), is_bin[i_s[i]].max()
                        obs_frac[is_bin[i_s[i]], j, k] = rs[i, k]
                        if increasing_bias[j]:
                            int_frac[is_bin[i_s[i]], j, k] = rs[0, k]
                            #sigma2 [j, k, q] = ((rs[:, k] - rs[0, k])**2).sum()/bins_obs
                        else:
                            int_frac[is_bin[i_s[i]], j, k] = rs[-1, k]
                            #sigma2 [j, k, q] = ((rs[:, k] - rs[-1, k])**2).sum()/bins_obs
        return int_frac, obs_frac
        

    def getLsBins (self, intrinsic, observable, y, labels, crit = True,
                   bins_in = (5,5), bins_ob = 20, equalN = True, 
                   minElementsBin = 100, increasing_bias = True):

        # create array of indices for the data
        if np.isscalar (crit):
            crit = np.arange(intrinsic.shape[0])

        Ls = np.zeros (bins_in)
        Ns = np.zeros (bins_in, dtype = int)

        # Check if there is enough data for minElementsBin objects per
        # bin
        N_data = intrinsic[crit].shape[0]
        N_objs_per_bin = bins_in.prod() * bins_ob * minElementsBin
        if (N_data < N_objs_per_bin):
            print "BiasAnalyser.getLsBins: not enough data", (Ls + 1).sum()
            print " ", N_data, " < ", N_objs_per_bin
            return Ls + 1, np.zeros (bins_in), np.zeros (bins_in), Ns

        # split tbdata into bins[0] vectors sorted in terms of field_r
        i_r = self.getEqualNumberSplit (intrinsic[:, 0][crit], bins_in[0])

        # equalN = True: same number of objs. per bin.
        if equalN:
            bins_r_mean = np.zeros (bins_in)
            bins_c_mean = np.zeros (bins_in)
        else:
            bins_r_mean = []
            bins_c_mean = []
            #calculate means on bins of field_r 
            for i in range(bins_in[0]):
                bins_r_mean.append(intrinsic[:, 0][crit][i_r[i]].mean())
    
        if not equalN:
            iBins = self.getEqualNumberSplit (intrinsic[:, 1], bins_in[1])
            fieldBins = []
            for i in range(bins[1]):
                fieldBins.append (intrinsic[:, 1][iBins[i]].min())
                bins_c_mean.append(intrinsic[:, 1][iBins[i]].mean())
            fieldBins.append(intrinsic[:, 1][iBins[-1]].max())
    
        for i in range(bins_in[0]):
            if equalN:
                # get bins in i_r[i] with equal number of objs in field_c
                i_c = self.getEqualNumberSplit (intrinsic[:, 1][crit][i_r[i]], 
                                                bins_in[1])
            else:
                i_c = []
                for j in range (bins[1]):
                    crit1 = ((intrinsic[:, 1][crit][i_r[i]] >= fieldBins[j]) & 
                             (intrinsic[:, 1][crit][i_r[i]] < fieldBins[j+1]))
                    i_c.append(np.arange(len(intrinsic[:, 1][crit][i_r[i]]))[crit1])
            for j in range(bins_in[1]):
                Ns[i][j] = len(observable[crit][i_r[i]][i_c[j]])

                # calculate mean field_r and field_c values
                if equalN:
                    bins_r_mean [i][j] = (intrinsic[:, 0][crit][i_r[i]][i_c[j]]).mean()
                    bins_c_mean [i][j] = (intrinsic[:, 1][crit][i_r[i]][i_c[j]]).mean()
                # if number of objects per bins < minElementsBin, L = -0.00001
                if len(observable[crit][i_r[i]][i_c[j]]) < minElementsBin / bins_ob:
                    Ls[i][j] = -0.00001
                    continue
                crit1 = np.arange(len(observable))[crit][i_r[i]][i_c[j]]

                # fs = means of field_plot per bin, rs = fractions of S
                fs, rs = self.getFractions (observable, y, labels, bins_ob,
                                            ini_crit = crit1)
                L_ij = 0
                for k in range(len(labels)):
                    if increasing_bias:
                        L_ij += np.sqrt(((rs[:, k] - rs[0, k])**2).sum()/bins_ob)
                    else:
                        L_ij += np.sqrt(((rs[:, k] - rs[-1, k])**2).sum()/bins_ob)
                Ls [i][j] = L_ij / len(labels)
        return Ls, bins_r_mean, bins_c_mean, Ns

    def L (self, intrinsic, observables, y, labels, increasing_bias, 
           log2_bins_int, bins_ob = 20, minElementsBin = 10, N = -1, 
           bootstrap = False, kd_tree = "iterative"):

        # define a criteria for choosing objects with labels in
        # "labels"
        crit = np.where(np.in1d (y, labels))[0]

        if N > len(crit):
            print "N > len(crit)", N, len(crit)
            return 0., 0.

        if N > 0 and N < len(crit):
            i_s = np.arange (len(crit))
            np.random.shuffle(i_s)
            crit = crit[i_s[:N]]     # N shuffled indices of class in
                                     # "labels"

        N = len(crit)
        if bootstrap:
            crit = np.random.choice (crit, size = N, replace = True)

        N_bins = 2**log2_bins_int*bins_ob
        if N < N_bins*minElementsBin:
            print "N < N_bins*minElementsBin", N, N_bins*minElementsBin
            return 0., 0.

        Npl = observables.shape[1]
        L = 0
        sigma2 = self.calculateSigma2 (intrinsic[crit], observables[crit], y[crit], 
                                       labels, log2_bins_int, bins_ob, 
                                       increasing_bias, kd_tree = kd_tree)

        return np.sqrt(sigma2.sum()/ (np.prod(sigma2.shape))), N

    def getRandomL (self, intrinsic, observables, y, labels, increasing_bias,
                    N_calc, bins_in, bins_ob = 20, minElementsBin = 10, 
                    N_objs = -1, bootstrap = False, kd_tree = "iterative"):
        Ls = np.zeros(N_calc)
        Ns = np.zeros(N_calc)
        
        for i in range(N_calc):
            Ls[i], Ns[i] = self.L (intrinsic, observables, y, labels, 
                                   increasing_bias, bins_in, bins_ob, 
                                   minElementsBin, N_objs, bootstrap = bootstrap,
                                   kd_tree = kd_tree)
                                  # (intrinsic, observables, y, labels, 
                                  #  bins_in, bins_ob, minElementsBin, 
                                  #  N_objs, increasing_bias)
        return Ls, Ns
