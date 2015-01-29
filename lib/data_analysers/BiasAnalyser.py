
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

        crit = np.where(np.in1d (classes, labels))[0][ini_crit]

        x = x[crit]
        classes = classes[crit]

        N = len(x)
        if N == 0:
            print "BiasAnalyzer.getFractions N = 0"
            ret = np.array([]), np.array([]), np.array([]), np.array([])
            if get_is:
                ret = ret + (np.array([]),)
            if get_N:
                ret += (N,)
            return ret
        ns = np.array(np.round(np.linspace (0, N, nx + 1)), dtype = int)
        i_s = x.argsort()
        x_means = []
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
                i1 = i_s [ns[i]:ns[i+1]]
                Ni = ns[i+1] - ns[i]
            x_means.append(x[i1].mean())
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
        if ( N_data < N_objs_per_bin):
            print "not enough data", (Ls + 1).sum()
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

    def L (self, intrinsic, observables, y, labels, 
           bins_in = (5,5), bins_ob = 20, minElementsBin = 10, 
           N = -1, increasing_bias = True):

        # define a criteria for choosing objects with labels in
        # "labels"
        crit = np.where(np.in1d (y, labels))[0]

        if N > 0 and N < len(crit):
            i_s = np.arange (len(crit))
            np.random.shuffle(i_s)
            crit = crit[i_s[:N]] # N shuffled indices of class == 1
                                 # and class == 2 objects

        N = len(crit)

        Npl = observables.shape[1]
        L = 0
        print observables.shape
        for k in range (Npl):
            print k
            observable = observables[:, k]
            if np.isscalar(increasing_bias):
                inc = increasing_bias
            else:
                inc = increasing_bias[k]
            Ls, rs, cs, Ns = self.getLsBins(intrinsic, observable, 
                                            y, labels, bins_in = bins_in, 
                                            minElementsBin = minElementsBin,
                                            bins_ob = bins_ob, 
                                            crit = crit, increasing_bias = inc)
            if (Ls.sum() == bins_in.prod()):
                return 0., 0
            L += (Ls**2).sum()
        return np.sqrt(L/(bins_in.prod() * Npl)), N

    def getRandomL (self, intrinsic, observables, y, labels, N_calc,
                    bins_in = (5,5), bins_ob = 20, minElementsBin = 10, 
                    N = -1, increasing_bias = True):
        Ls = np.zeros(N_calc)
        Ns = np.zeros(N_calc) 
        
        for i in range(N_calc):
            Ls[i], Ns[i] = self.L (intrinsic, observables, y, labels, 
                                   bins_in, bins_ob, minElementsBin, 
                                   N, increasing_bias)
        return Ls, Ns
