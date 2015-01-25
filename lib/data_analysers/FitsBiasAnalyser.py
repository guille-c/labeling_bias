
import sys
import pylab as pl
import pyfits as pf
import numpy as np
import matplotlib as mpl
from scipy.interpolate import *
from matplotlib import cm, colors

from BiasAnalyser import *

class FitsBiasAnalyser (BiasAnalyser):

    def __init__ (self, filename, N = -1):
        self.tbdata = pf.open(filename)[1].data
        if N > 0 and N < len(self.tbdata):
            i_s = np.arange (len(self.tbdata))
            np.random.shuffle(i_s)
            self.tbdata = self.tbdata[i_s[:N]]

    def getFractions (self, x_name, classes_name, nx, ini_crit = True,
                      stars = True, get_is = False, get_N = False, equalBins = False):
        #print self.tbdata.field(x_name).shape
        x = self.tbdata.field(x_name)[ini_crit]
        #print x.shape
        classes = self.tbdata.field(classes_name)[ini_crit]

        crit = (classes < 3)
        if not stars:
            crit = crit & (classes > 0)
        x = x[crit]
        classes = classes[crit]
        N = len(x)
        if N == 0:
            print "GZDataAnalyzer.getFractions N = 0"
            ret = np.array([]), np.array([]), np.array([]), np.array([])
            if get_is:
                ret = ret + (np.array([]),)
            if get_N:
                ret += (N,)
            return ret
        dn = 1. * N / nx
        ns = np.array(np.round(np.arange(0, N, dn)), dtype = "int")
        i_s = x.argsort()
        fs = []
        Es = []
        Ss = []
        if stars:
            stars = []
        i_ret = []
        dx = (x[-1]-x[0]) / nx
        for i in range(len(ns) - 1):
            if equalBins:
                print i, len(equalBins)
                x1 = equalBins[i]
                x2 = equalBins[i+1]
                i1 = i_s [(x[i_s] >= x1) & (x[i_s] < x2)]
                Ni = len(i1)
            else:
                i1 = i_s [ns[i]:ns[i+1]]
                Ni = ns[i+1] - ns[i]
            fs.append(x[i1].mean())
            if stars:
                stars.append(1.*(classes[i1] == 0).sum()/Ni)
            Es.append(1.*(classes[i1] == 1).sum()/Ni)
            Ss.append(1.*(classes[i1] == 2).sum()/Ni)
            i_ret.append(i1)
        i1 = i_s [ns[-1]:]
        Ni = N - ns[-1]
        fs.append(x[i1].mean())
        if stars:
            stars.append(1.*(classes[i1] == 0).sum()/Ni)
        Es.append(1.*(classes[i1] == 1).sum()/Ni)
        Ss.append(1.*(classes[i1] == 2).sum()/Ni)
        i_ret.append(i1)
        if stars:
            ret = np.array(fs), np.array(stars), np.array(Es), np.array(Ss)
        else:
            ret = np.array(fs), np.array(Es), np.array(Ss)
        if get_is:
            ret = ret + (np.array(i_ret),)
        if get_N:
            ret += (N,)
        return ret
    
    def get_rs (self, x_name, classes_name, nx, ini_crit = True,
                get_is = False, get_N = False, equalBins = False, log = False):
        x = self.tbdata.field(x_name)
        classes = self.tbdata.field(classes_name)

        crit = (classes < 3)
        crit = crit & (classes > 0)
        x = x[crit]
        classes = classes[crit]
        N = len(x)
        if N == 0:
            ret = np.array([]), np.array([]), np.array([]), np.array([])
            if get_is:
                ret = ret + (np.array([]),)
            if get_N:
                ret += (N,)
            return ret
        dn = 1. * N / nx
        ns = np.array(np.round(np.arange(0, N, dn)), dtype = "int")
        i_s = x.argsort()
        fs = []
        rs = []
        i_ret = []
        dx = (x[-1]-x[0]) / nx
        for i in range(len(ns) - 1):
            if equalBins:
                print i, len(equalBins)
                x1 = equalBins[i]
                x2 = equalBins[i+1]
                i1 = i_s [(x[i_s] >= x1) & (x[i_s] < x2)]
                Ni = len(i1)
            else:
                i1 = i_s [ns[i]:ns[i+1]]
                Ni = ns[i+1] - ns[i]
            fs.append(x[i1].mean())
            rs.append(1.*(classes[i1] == 1).sum()/(classes[i1] == 2).sum())
            i_ret.append(i1)
        i1 = i_s [ns[-1]:]
        Ni = N - ns[-1]
        fs.append(x[i1].mean())
        rs.append(1.*(classes[i1] == 1).sum()/(classes[i1] == 2).sum())
        i_ret.append(i1)
        if log:
            ret = np.array(fs), np.log(np.array(rs))
        else:
            ret = np.array(fs), np.array(rs)
        if get_is:
            ret = ret + (np.array(i_ret),)
        if get_N:
            ret += (N,)
        return ret
        
    def getPercentiles (self, x, ps, getIndices = False):
        i_s = x.argsort()
        if getIndices:
            return i_s[np.array(np.round(ps*len(i_s)), dtype = "int")]
        else:
            return x[i_s[np.array(np.round(ps*len(i_s)), dtype = "int")]]
           
    def getLsBins (self, field_r, field_c, field_plot, classf, crit = True,
                   bins = (5,5), nx = 20, fractions = False, equalN = True, 
                   minElementsBin = 100):
        tbdata = self.tbdata

        # create array of indices for the data
        if np.isscalar (crit):
            crit = np.arange(len(tbdata.field(field_r)))
        Ls = np.zeros (bins)
        Ns = np.zeros (bins, dtype = int)

        # Check if there is enough data for minElementsBin objects per
        # bin
        if (len (tbdata.field(field_r)[crit]) < bins[0] * bins[1] * nx * minElementsBin):
            print "not enough data", (Ls + 1).sum()
            print " ", len (tbdata.field(field_r)[crit]), " < ", bins[0] * bins[1] * nx * minElementsBin
            return Ls + 1, np.zeros (bins), np.zeros (bins), Ns

        # split tbdata into bins[0] vectors sorted in terms of field_r
        i_r = self.getEqualNumberSplit (tbdata.field(field_r)[crit], bins[0])
    
        # equalN = True: same number of objs. per bin.
        if equalN:
            bins_r_mean = np.zeros (bins)
            bins_c_mean = np.zeros (bins)
        else:
            bins_r_mean = []
            bins_c_mean = []
            #calculate means on bins of field_r 
            for i in range(bins[0]):
                bins_r_mean.append(tbdata.field(field_r)[crit][i_r[i]].mean())
    
        if not equalN:
            iBins = self.getEqualNumberSplit (tbdata.field(field_c), bins[1])
            fieldBins = []
            for i in range(bins[1]):
                fieldBins.append (tbdata.field(field_c)[iBins[i]].min())
                bins_c_mean.append(tbdata.field(field_c)[iBins[i]].mean())
            fieldBins.append(tbdata.field(field_c)[iBins[-1]].max())
    
        for i in range(bins[0]):
            if equalN:
                # get bins in i_r[i] with equal number of objs in field_c
                i_c = self.getEqualNumberSplit (tbdata.field(field_c)[crit][i_r[i]], bins[1])
            else:
                i_c = []
                for j in range (bins[1]):
                    crit1 = ((tbdata.field(field_c)[crit][i_r[i]] >= fieldBins[j]) & 
                            (tbdata.field(field_c)[crit][i_r[i]] < fieldBins[j+1]))
                    i_c.append(np.arange(len(tbdata.field(field_c)[crit][i_r[i]]))[crit1])
    
            for j in range(bins[1]):
                # f21 = np.round(fieldBins[j], 1)
                # f22 = np.round(fieldBins[j+1], 1) Not used?
    
                Ns[i][j] = len(tbdata.field(field_plot)[crit][i_r[i]][i_c[j]])
                # calculate mean field_r and field_c values
                if equalN:
                    bins_r_mean [i][j] = (tbdata.field(field_r)[crit][i_r[i]][i_c[j]]).mean()
                    bins_c_mean [i][j] = (tbdata.field(field_c)[crit][i_r[i]][i_c[j]]).mean()
                if fractions:
                    # if number of objects per bins < minElementsBin, L = -0.00001
                    if len(tbdata.field(field_plot)[crit][i_r[i]][i_c[j]]) < minElementsBin / nx:
                        Ls[i][j] = -0.00001
                        continue
                    crit1 = np.arange(len(tbdata.field(field_plot)))[crit][i_r[i]][i_c[j]]
                    # fractions = E/(E+S)
                    # fs = means of field_plot per bin, rs = fractions of S
                    fs, Es, rs = self.getFractions (field_plot, classf,nx, 
                                                    ini_crit = crit1, stars = False)
                else:
                    # rs = S/E 
                    fs, rs, i_s, N = get_rs (tbdata.field(field_plot)[crit][i_r[i]][i_c[j]], 
                                             tbdata.field(classf)[crit][i_r[i]][i_c[j]], 
                                             nx, get_is = True, get_N = True)
                #mr = rs.mean()
                #Ls [i][j] = np.sqrt(((rs - mr)**2/nx).sum())
                #Ls [i][j] = rs.std()
                if field_plot == "petroRad_r_psf":
                    Ls [i][j] = np.sqrt(((rs - rs[-1])**2).sum()/nx)
                else:
                    Ls [i][j] = np.sqrt(((rs - rs[0])**2).sum()/nx)
        return Ls, bins_r_mean, bins_c_mean, Ns
    
    def getLFractions(self, field_r, field_c, fields_beta, classf, 
                      field_cut = False, cut = 0, bins = (5,5), nx = 20, 
                      minElementsBin = 10, N_bins = -1):
        tbdata = self.tbdata
        if not field_cut:
            crit = True
        elif field_cut[0] == '-':
            crit = tbdata.field(field_cut[1:]) < cut
        else:
            crit = tbdata.field(field_cut) >= cut

        crit = np.arange(len(tbdata.field(fields_beta[0])))[crit]
        N = N_bins*bins[0]*bins[1]*nx
        if N_bins > 0 and N < len(crit):
            i_s = np.arange (len(crit))
            np.random.shuffle(i_s)
            crit = crit[i_s[:N]]

        Npl = len(fields_beta)
        vmax = 0
        L = 0
        for k in range (Npl):
            field_plot = fields_beta[k]
            Ls, rs, cs, Ns = self.getLsBins(field_r, field_c, field_plot, classf, 
                                            bins = bins, minElementsBin = minElementsBin,
                                            nx = nx, fractions = True, crit = crit)
            if (Ls.sum() == bins[0]*bins[1]):
                return 0.
            L += (Ls**2).sum()
            if Ls.max() > vmax:
                vmax = Ls.max()
        return np.sqrt(L/(bins[0] * bins[1] * Npl))

    def getCutLFractions(self, field_r, field_c, fields_beta, classf, 
                         fields_cut, cuts, bins = (5,5), nx = 20, 
                         minElementsBin = 10, N_bins = -1):
        tbdata = self.tbdata

        crit = True
        crit = np.ones(len(tbdata.field(field_r)), dtype = "bool")
        crit = ((tbdata.field(classf) == 1) | (tbdata.field(classf) == 2))
        for i in range (len (fields_cut)):
            if fields_cut[i][0] == '-':
                crit = crit & (tbdata.field(fields_cut[i][1:]) < cuts[i])
            else:
                crit = crit & (tbdata.field(fields_cut[i]) >= cuts[i])
                
        #printing for debugging
        # print crit
        # for i in range (len (fields_cut)):
        #     if fields_cut[i][0] == '-':
        #         print fields_cut[i], cuts[i], tbdata.field(fields_cut[i][1:])[crit].min(), tbdata.field(fields_cut[i][1:])[crit].max()
        #     else:
        #         print fields_cut[i], cuts[i], tbdata.field(fields_cut[i])[crit].min(), tbdata.field(fields_cut[i])[crit].max()
        # print "-----"

        crit = np.arange(len(tbdata.field(fields_beta[0])))[crit]
        # print crit

        N = N_bins*bins[0]*bins[1]*nx
        if N_bins > 0 and N < len(crit):
            i_s = np.arange (len(crit))
            np.random.shuffle(i_s)
            crit = crit[i_s[:N]] #this is wrong
            #crit.sort()
            #crit = crit[:N]
        #else:
            # i_s = np.arange (len(crit))
            # np.random.shuffle(i_s)
            #np.random.shuffle(crit)
            #crit.sort()
            

        N = len(crit)
        print "N cut = ", N

        Npl = len(fields_beta)
        vmax = 0
        L = 0
        for k in range (Npl):
            field_plot = fields_beta[k]
            # print k
            Ls, rs, cs, Ns = self.getLsBins(field_r, field_c, field_plot, classf, 
                                            bins = bins, minElementsBin = minElementsBin,
                                            nx = nx, fractions = True, crit = crit)
            if (Ls.sum() == bins[0]*bins[1]):
                return 0., 0
            L += (Ls**2).sum()
            if Ls.max() > vmax:
                vmax = Ls.max()
        return np.sqrt(L/(bins[0] * bins[1] * Npl)), N

    def getCutLFractionsMean(self, field_r, field_c, fields_beta, classf, 
                         fields_cut, cuts, bins = (5,5), nx = 20, 
                         minElementsBin = 10, N_bins = -1, N_calc = 1):
        Ls = np.zeros(N_calc)
        Ns = np.zeros(N_calc) 
        
        for i in range(N_calc):
            print i
            Ls[i], Ns[i] = self.getCutLFractions(field_r, field_c, fields_beta, classf, 
                                                 fields_cut, cuts, bins, nx, 
                                                 minElementsBin, N_bins)
           
        return Ls, Ns

    def L (self, field_r, field_c, fields_beta, classf, 
           bins = (5,5), nx = 20, minElementsBin = 10, N = -1):
        tbdata = self.tbdata

        crit = ((tbdata.field(classf) == 1) | (tbdata.field(classf) == 2))
        crit = np.arange(len(tbdata.field(fields_beta[0])))[crit]
        # crit are indices of class == 1 and class == 2 objects

        #N = N_bins*bins[0]*bins[1]*nx
        if N > 0 and N < len(crit):
            i_s = np.arange (len(crit))
            np.random.shuffle(i_s)
            crit = crit[i_s[:N]] # N shuffled indices of class == 1
                                 # and class == 2 objects

        N = len(crit)

        Npl = len(fields_beta)
        L = 0
        for k in range (Npl):
            field_plot = fields_beta[k]
            Ls, rs, cs, Ns = self.getLsBins(field_r, field_c, field_plot, 
                                            classf, bins = bins, 
                                            minElementsBin = minElementsBin,
                                            nx = nx, fractions = True, crit = crit)
            if (Ls.sum() == bins[0]*bins[1]):
                return 0., 0
            L += (Ls**2).sum()
        return np.sqrt(L/(bins[0] * bins[1] * Npl)), N
