import sys
import os
#sys.path.append(sys.argv[1])
import numpy as np

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

field_r = "petroRad_r_kpc"
field_c = "absPetroMag_r"
fields_beta = ["z", "corrMag_r", "petroRad_r_psf"]
classf = "class_SVM"
classf = sys.argv[2]
N = -1843
N = 50000

dataAn = BiasAnalyser (sys.argv[1])

print "E = ", (dataAn.tbdata.field(classf) == 1).sum()
print "S = ", (dataAn.tbdata.field(classf) == 2).sum()
print "Total = ", (dataAn.tbdata.field(classf) == 1).sum() + (dataAn.tbdata.field(classf) == 2).sum()

bins = (10, 10)
nx = 10

N_iter = 10

fields_cut = []
cuts = []
N_calc = 5
minElementsBin = 5
saveLs (dataAn, field_r, field_c, fields_beta, classf, fields_cut, cuts, np.array ([[10, 10, 10]]), [50],
        minElementsBin, N_calc, "./", txt = "")

Ls = np.zeros (N_iter)
for i in range (N_iter):
    L, N = dataAn.L (field_r, field_c, fields_beta, classf, bins = bins, nx = nx, N = N)
    print i, "------"
    print "L = ", L
    print "N per bin = ", 1. *N / (nx * bins[0] * bins[1])
    print "N = ", N
    Ls[i] = L
print "-------"
print "L = ", Ls.mean(), " +- ", Ls.std()
