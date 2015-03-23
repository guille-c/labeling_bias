import numpy as np
import os

def saveLs (dataAn, intrinsic, observables, y, labels, bins,
            N_objs, N_calc, increasing_bias = True, out_path = "./", txt = ""):
    print bins, N_objs
    Ls = np.zeros((len(N_objs), bins.shape[0]))
    Ls_std = np.zeros((len(N_objs), bins.shape[0]))
    for i in np.arange (len(bins)):
        for j in np.arange (len(N_objs)):
            print N_objs[j]
            filenameLN = (out_path + "L_N" + txt + str(bins[i]) + "_" + 
                          str(N_objs[j]) + "_" + str(N_calc) + ".dat")
            if os.path.exists(filenameLN):
                # read from file
                file = open (filenameLN, "r")
                file.readline()
                L1, L1_s = file.readline().split()[:2]
                print bins[i], N_objs[j], L1
                Ls[j, i], Ls_std[j, i] = float(L1), float(L1_s)
            else:
                print "Calculating L for ", N_objs[j], bins[i]
                L1, N1 = dataAn.getRandomL (intrinsic, observables, y, labels, 
                                            increasing_bias, N_calc, bins[i][0],
                                            bins[i][1], minElementsBin = N_objs[j], 
                                            N_objs = N_objs[j]*2**bins[i][0]*bins[i][1],
                                            bootstrap = True, 
                                            kd_tree = "highest_fraction_diference")

                print filenameLN, L1
                file = open (filenameLN, "w")
                file.write ("L L_std N\n" + str(L1.mean()) + " " + str(L1.std()) + 
                            " " + str(N1.mean()) + "\n")
                file.close ()
                Ls[j, i], Ls_std[j, i] = L1.mean(), L1.std()
    return Ls, Ls_std

# Get labels / probabilities column.
def createLabels (tbdata, classf, pbb_thresholds):
    N_tot = len(tbdata)
    if not np.isscalar(pbb_thresholds):
        if not (pbb_thresholds is None) and len (classf) != len (pbb_thresholds):
            print "ERROR: the number of pbb fields does not match the number of thresholds."
            exit()
        y = np.zeros(N_tot)
        for i in range(len(classf)):
            print classf[i], classf, tbdata.field (classf[i]).shape
            crit = (tbdata.field (classf[i]) > pbb_thresholds[i])
            print i, crit.shape
            y[crit] = i + 1
    else:
        if len(classf) != 1:
            #print "ERROR: classf can only be an array when using --pbb_thresholds."
            print "ERROR: Can use more than one 'classf' only when using --pbb_thresholds."
            exit()
        y = tbdata.field (classf[0])
    return y
