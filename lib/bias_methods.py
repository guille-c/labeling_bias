import numpy as np
import os
from glob import glob

def saveLs (dataAn, intrinsic, observables, y, labels, bins,
            N_objs, N_calc, increasing_bias = True, out_path = "./", txt = ""):
    print bins, N_objs
    Ls = np.zeros((len(N_objs), bins.shape[0]))
    Ls_std = np.zeros((len(N_objs), bins.shape[0]))
    for i in np.arange (len(bins)):
        for j in np.arange (len(N_objs)):
            print N_objs[j]            
            filenameLN = (out_path + "L_N" + txt + str(bins[i]) + "_" + 
                          str(N_objs[j]) + "_" + str(N_calc))
            print "filenameLN = ", filenameLN
            max_calc = 0
            if (os.path.exists(filenameLN + ".dat") and 
                os.path.exists(filenameLN + ".npy")):
                # read from file
                file = open (filenameLN + ".dat", "r")
                file.readline()
                L1, L1_s = file.readline().split()[:2]
                print bins[i], N_objs[j], L1
                Ls[j, i], Ls_std[j, i] = float(L1), float(L1_s)
            else:
                # find largest file already created
                prefix = out_path + "L_N" + txt + "[[]" + str(bins[i]).split("]")[0][1:] + \
                    "[]]_" + str(N_objs[j]) + "_"
                files = glob(prefix + "*.npy")
                print files
                if len(files) > 0:
                    for f in files:
                        n_cand = int(f.split("_")[-1][:-4])
                        if n_cand > max_calc:
                            max_calc = n_cand
                    print "max_calc = ", max_calc
                    filenameMaxCalc = (out_path + "L_N" + txt + str(bins[i]) + "_" + 
                                       str(N_objs[j]) + "_" + str(max_calc)) + ".npy"
                    calcLs = np.load(filenameMaxCalc)
                    print calcLs

                print "Calculating L for ", N_objs[j], bins[i]
                L1, N1 = dataAn.getRandomL (intrinsic, observables, y, labels, 
                                            increasing_bias, N_calc - max_calc, 
                                            bins[i][0], bins[i][1],
                                            minElementsBin = N_objs[j], 
                                            N_objs = N_objs[j]*2**bins[i][0]*bins[i][1],
                                            bootstrap = True, 
                                            kd_tree = "highest_fraction_difference")
                if max_calc > 0:
                    L1 = np.concatenate((calcLs[0], L1))
                    N1 = np.concatenate((calcLs[1], N1))
                print filenameLN, L1, len(L1)
                file = open (filenameLN + ".dat", "w")
                file.write ("L L_std N\n" + str(L1.mean()) + " " + str(L1.std()) + 
                            " " + str(N1.mean()) + "\n")
                file.close ()
                np.save(filenameLN, np.array([L1, N1]))
                Ls[j, i], Ls_std[j, i] = L1.mean(), L1.std()
    return Ls, Ls_std

def readLs (bins, N_objs, N_calc, out_path = "./", txt = ""):
    print bins, N_objs
    Ls = np.zeros(bins.shape[0])
    Ls_std = np.zeros(bins.shape[0])
    for i in np.arange (len(bins)):
        filenameLN = (out_path + "L_N" + txt + str(bins[i]) + "_" + 
                      str(N_objs) + "_" + str(N_calc) + ".dat")
        if os.path.exists(filenameLN):
            # read from file
            file = open (filenameLN, "r")
            file.readline()
            L1, L1_s = file.readline().split()[:2]
            print bins[i], N_objs, L1
            Ls[i], Ls_std[i] = float(L1), float(L1_s)
        else:
            print "ERROR bias_methods.readLs: file ", filenameLN, "does not exist."
            exit()
    return Ls, Ls_std

def readAllLsNpy (bins, N_objs, N_calc, out_path, txt = ""):
    Ls = []
    for i in np.arange (len(bins)):
        filenameLN = (out_path + "L_N" + txt + str(bins[i]) + "_" + 
                       str(N_objs) + "_" + str(N_calc) + ".npy")
        if os.path.exists(filenameLN):
            # read from file
            Ls1, N1 = np.load (filenameLN)
            Ls.append(Ls1)
        else:
            print "ERROR bias_methods.readLsDifferences: file ", filenameLN1, filenameLN2, "does not exist."
            Ls.append (np.array([0]))
    return np.array(Ls)

def readLsDifferences (bins, N_objs, N_calc, out_path1, out_path2, txt = ""):
    print bins, N_objs
    Ls_diff = np.zeros(bins.shape[0])
    Ls_diff_std = np.zeros(bins.shape[0])
    for i in np.arange (len(bins)):
        filenameLN1 = (out_path1 + "L_N" + txt + str(bins[i]) + "_" + 
                       str(N_objs) + "_" + str(N_calc) + ".npy")
        filenameLN2 = (out_path2 + "L_N" + txt + str(bins[i]) + "_" + 
                       str(N_objs) + "_" + str(N_calc) + ".npy")
        if os.path.exists(filenameLN1) and os.path.exists(filenameLN2):
            # read from file
            Ls1, N1 = np.load (filenameLN1)
            Ls2, N2 = np.load (filenameLN2)
            Ls_diffs = []
            for i1 in range (len(Ls1)):
                for i2 in range (len(Ls2)):
                    Ls_diffs.append([Ls1[i1] - Ls2[i2]])
            Ls_diff[i], Ls_diff_std[i] = np.mean(Ls_diffs), np.std(Ls_diffs)
            #Ls_diff[i], Ls_diff_std[i] = (Ls1 - Ls2).mean(), (Ls1 - Ls2).std()
        else:
            print "ERROR bias_methods.readLsDifferences: file ", filenameLN1, filenameLN2, "does not exist."
            Ls_diff[i] = 0
            Ls_diff_std[i] = 1e-6
            #exit()
    return Ls_diff, Ls_diff_std

# Get labels / probabilities column.
def createLabels (tbdata, classf, pbb_thresholds):
    N_tot = len(tbdata)
    if not np.isscalar(pbb_thresholds):
        if not (pbb_thresholds is None) and len (classf) != len (pbb_thresholds):
            print "ERROR: the number of pbb fields does not match the number of thresholds."
            exit()
        y = np.zeros(N_tot)
        if len(classf) == 1:
            y[tbdata.field (classf[0]) >= pbb_thresholds[0]] = 1
            y[1 - tbdata.field (classf[0]) > pbb_thresholds[0]] = 2
        else:
            for i in range(len(classf)):
                print classf[i], classf, tbdata.field (classf[i]).shape
                crit = (tbdata.field (classf[i]) >= pbb_thresholds[i])
                print i, crit.sum()
                y[crit] = i + 1
    else:
        if len(classf) != 1:
            #print "ERROR: classf can only be an array when using --pbb_thresholds."
            print "ERROR: Can use more than one 'classf' only when using --pbb_thresholds."
            exit()
        y = tbdata.field (classf[0])
    return y
