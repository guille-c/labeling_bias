import numpy as np
#import pylab as pl
import pyfits as pf
import argparse
import os

from time import time
from sklearn.neighbors import KernelDensity
from sklearn.grid_search import GridSearchCV
from sklearn.preprocessing import StandardScaler

parser = argparse.ArgumentParser(description='Create simulations following fits file parameters distribution.')

parser.add_argument ("in_file", help = "input fits file.")
parser.add_argument ("out_dir", help = "output fits file.")
parser.add_argument ("parameters", nargs = "+", help = "columns with parameters to follow.")
parser.add_argument ("--N_objs_in", type = int, default = -1,
                     help = "number of objects used for the KDE.")
parser.add_argument ("--N_objs_out", type = int, default = -1,
                     help = "number of objects of simulation.")

args = parser.parse_args()

if not os.path.exists(args.out_dir):
    os.makedirs(args.out_dir)
out_file = args.out_dir + "Simulations_" + str(args.N_objs_in) + "_" + str(args.N_objs_out) + ".fits"

if out_file == args.in_file:
    print "ERROR: in and out files are the same!"
    exit()
hdu_in = pf.open(args.in_file)[1]
tbdata = hdu_in.data # Open fits table file.
# Remove NANs
print tbdata.shape
tbdata = tbdata[~np.isnan(np.asarray(tbdata.tolist())).any(axis=1)]
# Remove INFs
tbdata = tbdata[~np.isinf(np.asarray(tbdata.tolist())).any(axis=1)]
print tbdata.shape

if args.N_objs_in > 0:
    i_s = np.arange (len(tbdata))
    np.random.shuffle(i_s)
    tbdata = tbdata[i_s[:args.N_objs_in]]
    hdu_in.data = tbdata
    hdu_in.writeto(args.out_dir + "input_" + str(args.N_objs_in) + "_" + str(args.N_objs_out) + ".fits")

if args.N_objs_out < 1:
    N_objs_out = tbdata.data.shape[0]
else:
    N_objs_out = args.N_objs_out

data = []
for par in args.parameters:
    data.append (tbdata.field (par))
data = np.array(data).transpose()
print "data shape = ", data.shape

#Mean removal and variance scaling
scaler = StandardScaler().fit(data)

# use grid search cross-validation to optimize the bandwidth
params = {'bandwidth': np.logspace(-2, 0, 10)}
grid = GridSearchCV(KernelDensity(), params)
print "Fitting KDE"
t1 = time()
grid.fit(scaler.transform(data))
t2 = time() - t1
print "fitting time = ", t2
print("best bandwidth: {0}".format(grid.best_estimator_.bandwidth))

# use the best estimator to compute the kernel density estimate
kde = grid.best_estimator_

# sample N_objs new points from the data
new_data = kde.sample(N_objs_out, random_state=0)
print "new data = ", new_data 
new_data = new_data * scaler.std_ + scaler.mean_
print scaler.std_, scaler.mean_
print "new data scaled = ", new_data

cols = []
for i in range (data.shape[1]):
    cols.append (pf.Column (name = args.parameters[i], array = new_data[:, i],
                            format = "D"))
tbhdu = pf.new_table (pf.ColDefs(cols))
tbhdu.writeto (out_file, clobber = True)
