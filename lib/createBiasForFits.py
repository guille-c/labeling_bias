import numpy as np
#import pylab as pl
import pyfits as pf
import astropy as ap
import astropy.cosmology as cosmo
import argparse
import os

def processSelection (hdu, i_ini, i_end, field_r, field_m, field_z, field_psf, 
                      p_class1, biases, out_dir):
    out_fn = out_dir + "sim_bias_" + str(i_ini) + "_" + str(i_end) + ".fits"
    if os.path.exists (out_fn):
        return
    print "selecting data"
    tbdata = hdu.data[i_ini:i_end]
    print "selected"
    data_aux = tbdata.tolist()
    print "tbdata.tolist() done"
    data_aux = np.asarray(tbdata.tolist())
    print "creating rm_crit"
    rm_crit = np.isnan (data_aux).any(axis=1) | np.isinf (data_aux).any(axis=1) | (tbdata.field(field_z) <= 0) | (tbdata.field(field_r) <= 0)
    print "is_nan = ", np.isnan (data_aux).any(axis=1).sum()
    print "is_inf = ", np.isinf (data_aux).any(axis=1).sum()
    print field_z, " <= 0 = ", (tbdata.field(field_z) <= 0).sum()
    print field_r, " <= 0 = ", (tbdata.field(field_r) <= 0).sum(), tbdata.field(field_r).max()
    print "rm_crit.sum = ", rm_crit.sum()
    del data_aux
    print tbdata.shape
    tbdata = tbdata[~rm_crit]
    print tbdata.shape

    r = tbdata.field(field_r)
    m = tbdata.field(field_m)


    # Uniform distributions to test
    # r = np.random.random (tbdata.shape[0])
    # m = np.random.random (tbdata.shape[0])*30

    z = tbdata.field(field_z)
    psf = tbdata.field(field_psf)
    y = (np.random.random(tbdata.shape[0]) < p_class1)

    d = cosmo.angular_diameter_distance(z).to(ap.units.kpc) # distance in kpc
    scale = d*ap.units.arcsec.to(ap.units.radian) # scale [kpc/arcsec]
    
    R = r * scale
    M = m - 5*np.log10(4.28E8*z)
    r_psf = 2.*np.sqrt(2*np.log(2)) * r/psf
    
    cols = [pf.Column (name = field_r, array = r, format = "D"), 
            pf.Column (name = field_m, array = m, format = "D"), 
            pf.Column (name = field_z, array = z, format = "D"), 
            pf.Column (name = field_psf, array = psf, format = "D"), 
            pf.Column (name = "absPetroMag_r", array = M, format = "D"), 
            pf.Column (name = "petroRad_r_psf", array = r_psf, format = "D"),
            pf.Column (name = "petroRad_r_kpc", array = R, format = "D"),
            pf.Column (name = "class", array = y, format = "B")]

    # Create biased labels
    r_m = np.median(r)
    m_m = np.median(m)-17.8
    crit1 = (y == 1)
    r1 = r[crit1]
    m1 = m[crit1]-17.8
    factor = np.exp(-(r1*r1/(2*r_m*r_m) + m1*m1/(2*m_m*m_m)))
    factor [m1 >= 0] = np.exp(-(r1[m1 >= 0]**2/(2*r_m*r_m)))
    uniform = np.random.random(crit1.sum())
    for bias in biases:
        print "bias = ", bias
        p_b = factor**(1./bias**2)
        y_b = y.copy()
        y_b[crit1] = y_b[crit1] * (uniform >= p_b)
        cols = cols + [pf.Column (name = "class_bias" + str(bias), 
                                  array = y_b, format = "B")]
    hdu = pf.BinTableHDU.from_columns(cols)
    hdu.writeto (out_fn, clobber = True)

def concatenate (iters, fields_names, out_dir):
    conc_data = []
    for i in range (len(iters) - 1):
        print "concatenate i = ", i
        i_ini, i_end = iters[i], iters[i + 1]
        fn = out_dir + "sim_bias_" + str(i_ini) + "_" + str(i_end) + ".fits"
        hdu = pf.open (fn)
        print hdu
        #if 'conc_data' in locals():
        conc_data = conc_data + hdu[1].data.tolist()
        cols = []
    conc_data = np.array(conc_data).transpose()
    print conc_data.shape
    for i in range(len(hdu[1].columns)):
        print i
        cols.append(pf.Column (name = hdu[1].columns.names[i], 
                               array = conc_data[i], format = "D"))
        #print new_hdu.data.field(i)[nrows:].shape, hdu[1].data.field(i).shape
        #new_hdu.data.field(i)[nrows:] = hdu[1].data.field(i)[:]
    new_hdu = pf.new_table(cols)
    new_hdu.writeto (out_dir + "all_sims_biased_" + str(conc_data.shape[1]) + ".fits", clobber = True)

cosmo.core.set_current(cosmo.WMAP9)

parser = argparse.ArgumentParser(description='Simulate morphology bias for fits table.')
parser.add_argument ("in_file", help = "input fits file.")
parser.add_argument ("out_dir", help = "output directory for fits files.")
parser.add_argument ("--N_iter_lines", type = int, default = 10000, 
                     help = "Aprox number of lines per sub-iteration.")
parser.add_argument ("--p1", type = float, default = 0.5, 
                     help = "Probability of an object of being class 1.")
parser.add_argument ("--biases", default = [0.1], type = float, nargs = "+",
                     help = "Probability of an object of being class 1.")
parser.add_argument ("--N_max", type = int,
                     help = "Maximum number of objects (rows) to process.")

args = parser.parse_args()

biases = np.array(args.biases, dtype = float)

field_r = "petroRad_r"
field_m = "corrMag_r"
field_z = "z"
field_psf = "psfWidth_r"

hdu_in = pf.open(args.in_file)[1] # Open fits table file.
#print hdu_in.data.field(field_r).max()
#exit()
if args.N_max:
    N_tot = args.N_max
else:
    N_tot = hdu_in.data.shape[0]
N_iters = int (np.round(1. * N_tot / args.N_iter_lines)) + 1
iters = np.array (np.round(np.linspace (0, N_tot, N_iters)), dtype = int)
for i in range (len(iters) - 1):
    print "PROCESSING ", iters[i], iters[i + 1]
    processSelection (hdu_in, iters[i], iters[i + 1], 
                      field_r, field_m, field_z, field_psf, 
                      args.p1, biases, args.out_dir)
    #(hdu, i_ini, i_end, field_r, field_m, field_z, field_psf, 
    # p_class1, biases, out_dir)
fields_names = [field_r, field_m, field_z, field_psf]
fields_names = fields_names + ["absPetroMag_r", "petroRad_r_psf", "petroRad_r_kpc"]

print "CONCATENATING"
concatenate (iters, fields_names, args.out_dir)
exit()

# Remove NANs, INFs, z <= 0 and r <= 0
print "Original data shape: ", hdu_in.data.shape
data_aux = np.asarray(hdu_in.data.tolist())
rm_crit = np.isnan (data_aux).any(axis=1) | np.isinf (data_aux).any(axis=1) | (hdu_in.data.field(field_z) <= 0) | (hdu_in.data.field(field_r) <= 0)
del data_aux
hdu_in.data = hdu_in.data[~rm_crit]
#tbdata = tbdata[~np.isnan(np.asarray(tbdata.tolist())).any(axis=1)]
tbdata = hdu_in.data 
print "Removed NANs and INFs: ", tbdata.shape

r = tbdata.field(field_r)
m = tbdata.field(field_m)
z = tbdata.field(field_z)
psf = tbdata.field(field_psf)

d = cosmo.angular_diameter_distance(z).to(ap.units.kpc) # distance in kpc
scale = d*ap.units.arcsec.to(ap.units.radian) # scale [kpc/arcsec]

R = r * scale
M = m - 5*np.log10(4.28E8*z)
r_psf = 2.*np.sqrt(2*np.log(2)) * r/psf

print hdu_in.data.field(field_z).shape
cols = [pf.Column (name = field_r, array = r, format = "D"), 
        pf.Column (name = field_m, array = m, format = "D"), 
        pf.Column (name = field_z, array = z, format = "D"), 
        pf.Column (name = field_psf, array = psf, format = "D"), 
        pf.Column (name = "absPetroMag_r", array = M, format = "D"), 
        pf.Column (name = "petroRad_r_psf", array = r_psf, format = "D"),
        pf.Column (name = "petroRad_r_kpc", array = R, format = "D")]
#t = hdu_in.columns + pf.ColDefs(cols)
hdu = pf.new_table(cols)
hdu.writeto (args.out_file, clobber = True)

