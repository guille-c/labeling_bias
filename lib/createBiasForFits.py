import numpy as np
import pylab as pl
import pyfits as pf
import astropy as ap
import astropy.cosmology as cosmo
import argparse

def processSelection (hdu, i_ini, i_end, field_r, field_m, field_z, field_psf):
    tbdata = hdu.data[i_ini:i_end]
    data_aux = np.asarray(hdu_in.data.tolist())
    rm_crit = np.isnan (data_aux).any(axis=1) | np.isinf (data_aux).any(axis=1) | (hdu_in.data.field(field_z) <= 0) | (hdu_in.data.field(field_r) <= 0)
    del data_aux

cosmo.core.set_current(cosmo.WMAP9)

parser = argparse.ArgumentParser(description='Simulate morphology bias for fits table.')
parser.add_argument ("in_file", help = "input fits file.")
parser.add_argument ("out_dir", help = "output fits file.")

args = parser.parse_args()

field_r = "petroRad_r"
field_m = "corrMag_r"
field_z = "z"
field_psf = "psfWidth_r"

hdu_in = pf.open(args.in_file)[1] # Open fits table file.
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

