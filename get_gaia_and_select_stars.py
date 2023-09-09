import os
import sys
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
import math
import matplotlib.cm as cm
import matplotlib.patches as patches
import mw_plot
from mw_plot import center_radec, anti_center_radec, northpole_radec, southpole_radec  # constants
from mw_plot import mw_radec # milkyway plane in RA/DEC
import astropy.units as u
import matplotlib.colors as clrs
from scipy.stats import gaussian_kde
import corner
import csv

#pylint: disable=wrong-import-position
#pylint: disable=import-error, C0103, R0914, W0311, C0114, C0301, R0903, C0304, C0200, R1705, R0911, R0912, R0915, R1702, R1710, W1401, C0209
import numpy as np
import pyvo as vo
from astroquery.ukidss import Ukidss
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS
import astropy.units as u


class GetData():

    def __init__(self, ra, dec):
        self.ra, self.dec = ra, dec


    def get_gaia_data(self):
        """
        `irgsctool.GetData.get_gaia_data()`

        <justify> This function sends a query to obtain GAIA DR3 data
        using the astroquery module.
        The ROW_LIMIT is set to -1 which implies that the query retriees all
        the rows in the given field.</justify>
        Raises:
                ValueError if there is no observed GAIA DR3 data for the
                            given set of input coordinates.
        Returns:   VOTable containing Gaia data.
        """
        ra_name = str(self.ra).replace('.','_')
        dec_name = str(self.dec).replace('.', '_')
        file_name = 'GAIA' + '_' + 'RA'+str(ra_name)\
                    + 'DEC' + str(dec_name) + '.csv'
        if file_name is True:
            table = np.genfromtxt(str(file_name),\
                        delimiter=',', skip_header=1)
        #tables = Gaia.load_tables(only_names=True)
        #for table in (tables):
            #print (table.get_qualified_name())
        from astroquery.gaia import Gaia
        job = Gaia.launch_job("select "
                      "phot_g_mean_mag,teff_gspphot,ra,dec,ag_gspphot,"
                      "parallax,mh_gspphot,logg_gspphot,"
                      "phot_bp_mean_mag, phot_rp_mean_mag"
                      "from gaiadr3.gaia_source order by source_id")
        r = job.get_results()
        r.save('data')
        return r

hdulist = fits.open('data.fits')
p = hdulist[1].data
gmag = p['phot_g_mean_mag']
teff = p['teff_gspphot']
parag = p['parallax']
parag = parag*1e-3
#eparag = p['parallax_error']
#eparag = eparag*1e-3
ag = p['ag_gspphot']
ra = p['ra']
dec = p['dec']
feh = p['mh_gspphot']
logg = p['logg_gspphot']
bp = p['phot_bp_mean_mag']
rp = p['phot_rp_mean_mag']
ind = np.where((parag>0) & (ag > 0))[0]
print('ra,dec=', ra, dec)
gmag = np.asarray(gmag[ind],dtype = np.float32)
teff = np.asarray(teff[ind],dtype = np.float32)
feh = np.asarray(feh[ind],dtype = np.float32)
logg = np.asarray(logg[ind],dtype = np.float32)
parag = np.asarray(parag[ind],dtype = np.float32)
#eparag = np.asarray(eparag[ind],dtype = np.float32)
ag = np.asarray(ag[ind],dtype = np.float32)
ra = np.asarray(ra[ind],dtype = np.float32)
dec = np.asarray(dec[ind],dtype = np.float32)
bp = np.asarray(bp[ind],dtype = np.float32)
rp = np.asarray(rp[ind],dtype = np.float32)


plt.clf()
cr = np.arange(len(teff_rgb))
xy = np.vstack([teff_rgb, (lumg_rgb)])
print('xy=', xy)
z = gaussian_kde(xy)(xy)
print('z=', z)
bins2 = np.arange(lumg_rgb.min(),(lumg_rgb).max()+.1, 0.1)
plt.clf()
fig = plt.figure(figsize=(8,8))
from matplotlib.gridspec import GridSpec
gs = GridSpec(4, 4)
ax_joint = fig.add_subplot(gs[1:4,0:3])
ax_marg_x = fig.add_subplot(gs[0,0:3])
ax_marg_y = fig.add_subplot(gs[1:4,3])
ax_joint.scatter(teff_rgb, (lumg_rgb), s=1, c = z,cmap = 'YlGn', norm= clrs.Normalize(), alpha = 0.5)
ax_joint.invert_xaxis()
#ax_joint.scatter(teff_rgb1, lumg_rgb1, s=1, c = z,cmap = 'red', norm= clrs.Normalize(), alpha = 0.5)
x1 = [4800, 4600]
y1 = [1.625, 1.775]
x2 = [4600, 4600]
y2 = [1.625, 1.775]
x3 = [4600, 4800]
y3 = [1.775, 1.775]
x4 = [4800, 4800]
y4 = [1.775, 1.625]
plt.grid()
#ax_joint.plot(x1, y1, linewidth=2, color='r', linestyle = 'dashed')
#ax_joint.plot(x2, y2, linewidth=2, color='r', linestyle = 'dashed')
#ax_joint.plot(x3, y3, linewidth=2, color='r', linestyle = 'dashed')
#ax_joint.plot(x4, y4, linewidth=2, color='r', linestyle = 'dashed')
ax_joint.grid()
#ax_joint.set_ylim(-2.5,5.0)
#ax_joint.set_xlim(6000,3500)
#ax.joint.set_yscale("log")
#ax_joint.legend(fontsize=8, loc = 4)
nx, bx, px = ax_marg_x.hist(teff_rgb, color = 'm', edgecolor = 'g', alpha = 0.5, label = '$T_{eff}$')
ny, by, px = ax_marg_y.hist((lumg_rgb), bins = bins2, orientation="horizontal", edgecolor = 'g', alpha = 0.5, facecolor = 'orange', label = 'Luminosity')
biny_max = np.where(ny == ny.max())
print('maxbin', "{:.2f}".format(by[biny_max][0]))
#plt.text(500, -1.5, str('mode at '"{:.2f}".format(by[biny_max][0])), rotation = 270)
#ax_marg_y.set_ylim(1.3,3.0)
ax_marg_x.grid()
#ax_marg_x.legend()
ax_marg_y.grid()
#ax_marg_y.legend(fontsize=8)
# Turn off tick labels on marginals
plt.setp(ax_marg_x.get_xticklabels(), visible=False)
plt.setp(ax_marg_y.get_yticklabels(), visible=False)
# Set labels on joint
ax_joint.set_xlabel('$T_{eff}$')
ax_joint.set_ylabel('log(Luminosity)')
# Set labels on marginals
ax_marg_y.set_xlabel('Counts')
ax_marg_x.set_ylabel('Counts')
plt.savefig('all_rc_stars_gaia_colour_hrd.png')
plt.clf()

teff_sun = np.asarray([5772], dtype = np.float32)
bolo_sun = np.asarray([4.74], dtype = np.float32)
logg_sun = np.asarray([4.44], dtype = np.float32)

a0 = np.asarray([0.06], dtype = np.float32)
a1 = np.asarray([6.731*1e-5], dtype = np.float32)
a2 = np.asarray([-6.647*1e-8], dtype = np.float32)
a3 = np.asarray([2.859*1e-11], dtype = np.float32)
a4 = np.asarray([-7.197*1e-15], dtype = np.float32)

with open('low_mass_red_giants_list.txt', 'w') as file0:
    for i in range(len(teff)):
        M = gmag[i] - 5*np.log10((1/parag[i])) + 5 - ag[i]
        bolo_g = np.sum([a1*(teff[i] - teff_sun)] + [a2*(teff[i] - teff_sun)**2] + [a3*(teff[i] - teff_sun)**3] + [a4*(teff[i] - teff_sun)**4] + a0)
        lum_g = 10**((-2.5**-1)*(M + bolo_g - bolo_sun))
        mass_g = 10**(np.log10(lum_g) + logg[i] - logg_sun + (4*np.log10(teff_sun/teff[i])))

        if mass_g <= 2.0:
                file0.write('%0.6f %0.6f %0.2f %0.4f %0.4f %0.1f %0.1f %0.2f %0.2f %0.2f %0.2f\n' %(ra[i], dec[i], teff[i], lum_g, np.log10(mass_g), feh[i], logg[i], bp[i], rp[i], gmag[i], M))

p0 = np.genfromtxt('low_mass_red_giants_list.txt')
teff_rgb = p0[:,2]
lumg_rgb = p0[:,3]
lumg_rgb = np.log10(lumg_rgb)
ra = p0[:,0]
dec = p0[:,1]
mass_rgb = p0[:,4]
feh_rgb = p0[:,5]
logg_rgb = p0[:,6]
bp_rgb = p0[:,7]
rp_rgb = p0[:,8]
gmag_rgb = p0[:,9]
Mg_rgb = p0[:,10]

pts0 = np.where((teff_rgb<5300) & (teff_rgb>4300) & ((lumg_rgb)>1.5) &(lumg_rgb<2.0) )[0]
#pts0 = np.where((teff_rgb>4400) & ((lumg_rgb)>2.0))[0]
print('Number of interesting stars=', len(pts0))

ra1 = ra[pts0]
dec1 = dec[pts0]
teff_rgb1 = teff_rgb[pts0]
lumg_rgb1 = lumg_rgb[pts0]
mass_rgb1 = mass_rgb[pts0]
feh_rgb1 = feh_rgb[pts0]
logg_rgb1 = logg_rgb[pts0]
bp_rgb1 = bp_rgb[pts0]
rp_rgb1 = rp_rgb[pts0]
gmag_rgb1 = gmag_rgb[pts0]
Mg_rgb1 = Mg_rgb[pts0]


ra1 = ra1*u.degree
dec1 = dec1*u.degree

ra1 = ra1/u.degree
dec1 = dec1/u.degree

header = ['ra', 'dec', '$T_{eff}$', 'Luminosity', 'Mass', '[Fe/H]', 'log(g)', '$Gaia_{bp}$', '$Gaia_{rp}$', '$M_{G}$']

with open('rc.csv', 'w', encoding='UTF8') as file1:
    writer=csv.writer(file1)
    writer.writerow(header)
    for i in range(len(ra1)):
        #file1.write('%0.6f %0.6f %0.2f %0.4f %0.4f %0.1f %0.2f %0.2f %0.2f %0.2f\n' %(ra1[i], dec1[i], teff_rgb1[i], lumg_rgb1[i], mass_rgb1[i], feh_rgb1[i], logg_rgb1[i], bp_rgb1[i], rp_rgb1[i], Mg_rgb1[i]))
        data = ra1[i], dec1[i], teff_rgb1[i], lumg_rgb1[i], mass_rgb1[i], feh_rgb1[i], logg_rgb1[i], bp_rgb1[i], rp_rgb1[i], Mg_rgb1[i]
        writer.writerow(data)

rgbf = np.genfromtxt('rc.csv', delimiter=',', skip_header=1)
bpf = rgbf[:,7]
rpf = rgbf[:,8]
Mgf = rgbf[:,9]
teff_rgbf = rgbf[:,2]
lumg_rgbf = rgbf[:,3]


plt.clf()
plt.scatter((bp_rgb - rp_rgb), Mg_rgb, s=1, color='green', alpha=0.1)
plt.scatter((bpf - rpf), Mgf, s=1, color='red', alpha=0.1, label='RC stars')
plt.xlabel('Bp-Rp colour', fontsize=18)
plt.ylabel('Absolute G magnitude', fontsize=18)
plt.grid()
plt.savefig('gaia_colour_hrd.png')
plt.clf()

x = teff_rgb
y = lumg_rgb

cr = np.arange(len(teff_rgbf))
xy = np.vstack([teff_rgbf,(lumg_rgbf)])
print('xy=', xy)
z = gaussian_kde(xy)(xy)
print('z=', z)
bins2 = np.arange((lumg_rgbf).min(), (lumg_rgbf).max()+.1, 0.1)
plt.clf()
fig = plt.figure()
from matplotlib.gridspec import GridSpec
gs = GridSpec(4, 4)
ax_joint = fig.add_subplot(gs[1:4,0:3])
ax_marg_x = fig.add_subplot(gs[0,0:3])
ax_marg_y = fig.add_subplot(gs[1:4,3])
ax_joint.scatter(teff_rgbf, (lumg_rgbf), s=1, c = z,cmap = 'YlGn', norm= clrs.Normalize(), alpha = 0.5)
#ax_joint.scatter(teff_rgb1, lumg_rgb1, s=1, c = z,cmap = 'red', norm= clrs.Normalize(), alpha = 0.5)
x1 = [4800, 4600]
y1 = [1.625, 1.775]
x2 = [4600, 4600]
y2 = [1.625, 1.775]
x3 = [4600, 4800]
y3 = [1.775, 1.775]
x4 = [4800, 4800]
y4 = [1.775, 1.625]
plt.grid()
ax_joint.invert_xaxis()
#ax_joint.plot(x1, y1, linewidth=2, color='r', linestyle = 'dashed')
#ax_joint.plot(x2, y2, linewidth=2, color='r', linestyle = 'dashed')
#ax_joint.plot(x3, y3, linewidth=2, color='r', linestyle = 'dashed')
#ax_joint.plot(x4, y4, linewidth=2, color='r', linestyle = 'dashed')
#ax_joint.grid()
#ax_joint.set_ylim(1.4,2.0)
#ax_joint.set_xlim(5150,4500)
#ax_joint.legend(fontsize=8, loc = 4)
nx, bx, px = ax_marg_x.hist(teff_rgbf, color = 'm', edgecolor = 'g', alpha = 0.5, label = '$T_{eff}$')
ny, by, px = ax_marg_y.hist((lumg_rgbf), bins = bins2, orientation="horizontal", edgecolor = 'g', alpha = 0.5, facecolor = 'orange', label = 'Luminosity')
biny_max = np.where(ny == ny.max())
print('maxbin', "{:.2f}".format(by[biny_max][0]))
#plt.text(500, -1.5, str('mode at '"{:.2f}".format(by[biny_max][0])), rotation = 270)
#ax_marg_y.set_ylim(1.3,3.0)
ax_marg_x.grid()
#ax_marg_x.legend()
ax_marg_y.grid()
#ax_marg_y.legend(fontsize=8)
# Turn off tick labels on marginals
plt.setp(ax_marg_x.get_xticklabels(), visible=False)
plt.setp(ax_marg_y.get_yticklabels(), visible=False)
# Set labels on joint
ax_joint.set_xlabel('$T_{eff}$(K)')
ax_joint.set_ylabel('Log(L/$L_{\odot}$)')
# Set labels on marginals
ax_marg_y.set_xlabel('Counts')
ax_marg_x.set_ylabel('Counts')
plt.savefig('scatter_with_hist.png')
plt.clf()

sys.exit(0)
with open('lum1.txt', 'w') as file0:
    for i in range(len(teff)):
        M = gmag[i] - 5*np.log10((1/parag[i])) + 5 - ag[i]
        bolo_g = np.sum([a1*(teff[i] - teff_sun)] + [a2*(teff[i] - teff_sun)**2] + [a3*(teff[i] - teff_sun)**3] + [a4*(teff[i] - teff_sun)**4] + a0)
        lum_g = 10**((-2.5**-1)*(M + bolo_g - bolo_sun))
        mass_g = 10**(np.log10(lum_g) + logg[i] - logg_sun + (4*np.log10(teff_sun/teff[i])))

        if mass_g < 2.0:
                file0.write('%0.6f %0.6f %0.2f %0.4f %0.4f %0.1f %0.1f %0.2f %0.2f %0.2f %0.2f\n' %(ra[i], dec[i], teff[i], lum_g, np.log10(mass_g), feh[i], logg[i], bp[i], rp[i], gmag[i], M))

print(len(teff), len(gmag))
