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
import matplotlib
import warnings
warnings.filterwarnings("ignore")

rgb_list = np.genfromtxt('lum2.txt')
ra = rgb_list[:,0]
#ra = np.arange(2609,ra[-1],1)

for i0 in range(len(ra)):
    print('i0=', i0)
    name1 = 'g'+str(i0)
    name2 = 'r'+str(i0)
    name3 = 'i'+str(i0)
    try:
        p0 = np.genfromtxt(name1+'.csv', delimiter=',', skip_header=1)
        p1 = np.genfromtxt(name2+'.csv', delimiter=',', skip_header=1)
        p3 = np.genfromtxt(name3+'.csv', delimiter=',', skip_header=1)
    except FileNotFoundError as not_found:
        print('not_found_file_name=', not_found.filename)
    #if (os.path.isfile(p0)) is not True or (os.path.isfile(p1)) is not True or  (os.path.isfile(p3)) is not True:
        continue
    if len(p0) ==0.0 and len(p1) ==0.0 and len(p3)==0.0:
        continue

    elif len(p0) == 0.0 and len(p1) !=0.0 and len(p3) == 0.0:
        if len(p1)>24:
            objid_r = p1[:,0];hjd_r = p1[:,2];mag_r = p1[:,4];emag_r = p1[:,5];hjd_r = hjd_r - 2450000
            ptsr = np.where((mag_r/emag_r)>5.0)[0]
            objid_r = objid_r[ptsr]; emag_r = emag_r[ptsr]; mag_r = mag_r[ptsr]; hjd_r = hjd_r[ptsr]
            plt.clf()
            #plt.errorbar(hjd_g, mag_g, yerr=emag_g, fmt='o', markersize=1, color='purple', alpha=0.3, label='g band')
            plt.errorbar(hjd_r, mag_r, yerr=emag_r, fmt='o', markersize=1, color='red', alpha=0.3, label='r band')
            plt.xlabel('Time (HJD)')
            plt.ylabel('${Mag}$')
            ax = plt.gca()
            ax.invert_yaxis()
            plt.grid()
            plt.legend(loc='best')
            #plt.xlim(8866.6,8867)
            plt.savefig('/mnt/c/Users/sshah/Documents/eswar_anohita_project/9_12/mstar_gt_2msun/plots/lc_'+str(name2)+'.png')
            plt.clf()

    elif len(p0)!= 0.0 and len(p1) == 0.0 and len(p3) ==0.0:
        print('lenp0=', len(p0))
        try:
            objid_g = p0[:,0];hjd_g = p0[:,2];mag_g = p0[:,4];emag_g = p0[:,5];hjd_g = hjd_g - 2450000
            ptsg  = np.where((mag_g/emag_g)>5.0)[0]
            objid_g = objid_g[ptsg]; emag_g = emag_g[ptsg]; mag_g = mag_g[ptsg]; hjd_g = hjd_g[ptsg]
            plt.clf()
            plt.errorbar(hjd_g, mag_g, yerr=emag_g, fmt='o', markersize=1, color='green', alpha=0.3, label='g band')
            #plt.errorbar(hjd_r, mag_r, yerr=emag_r, fmt='o', markersize=1, color='yellow', alpha=0.3, label='r band')
            plt.xlabel('Time (HJD)')
            plt.ylabel('${Mag}$')
            ax = plt.gca()
            ax.invert_yaxis()
            plt.grid()
            plt.legend(loc='best')
            #plt.xlim(8866.6,8867)
            plt.savefig('/mnt/c/Users/sshah/Documents/eswar_anohita_project/9_12/mstar_gt_2msun/plots/lc_'+str(name1)+'.png')
            plt.clf()
        except IndexError:
            print('Index Error at', i0)
    elif len(p0)== 0.0 and len(p1) == 0.0 and len(p3) !=0.0 and len(p3) > 24:
        #print('Empty file r and g band', name2, name1)
        objid_i = p3[:,0];hjd_i = p3[:,2];mag_i = p3[:,4];emag_i = p3[:,5];hjd_i = hjd_i - 2450000
        ptsi  = np.where((mag_i/emag_i)>5.0)[0]
        objid_i = objid_i[ptsi]; emag_i = emag_i[ptsi]; mag_i = mag_i[ptsi]; hjd_i = hjd_i[ptsi]
        plt.clf()
        plt.errorbar(hjd_i, mag_i, yerr=emag_i, fmt='o', markersize=1, color='brown', alpha=0.3, label='i band')
        #plt.errorbar(hjd_r, mag_r, yerr=emag_r, fmt='o', markersize=1, color='yellow', alpha=0.3, label='r band')
        plt.xlabel('Time (HJD)')
        plt.ylabel('${Mag}$')
        ax = plt.gca()
        ax.invert_yaxis()
        plt.grid()
        plt.legend(loc='best')
        #plt.xlim(8866.6,8867)
        plt.savefig('/mnt/c/Users/sshah/Documents/eswar_anohita_project/9_12/mstar_gt_2msun/lc_'+str(name3)+'.png')
        plt.clf()

    elif len(p0)!= 0.0 and len(p1) == 0.0 and len(p3) !=0.0:
        #print('Empty file r band', name2)
        objid_g = p0[:,0];hjd_g = p0[:,2];mag_g = p0[:,4];emag_g = p0[:,5];hjd_g = hjd_g - 2450000
        objid_i = p3[:,0];hjd_i = p3[:,2];mag_i = p3[:,4];emag_i = p3[:,5];hjd_i= hjd_i - 2450000
        ptsi  = np.where((mag_i/emag_i)>5.0)[0]
        objid_i = objid_i[ptsi]; emag_i = emag_i[ptsi]; mag_i = mag_i[ptsi]; hjd_i = hjd_i[ptsi]
        ptsg  = np.where((mag_g/emag_g)>5.0)[0]
        objid_g = objid_g[ptsg]; emag_g = emag_g[ptsg]; mag_g = mag_g[ptsg]; hjd_g = hjd_g[ptsg]
        plt.clf()
        plt.errorbar(hjd_g, mag_g, yerr=emag_g, fmt='o', markersize=1, color='green', alpha=0.3, label='g band')
        plt.errorbar(hjd_i, mag_i, yerr=emag_i, fmt='o', markersize=1, color='brown', alpha=0.3, label='i band')
        plt.xlabel('Time (HJD)')
        plt.ylabel('${Mag}$')
        ax = plt.gca()
        ax.invert_yaxis()
        plt.grid()
        plt.legend(loc='best')
        #plt.xlim(8866.6,8867)
        plt.savefig('/mnt/c/Users/sshah/Documents/eswar_anohita_project/9_12/mstar_gt_2msun/lc_'+str(name1)+str(name3)+'.png')
        plt.clf()

    elif len(p0)!=0.0 and len(p1)!=0.0 and len(p3) == 0.0 and len(p0)>24 and len(p1)>24:
        #print('Empty i band=', name3)
        objid_g = p0[:,0];hjd_g = p0[:,2];mag_g = p0[:,4];emag_g = p0[:,5];hjd_g = hjd_g - 2450000
        objid_r = p1[:,0];hjd_r = p1[:,2];mag_r = p1[:,4];emag_r = p1[:,5];hjd_r = hjd_r - 2450000
        ptsg  = np.where((mag_g/emag_g)>5.0)[0]
        objid_g = objid_g[ptsg]; emag_g = emag_g[ptsg]; mag_g = mag_g[ptsg]; hjd_g = hjd_g[ptsg]
        ptsr = np.where((mag_r/emag_r)>5.0)[0]
        objid_r = objid_r[ptsr]; emag_r = emag_r[ptsr]; mag_r = mag_r[ptsr]; hjd_r = hjd_r[ptsr]
        plt.clf()
        plt.errorbar(hjd_g, mag_g, yerr=emag_g, fmt='o', markersize=1, color='green', alpha=0.3, label='g band')
        plt.errorbar(hjd_r, mag_r, yerr=emag_r, fmt='o', markersize=1, color='red', alpha=0.3, label='r band')
        plt.xlabel('Time (HJD)')
        plt.ylabel('$Magnitude$')
        ax = plt.gca()
        ax.invert_yaxis()
        plt.grid()
        plt.legend(loc='best')
        #plt.xlim(8866.6,8867)
        plt.savefig('/mnt/c/Users/sshah/Documents/eswar_anohita_project/9_12/mstar_gt_2msun/lc_'+str(name1)+'_'+str(name2)+'.png')
        plt.clf()

    elif len(p0)!=0.0 and len(p1)==0.0 and len(p3) != 0.0:
        #print('Empty r band=', name2)
        objid_g = p0[:,0];hjd_g = p0[:,2];mag_g = p0[:,4];emag_g = p0[:,5];hjd_g = hjd_g - 2450000
        objid_i = p3[:,0];hjd_i = p3[:,2];mag_i = p3[:,4];emag_i = p3[:,5];hjd_i = hjd_i - 2450000
        ptsg  = np.where((mag_g/emag_g)>5.0)[0]
        objid_g = objid_g[ptsg]; emag_g = emag_g[ptsg]; mag_g = mag_g[ptsg]; hjd_g = hjd_g[ptsg]
        objid_i = objid_i[ptsi]; emag_i = emag_i[ptsi]; mag_i = mag_i[ptsi]; hjd_i = hjd_i[ptsi]
        ptsg  = np.where((mag_g/emag_g)>5.0)[0]
        plt.clf()
        plt.errorbar(hjd_g, mag_g, yerr=emag_g, fmt='o', markersize=1, color='green', alpha=0.3, label='g band')
        plt.errorbar(hjd_i, mag_i, yerr=emag_i, fmt='o', markersize=1, color='brown', alpha=0.3, label='i band')
        plt.xlabel('Time (HJD)')
        plt.ylabel('$Magnitude$')
        ax = plt.gca()
        ax.invert_yaxis()
        plt.grid()
        plt.legend(loc='best')
        #plt.xlim(8866.6,8867)
        plt.savefig('/mnt/c/Users/sshah/Documents/eswar_anohita_project/9_12/mstar_gt_2msun/lc_'+str(name1)+'_'+str(name3)+'.png')
        plt.clf()

    elif len(p0)==0.0 and len(p1)!=0.0 and len(p3) != 0.0 and len(p1)>24 and len(p3) >24:
        objid_i = p3[:,0];hjd_i = p3[:,2];mag_i = p3[:,4];emag_i = p3[:,5];hjd_i = hjd_i - 2450000
        objid_r = p1[:,0];hjd_r = p1[:,2];mag_r = p1[:,4];emag_r = p1[:,5];hjd_r = hjd_r - 2450000
        ptsr = np.where((mag_r/emag_r)>5.0)[0]
        ptsi = np.where((mag_i/emag_i)>5.0)[0]
        objid_i = objid_i[ptsi]; emag_i = emag_i[ptsi]; mag_i = mag_i[ptsi]; hjd_i = hjd_i[ptsi]
        objid_r = objid_r[ptsr]; emag_r = emag_r[ptsr]; mag_r = mag_r[ptsr]; hjd_r = hjd_r[ptsr]
        plt.clf()
        plt.errorbar(hjd_i, mag_i, yerr=emag_i, fmt='o', markersize=1, color='brown', alpha=0.3, label='i band')
        plt.errorbar(hjd_r, mag_r, yerr=emag_r, fmt='o', markersize=1, color='red', alpha=0.3, label='r band')
        plt.xlabel('Time (HJD)')
        plt.ylabel('$Magnitude$')
        ax = plt.gca()
        ax.invert_yaxis()
        plt.grid()
        plt.legend(loc='best')
        #plt.xlim(8866.6,8867)
        plt.savefig('/mnt/c/Users/sshah/Documents/eswar_anohita_project/9_12/mstar_gt_2msun/lc_'+str(name3)+'_'+str(name2)+'.png')
        plt.clf()

    elif len(p0)!=0.0 and len(p1)!=0.0 and len(p3) != 0.0:
        if len(p0) > 24 and len(p1)>24 and len(p3)>24:
        #print('No files empty')
            objid_g = p0[:,0];hjd_g = p0[:,2];mag_g = p0[:,4];emag_g = p0[:,5];hjd_g = hjd_g - 2450000
            objid_r = p1[:,0];hjd_r = p1[:,2];mag_r = p1[:,4];emag_r = p1[:,5];hjd_r = hjd_r - 2450000
            objid_i = p3[:,0];hjd_i = p3[:,2];mag_i = p3[:,4];emag_i = p3[:,5];hjd_i = hjd_i - 2450000
            ptsi  = np.where((mag_i/emag_i)>5.0)[0]
            ptsg  = np.where((mag_g/emag_g)>5.0)[0]
            ptsr = np.where((mag_r/emag_r)>5.0)[0]
            objid_i = objid_i[ptsi]; emag_i = emag_i[ptsi]; mag_i = mag_i[ptsi]; hjd_i = hjd_i[ptsi]
            objid_r = objid_r[ptsr]; emag_r = emag_r[ptsr]; mag_r = mag_r[ptsr]; hjd_r = hjd_r[ptsr]
            objid_g = objid_g[ptsg]; emag_g = emag_g[ptsg]; mag_g = mag_g[ptsg]; hjd_g = hjd_g[ptsg]
            plt.clf()
            plt.errorbar(hjd_g, mag_g, yerr=emag_g, fmt='o', markersize=1, color='green', alpha=0.3, label='g band')
            plt.errorbar(hjd_r, mag_r, yerr=emag_r, fmt='o', markersize=1, color='red', alpha=0.3, label='r band')
            plt.errorbar(hjd_i, mag_i, yerr=emag_i, fmt='o', markersize=1, color='brown', alpha=0.3, label='i band')
            plt.xlabel('Time (HJD)')
            plt.ylabel('$Magnitude$')
            ax = plt.gca()
            ax.invert_yaxis()
            plt.grid()
            plt.legend(loc='best')
            #plt.xlim(8866.6,8867)
            plt.savefig('/mnt/c/Users/sshah/Documents/eswar_anohita_project/9_12/mstar_gt_2msun/lc_'+str(name1)+'_'+str(name2)+'_'+str(name3)+'.png')
            plt.clf()
