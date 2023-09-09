import os
import sys
import numpy as np
import requests

p = np.genfromtxt('lum2.txt')

ra = p[:,0]
dec = p[:,1]

good_points = np.where((ra!=np.nan) & (dec!=np.nan))[0]
ra_good = ra[good_points]
dec_good = dec[good_points]

rng = np.linspace(40076, len(ra_good),(len(ra_good) - 40076), dtype=np.int64)
print(rng, len(ra_good))
import requests

for i in range(len(ra_good)):
    print('i=', i, rng[i])
    b=rng[i]
    #print(b,i)
    #b = i
    #url = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph-light_curves?POS=CIRCLE"+"%20"+str(ra_good[i])+"%20"+str(dec_good[i])+"%20"+"0.00028&BANDNAME=g&NOBS_MIN=5&BAD_CATFLAGS_MASK=32768&FORMAT=CSV"
    #url = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph-lightcurves?POS=CIRCLE%20"+str(ra_good[i])+"%20"+str(dec_good[i])+"%20" + "0.00028&BANDNAME=g&NOBS_MIN=3&TIME=-Inf+Inf&BAD_CATFLAGS_MASK=32768&FORTMAT=CSV"
    #print(url)
    #file_i = wget.download(url)
    url = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?POS=CIRCLE%20"+str(ra_good[b])+"%20"+str(dec_good[b])+"%200.00028&BANDNAME=g&TIME=-Inf+Inf&NOBS_MIN=3&BAD_CATFLAGS_MASK=32768&FORMAT=CSV"
    print('url=', url)
    response = requests.get(url)
    print('response=', response.text)
    open("/mnt/c/Users/sshah/Documents/eswar_anohita_project/9_12/mstar_gt_2msun/ztf_data/g"+str(b)+".csv", "wb").write(response.content)

    url = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?POS=CIRCLE%20"+str(ra_good[b])+"%20"+str(dec_good[b])+"%200.00028&BANDNAME=r&TIME=-Inf+Inf&NOBS_MIN=3&BAD_CATFLAGS_MASK=32768&FORMAT=CSV"
    print('url=', url)
    response = requests.get(url)
    print('response=', response.text)
    open("/mnt/c/Users/sshah/Documents/eswar_anohita_project/9_12/mstar_gt_2msun/ztf_data/r"+str(b)+".csv", "wb").write(response.content)

    url = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?POS=CIRCLE%20"+str(ra_good[b])+"%20"+str(dec_good[b])+"%200.00028&BANDNAME=i&TIME=-Inf+Inf&NOBS_MIN=3&BAD_CATFLAGS_MASK=32768&FORMAT=CSV"
    print('url=', url)
    response = requests.get(url)
    print('response=', response.text)
    open("/mnt/c/Users/sshah/Documents/eswar_anohita_project/9_12/mstar_gt_2msun/ztf_data/i"+str(b)+".csv", "wb").write(response.content)
