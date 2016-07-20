
# coding: utf-8
## An example of producing light curve with errorbars using my library
data = rl.time_series(167, 161, 80, hdu_filenames)
errors = rl.norm_flux_error(data['Flux'], 2.2, hdu_filenames)
rl.light_curve(data['Time'], data['Flux'], 0, errors, '.-')
# In[66]:

from astropy.io import fits
from astropy.table import Table, Column, hstack
import numpy as np
import rafias_lib as rl
import glob
import matplotlib
import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib inline')
plt.rcParams['figure.figsize'] = (10.0, 8.0)


# In[49]:

the_type = ['MMM', 'MMP', 'MPM', 'PMM', 'PPM', 'PMP', 'MPP', 'PPP']
the_nmbr = ['481', '489']


# In[55]:

a1_files = []
for each in the_type:
    a1_files.append(glob.glob('/data/External/ISIMCV3_unzipped/NRCN821/fitsfilesonly/raw_separated_'+each+
                                 '/NRCN821WLP8SUB-6012134600_1_'+the_nmbr[0]+'_SE_2016-01-12T16h42m53/*.slp.fits'))
print len(a1_files), len(a1_files[1])


# In[63]:

b4_files = []
for each in the_type:
    b4_files.append(glob.glob('/data/External/ISIMCV3_unzipped/NRCN821/fitsfilesonly/raw_separated_'+each+
                                 '/NRCN821WLP8SUB-6012134957_1_'+the_nmbr[1]+'_SE_2016-01-12T16h43m29/*.slp.fits'))
print len(b4_files), len(b4_files[0])


# In[33]:

def get_stdev(hdu_filenames_s1, hdu_filenames_s2):
    data = rl.average_residual_flux(167, 161, 80, 90, 100, hdu_filenames_s1, hdu_filenames_s2)
    norm_res_flux = data['a1_b4_res_flux']/np.median(data['a1_b4_res_flux'])
    return np.std(norm_res_flux)


# In[71]:

# I tried to get std_devs seperately, they still give me the same answer
print get_stdev(a1_files[1],b4_files[1])


# In[70]:

# getting std_devs with for loop
stdev_comparison = Table(names = ('Index','Stdev'))
for index, (a1, b4) in enumerate(zip(a1_files, b4_files)):
    stdev = get_stdev(a1, b4)
    stdev_comparison.add_row([index,stdev])
stdev_comparison


# In[ ]:



