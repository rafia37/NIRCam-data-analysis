
# coding: utf-8

# In[1]:

#Import everything you need
from astropy.io import fits
from astropy.table import Table, Column
from photutils import CircularAperture, aperture_photometry
import numpy as np
import glob
import matplotlib
import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib inline')
plt.rcParams['figure.figsize'] = (10.0, 8.0)


# In[23]:

# Import all the files (A1 data)
hdu_filenames = glob.glob('//data1//tso_analysis//wlp8_sub_data//*.fits')
hdu_filenames


# Random comment

# In[24]:

# Run the code for the best aperture radius
center = (167,161)
single_rad_data = Table(names=('Flux','Time'))
for hdus in hdu_filenames:
    hdu = fits.open(hdus)
    image = hdu[0].data
    image2d = image[0]
    mask = np.isnan(image2d) == True
    aperture = CircularAperture(center, r = 39)
    phot_table = aperture_photometry(image2d, aperture, mask = mask)
    header = hdu[0].header
    time = [(header["NGROUP"] + 1) * header["TGROUP"] * (header["ON_NINT"] - 1)]
    a = [phot_table[0][0]]
    b = time
    single_rad_data.add_row([a,b])
    hdu.close()
single_rad_data


# In[9]:

plt.imshow(image2d)  #ignore this line


# In[25]:

# Light curve of first data set without error bars
x = single_rad_data['Time']
y = single_rad_data['Flux']
first_data_plot = plt.plot(x, y/np.median(y))
plt.xlabel('Time[sec]')
plt.ylabel('Flux[DN/s]')


# In[26]:

# Calculate error (first data set)
gain = 2.2
errors_DNps = (np.sqrt((single_rad_data['Flux'])*header['INTTIME']*gain))/(gain*header['INTTIME'])
errors_normalized = errors_DNps/single_rad_data['Flux']
errors_normalized


# In[19]:

# Light curve of first data set with error bars
x = single_rad_data['Time']
y = single_rad_data['Flux']
plt.errorbar(x, y/np.median(y), xerr = 0, yerr = errors_normalized, fmt='.')
plt.xlabel('Time[sec]')
plt.ylabel('Flux[DN/s]')

#Generating white noise
# sigma = errors.mean()
# med = np.median(y)
# white = np.random.normal(med,sigma,y.size)
# plt.plot(x,white, 'o')


# In[ ]:

#Calculate Aperture sum error with the 2nd image plane
center = (167,161)
t_2 = Table(names=('Flux','Time','Aperture_sum_err'))
for hdus in hdu_filenames:
    hdu = fits.open(hdus)
    image = hdu[0].data
    error2d = image[1]
    mask = np.isnan(image2d) == True
    aperture = CircularAperture(center, r = 39)
    phot_table = aperture_photometry(image2d, aperture, error = error2d, mask = mask)
    header = hdu[0].header
    time = [(header["NGROUP"] + 1) * header["TGROUP"] * (header["ON_NINT"] - 1)]
    a = [phot_table[0][0]]
    b = time
    c = [phot_table[0][1]]
    t_2.add_row([a,b,c])
    hdu.close()
t_2


# In[ ]:

# Compare error by previous method to this one
err_comparison_table = Table(names = ('Error by formula','Error by image plane'))
for index in range(306):
    err_comparison_table.add_row([errors_normalized[index],t_2['Aperture_sum_err'][index]])
err_comparison_table


# In[ ]:

plt.subplot(2,1,1)
x = single_rad_data['Time']
y = single_rad_data['Flux']
plt.errorbar(x, y/np.median(y), xerr = 0, yerr = errors_normalized, fmt='.')
plt.xlabel('Time[sec]')
plt.ylabel('Flux[DN/s]')
plt.title('Normalized error by formula')

plt.subplot(2,1,2)
x = single_rad_data['Time']
y = single_rad_data['Flux']
plt.errorbar(x, y, xerr = 0, yerr = t_2['Aperture_sum_err'], fmt='.')
plt.xlabel('Time[sec]')
plt.ylabel('Flux[DN/s]')
plt.title('Error from error image plane')


# In[7]:

# ** B4 data **


# In[18]:

#Importing second data set
hdu_filenames_b4 = glob.glob('//data//External//ISIMCV3_unzipped//NRCN821//fitsfilesonly//raw_separated_MMM//NRCN821WLP8SUB-6012134957_1_489_SE_2016-01-12T16h43m29//*slp.fits')
hdu_filenames_b4


# In[10]:

#Generate flux & time of each image for all aperture radii (second data set)
center = (161,163)
t_new = Table(names=('Flux','Time','Aperture Radius'))
for hdus_new in hdu_filenames_b4:
    hdu_new = fits.open(hdus_new)
    image_new = hdu_new[0].data
    image2d_new = image_new[0,:,:]
    mask = np.isnan(image2d_new) == True
    radius = np.arange(5,200,1)
    for r in radius:
        aperture = CircularAperture(center, r)
        phot_table = aperture_photometry(image2d_new, aperture, mask = mask)
        header_new = hdu_new[0].header
        time = [(header_new["NGROUP"] + 1) * header_new["TGROUP"] * (header_new["ON_NINT"] - 1)]
        a = [phot_table[0][0]]
        b = time
        c = r
        t_new.add_row([a,b,c])
    hdu_new.close()
t_new


# In[11]:

# Second data set
# Generate Standard deviation of all images for each aperture radii
# Generate Median flux of all images for each aperture radii (so that you can normalize standard deviation)
s_new = Table(names=('Median_Flux','St_Dev', 'Ap_Rad'))
for r in radius:
    indices = t_new['Aperture Radius'] == r
    st_dev = np.std(t_new["Flux"][indices])
    median_flux = np.median(t_new["Flux"][indices])
    s_new.add_row([median_flux,st_dev,r])
s_new


# In[12]:

# Plot normalized standard deviation (second data set)
u_new = s_new['Ap_Rad']
v_new = s_new['St_Dev']/s_new['Median_Flux']
plt.plot(u_new, v_new)
plt.ylabel('Standard Deviation (DN/S)')
plt.xlabel('Aperture Radius (pixels)')


# In[13]:

# Find minimum deviation and the corresponding aperture radius (second data set)
min_stdev = np.amin(v_new)
best_ap_rad = u_new[v_new.argmin()]
print "The minimum standard deviation is %f" % min_stdev
print "It occurs for the radius r = %f" % best_ap_rad


# In[14]:

plt.imshow(image2d_new) #Second data set 


# In[15]:

# Figuring out the center for image2d_new
plt.subplot(2,1,1)
plt.title("Flux Vs. x")
plt.plot(image2d_new[:,163],'r-')
plt.grid(True)
plt.subplot(2,1,2)
plt.title("Flux Vs. y")
plt.plot(image2d_new[161,:],'r-')
plt.grid(True)


# In[19]:

#Flux Vs. time B4 data
center = (161,163)
single_rad_data_b4 = Table(names=('Flux','Time'))
for hdus_b4 in hdu_filenames_b4:
    hdu_b4 = fits.open(hdus_b4)
    image_b4 = hdu_b4[0].data
    image2d_b4 = image_b4[0,:,:]
    mask = np.isnan(image2d_b4) == True
    aperture = CircularAperture(center, r = 28)
    phot_table = aperture_photometry(image2d_b4, aperture, mask = mask)
    header_b4 = hdu_b4[0].header
    time = [(header_b4["NGROUP"] + 1) * header_b4["TGROUP"] * (header_b4["ON_NINT"] - 1)]
    a = [phot_table[0][0]]
    b = time
    single_rad_data_b4.add_row([a,b])
    hdu_b4.close()
single_rad_data_b4


# In[20]:

# Creating light curve
x_b4 = single_rad_data_b4['Time']
y_b4 = single_rad_data_b4['Flux']
second_data_plot = plt.plot(x_b4, y_b4/np.median(y_b4), 'r-')
plt.xlabel('Time[sec]')
plt.ylabel('Flux[DN/s]')


# In[21]:

# ** combining both data sets **


# In[45]:

# Averaging A1 data with b4
a1_flux = single_rad_data['Flux']
b4_flux = single_rad_data_b4['Flux'][:306]
average_flux_array = (a1_flux + b4_flux)/2
plt.plot(average_time_array,average_flux_array/np.median(average_flux_array))
np.std(average_flux_array)/np.median(average_flux_array)


# In[ ]:

# stdev of the average flux
center = (161,163)
a1_time = single_rad_data['Time']
b4_time = single_rad_data_b4['Time'][:306]
average_time_array = (a1_time + b4_time)/2 
image_average = (image2d + image2d_b4)/2
t_combo = Table(names=('Flux','Time','Aperture Radius'))
for index, hdus in enumerate(hdu_filenames):
    hdus_b4 = hdu_filenames_b4[index]
    hdu = fits.open(hdus)
    hdu_b4 = fits.open(hdus_b4)
    image = hdu[0].data
    image_b4 = hdu_b4[0].data
    image2d = image[0,:,:]
    image2d_b4 = image_b4[0,:,:]
    mask_a1 = np.isnan(image2d) == True
    mask_b4 = np.isnan(image2d_b4) == True
    radius = np.arange(5,200,1)
    for r in radius:
        aperture = CircularAperture(center, r)
        phot_table_a1 = aperture_photometry(image2d, aperture, mask = mask_a1)
        phot_table_b4 = aperture_photometry(image2d_b4, aperture, mask = mask_b4)
        a = (phot_table_a1[0][0] + phot_table_b4[0][0])/2
        b = average_time_array[index]
        c = r
        t_combo.add_row([a,b,c])
    hdu.close()
    hdu_b4.close()
t_combo


# In[ ]:

# Generate Standard deviation of all images for each aperture radii
# Generate Median flux of all images for each aperture radii (so that you can normalize standard deviation)
stdev_table = Table(names=('Median_Flux','St_Dev', 'Ap_Rad'))
for r in radius:
    indices = t_combo['Aperture Radius'] == r
    st_dev = np.std(t_combo["Flux"][indices])
    median_flux = np.median(t_combo["Flux"][indices])
    s_new.add_row([median_flux,st_dev,r])
stdev_table


# In[ ]:

x1 = stdev_table['Ap_Rad']
y1 = stdev_table['St_Dev']/stdev_table['Median_Flux']
plt.plot(x1, y1)
plt.title('StDev Vs. ApRad for the average of A1 & b4')
plt.ylabel('Standard Deviation (DN/S)')
plt.xlabel('Aperture Radius (pixels)')


# In[ ]:

# Find minimum deviation and the corresponding aperture radius (average)
min_stdev = np.amin(x1)
best_ap_rad = x1[y1.argmin()]
print "The minimum standard deviation is %f" % min_stdev
print "It occurs for the radius r = %f" % best_ap_rad


# In[ ]:



