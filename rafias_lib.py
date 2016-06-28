from astropy.io import fits
from astropy.table import Table, Column, hstack
from photutils import CircularAperture, CircularAnnulus, aperture_photometry
from astropy.modeling import models, fitting
import numpy as np
import glob
import matplotlib
import matplotlib.pyplot as plt



""" 
    time_series prameters:
        center_x = x coordinate of the circular aperture, type = int
        center_y = y coordinate of the circular aperture, type = int
        radius = radius of circular aperture, type = int
        hdu_filenames = list of fits filenames (not the files themselves, just the names. i.e. strings), type = list
    About the function: 
        The function uses circular aperture to calculate the flux of each image provided to it. It then calculates the time from each
        image's metadata and creates a table column names: 'Flux' and 'Time' of length = number of files given to the function. 
                       
"""
def time_series(center_x,center_y, radius, hdu_filenames):
    center = (center_x,center_y)
    single_rad_data = Table(names=('Flux','Time'))
    for hdus in hdu_filenames:
        hdu = fits.open(hdus)
        image = hdu[0].data
        image2d = image[0]
        mask = np.isnan(image2d) == True
        aperture = CircularAperture(center, r = radius)
        phot_table = aperture_photometry(image2d, aperture, mask = mask)
        header = hdu[0].header
        time = [(header["NGROUP"] + 1) * header["TGROUP"] * (header["ON_NINT"] - 1)]
        a = [phot_table[0][0]]
        b = time
        single_rad_data.add_row([a,b])
        hdu.close()
    return single_rad_data





""" 
    light_curve prameters:
        x = x data of your plot. i.e. the time array, type = array
        y = y data of your plot. i.e. the flux array, type = array
        x_err = set of errors in the x direction, type = array
        y_err = set of errors in the y direction, type = array
        style = fmt. ie. the color and style of your curve, type = string
    About the function: 
        This function takes in flux and time data and plots them to create a light curve. All of the arrays should have the same dimension.
        The function also adds in error bars to every point. If there is no error, simply put x_err = 0 and y_err = 0.
"""
def light_curve(x, y, x_err, y_err, style):
    plt.errorbar(x, y/np.median(y), xerr = x_err, yerr = y_err, fmt= style)
    plt.xlabel('Time[sec]')
    plt.ylabel('Normalized Flux[DN/s]')
    
    
    
    

""" 
    light_curve prameters:
        x = x data of your plot. i.e. the time array, type = array
        y = y data of your plot. i.e. the flux array, type = array
        bin_size_low = lower limit of bin size, type = int
        bin_size_up = upper limit of bin size, type = int
        bin_size_inc = The increment by which bin size will be increasing, type = int
        num_points = Number of points in the data/ length of any of the data arrays
        style = fmt. ie. the color and style of your curve, type = string
    About the function: 
        This function takes in your data set, number of bins you want and their sizes to create a rms vs. bins plot.
"""
def rms_vs_bin(x, y, bin_size_low, bin_size_up, bin_size_inc, num_points, style):
    stdev_array = []
    time_array = []
    bin_size_array = np.arange(bin_size_low, bin_size_up, bin_size_inc)
    for bin_size in bin_size_array:
        flux_array = []
        for bins in range(0, num_points, bin_size):
            bin_start = bins
            bin_end   = bins + bin_size
            flux_in_one_bin = np.average(y[bin_start:bin_end])
            flux_array.append(flux_in_one_bin)
        norm_flux_array = flux_array/np.median(y[bin_start:bin_end])
        stdev_in_one_bin = np.std(norm_flux_array)
        stdev_array.append(stdev_in_one_bin)
        time_point = x[bin_size] - x[0]
        time_array.append(time_point)

    model = stdev_array[0]/np.sqrt(bin_size_array)       
    plt.loglog(time_array,stdev_array, style)
    plt.loglog(time_array, model, 'k--')
    plt.xlabel('Bin size')
    plt.ylabel('Standard Deviation (DN/s)')
    
    
    
    
def norm_flux_error(flux, gain, header):    
    errors_DNps = (np.sqrt(flux*header['INTTIME']*gain))/(gain*header['INTTIME'])
    errors_normalized = errors_DNps/flux
    return errors_normalized
  

    
    