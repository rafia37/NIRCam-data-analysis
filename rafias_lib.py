from astropy.io import fits
from astropy.table import Table, Column, hstack
from photutils import CircularAperture, CircularAnnulus, aperture_photometry
from astropy.modeling import models, fitting
import numpy as np
import pdb
import glob
import matplotlib
import matplotlib.pyplot as plt



def time_series(center_x,center_y, radius, hdu_filenames):
    
    """ 
    time_series prameters:
        center_x = x coordinate of the circular aperture, type = float
        center_y = y coordinate of the circular aperture, type = float
        radius = radius of circular aperture, type = float
        hdu_filenames = list of fits filenames (not the files themselves, just the names. i.e. strings), type = list
    About the function: 
        Returns a simple time series without background subtraction, fitted centers or other stuff
                       
    """
    
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





def light_curve(x, y, x_err, y_err, style):
    
    """ 
    light_curve prameters:
        x = x data of your plot. i.e. the time array, type = array
        y = y data of your plot. i.e. the flux array, type = array
        x_err = set of errors in the x direction, type = array
        y_err = set of errors in the y direction, type = array
        style = fmt. ie. the color and style of your curve, type = string
    About the function: 
        Produces a simple light curve with errorbars. The flux array will be normalized.
    """
    
    plt.errorbar(x, y/np.median(y), xerr = x_err, yerr = y_err, fmt= style)
    plt.xlabel('Time[sec]')
    plt.ylabel('Normalized Flux[DN/s]')
    
    
    
    

def rms_vs_bin(x, y, bin_size_low, bin_size_up, bin_size_inc, num_points, style):
    
    """ 
    rms_vs_bin prameters:
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
    
    
    
    

def norm_flux_error(flux, gain, hdu_filenames):
    
    """ 
    norm_flux_error prameters:
        flux = Data i.e. the flux array, type = array/table
        gain = detector's gain parameter, type = float
        hdu_filenames = list of fits filenames (not the files themselves, just the names. i.e. strings), type = list
    About the function: 
        This function takes in your flux data and calculates normalized error at each point and returns the error table/array. Expect the table length to be equal to array/table length of your data.
    """
    
    hdu = fits.open(hdu_filenames[2])
    header = hdu[0].header
    errors_DNps = (np.sqrt(flux*header['INTTIME']*gain))/(gain*header['INTTIME'])
    errors_normalized = errors_DNps/flux
    return errors_normalized
  
    




    

def gen_center_g2d(center_x, center_y, box_width, amp, x_std, y_std, Theta, hdu_filenames):
    
    """
    gen_center_g2d prameters:
        center_x = x coordinate of the circular aperture, type = float
        center_y = y coordinate of the circular aperture, type = float
        amp = amplitude of the gaussian.  Find from the projection curve along the center, type = float
        x_std = Standard deviation of the Gaussian in x before rotating by theta, type = float
        y_std = Standard deviation of the Gaussian in y before rotating by theta, type = float
        Theta = Rotation angle in radians. The rotation angle increases counterclockwise, type = float
        hdu_filenames = list of fits filenames (not the files themselves, just the names. i.e. strings), type = list
    About the function: 
        The function does a Gaussian fitting to your 2d data using your guess parameters and then uses the fitted value of those parameters to calculate centers of each image. It returns an array named 'separate_centers' that contain centers of each image. 
    """
    
    x_values = []
    y_values = []
    hdu = fits.open(hdu_filenames[2])
    image = hdu[0].data
    image2d = image[0]
    y_pos, x_pos = np.mgrid[:image2d.shape[0],:image2d.shape[1]]
    for hdus in hdu_filenames:
        hdu = fits.open(hdus)
        image = hdu[0].data
        image2d = image[0]
        fit_g = fitting.LevMarLSQFitter()
        gauss2D = models.Gaussian2D(amplitude = amp, x_mean = center_x, y_mean = center_y, x_stddev = x_std, y_stddev = y_std, theta = Theta)
        g = fit_g(gauss2D,x_pos[center_y-box_width:center_y+box_width,center_x-box_width:center_x+box_width],y_pos[center_y-box_width:center_y+box_width,center_x-box_width:center_x+box_width],image2d[center_y-box_width:center_y+box_width,center_x-box_width:center_x+box_width])
        g2 = fit_g(g,x_pos[center_y-box_width:center_y+box_width,center_x-box_width:center_x+box_width],y_pos[center_y-box_width:center_y+box_width,center_x-box_width:center_x+box_width],image2d[center_y-box_width:center_y+box_width,center_x-box_width:center_x+box_width])
        x_values.append(g2.x_mean)
        y_values.append(g2.y_mean)
        hdu.close()
    separate_centers = zip(x_values,y_values)
    return separate_centers








def radius_testing(centers, r_src_low, r_src_up, r_src_inc, r_in_low, r_in_up, r_in_inc, r_out_low, r_out_up, r_out_inc, hdu_filenames):
    """ 
    radius_testing prameters:
        centers = center of every image of your data set, type = array of tuples
        r_src_low = lower limit of source radius (circular aperture radius) array, type = float 
        r_src_up = upper limit of source radius (circular aperture radius) array, type = float
        r_src_inc = increment of source radius (circular aperture radius) array, type = float
        r_in parameters = Same as the r_src parameters, just for inner radius of annular apperture (this is for background subtraction)
        r_out parameters = Same as the r_src parameters, just for outer radius of annular apperture (this is for background subtraction)
        hdu_filenames = list of fits filenames (not the files themselves, just the names. i.e. strings), type = list
        About the function: 
        This function tests different radii to see which combination gives the least standard deviation. It prints the least standard deviation, the radii combination for which it occurs and return a table with table names: 'Median_Res_Flux','St_Dev', 'norm_stdev', 'r_source','r_in','r_out'.  
    """
    r_source = np.arange(r_src_low,r_src_up,r_src_inc)
    r_inner = np.arange(r_in_low,r_in_up,r_in_inc)
    r_outer = np.arange(r_out_low,r_out_up,r_out_inc)
    flux_and_parameters = Table(names=('residual_aperture_sum', 'r_source', 'r_in','r_out'))
    for index, hdus in enumerate(hdu_filenames):
        hdu = fits.open(hdus)
        image = hdu[0].data
        image2d = image[0,:,:]
        mask = np.isnan(image2d) == True
        for r in r_source:
            for r_in in r_inner:
                for r_out in r_outer:
                    if (r<r_in) and (r<r_out) and (r_in<r_out):
                        aperture = CircularAperture(centers[index], r)
                        annular_apperture =CircularAnnulus(centers[index], r_in, r_out)
                        rawflux_table = aperture_photometry(image2d, aperture, mask = mask)
                        bkgflux_table = aperture_photometry(image2d, annular_apperture, mask = mask)
                        phot_table = hstack([rawflux_table, bkgflux_table], table_names = ['raw','bkg'])
                        bkg_mean = phot_table['aperture_sum_bkg']/annular_apperture.area()
                        bkg_sum = bkg_mean*aperture.area()
                        final_sum = phot_table['aperture_sum_raw'] - bkg_sum
                        phot_table['residual_aperture_sum'] = final_sum
                        flux_and_parameters.add_row([final_sum,r,r_in,r_out])
    
    #Generating median flux and standard deviation at each r_source
    s1 = Table(names=('Median_Res_Flux','St_Dev', 'norm_stdev', 'r_source','r_in','r_out'))
    for r in r_source:
        for r_in in r_inner:
            for r_out in r_outer:
                if (r<r_in) and (r<r_out) and (r_in<r_out):
                    indices = ((flux_and_parameters['r_source'] == r) & (flux_and_parameters['r_in'] == r_in) & (flux_and_parameters['r_out'] == r_out)) 
                    st_dev = np.std(flux_and_parameters["residual_aperture_sum"][indices])
                    median_flux = np.median(flux_and_parameters["residual_aperture_sum"][indices])
                    norm_stdev = st_dev/median_flux
                    s1.add_row([median_flux,st_dev,norm_stdev,r,r_in,r_out])
    
    #Finding the best combination
    r1 = s1['r_source']
    r_in1 = s1['r_in']
    r_out1 = s1['r_out']
    min_std_dev = np.amin(s1['norm_stdev'])
    best_r = r1[np.argmin(s1['norm_stdev'])]
    best_r_in = r_in1[np.argmin(s1['norm_stdev'])]
    best_r_out = r_out1[np.argmin(s1['norm_stdev'])]
    print "The minimum Standard deviation is %f" % min_std_dev
    print "It occurs for the radius r = %f" % best_r
    print "It occurs for the inner radius r_in = %f" % best_r_in
    print "It occurs for the outer radius r_out = %f" % best_r_out
    return s1
    
    
    
    
    

    
    
def average_residual_flux(centers_a1, centers_b4, R, R_in, R_out, hdu_filenames, hdu_filenames_b4):
    
    """ 
    average_residual_flux prameters:
        center_a1 = center tuples of every image for 1st data set, type = array of tuples
        center_b4 = center tuples of every image for 2nd data set, type = array of tuples
        R = radius of circular aperture, type = float
        R_in = radius of inner annular aperture, type = float
        R_out = radius of outer annular aperture, type = float
        hdu_filenames = list of fits filenames (not the files themselves, just the names. i.e. strings), type = list
        hdu_filenames_b4 = list of fits filenames [2nd dataset] (not the files themselves, just the names. i.e. strings), type = list
    About the function: 
        The function uses circular aperture & annular aperture to calculate the average [raw flux, bkg flux, res flux] of each image from 2 different data sets provided to it. It then calculates the time from each image's metadata & averages them from botha data sets. It finally creates a table with column names: 'a1_b4_raw_flux', 'a1_b4_bkg_flux', 'a1_res_flux', 'b4_res_flux', 'a1_b4_res_flux' and 'Time' of length = number of files given to the function.                    
    """

    a1_b4_flux = Table(names=('a1_b4_raw_flux','a1_b4_bkg_flux','a1_res_flux','b4_res_flux','a1_b4_res_flux','Time'))
    for index, (hdus, hdus_b4) in enumerate(zip(hdu_filenames, hdu_filenames_b4)):
        # opening, slicing and creating mask for a1 data
        hdu = fits.open(hdus)
        header = hdu[0].header
        time = (header["NGROUP"] + 1) * header["TGROUP"] * (header["ON_NINT"] - 1)
        image = hdu[0].data
        image2d = image[0,:,:]
        mask_a1 = np.isnan(image2d) == True
        # opening, slicing and creating mask for b4 data
        hdu_b4 = fits.open(hdus_b4)
        header_b4 = hdu_b4[0].header
        time_b4 = (header_b4["NGROUP"] + 1) * header_b4["TGROUP"] * (header_b4["ON_NINT"] - 1)
        image_b4 = hdu_b4[0].data
        image2d_b4 = image_b4[0,:,:]
        mask_b4 = np.isnan(image2d_b4) == True
        # Defining circular & annular aperture for a1 data
        aperture_a1 = CircularAperture(centers_a1[index], r = R)
        annular_apperture_a1 =CircularAnnulus(centers_a1[index], r_in = R_in, r_out = R_out)
        # Defining circular & annular aperture for b4 data
        aperture_b4 = CircularAperture(centers_b4[index], r = R)
        annular_apperture_b4 =CircularAnnulus(centers_b4[index], r_in = R_in, r_out = R_out)
        # Photometric analysis of a1 data
        rawflux_table_a1 = aperture_photometry(image2d, aperture_a1, mask = mask_a1)
        bkgflux_table_a1 = aperture_photometry(image2d, annular_apperture_a1, mask = mask_a1)
        phot_table_a1 = hstack([rawflux_table_a1, bkgflux_table_a1], table_names = ['raw','bkg'])
        bkg_mean_a1 = phot_table_a1['aperture_sum_bkg']/annular_apperture_a1.area()
        bkg_sum_a1 = bkg_mean_a1*aperture_a1.area()
        final_sum_a1 = phot_table_a1['aperture_sum_raw'] - bkg_sum_a1
        phot_table_a1['residual_aperture_sum'] = final_sum_a1
        # Photometric analysis of b4 data
        rawflux_table_b4 = aperture_photometry(image2d_b4, aperture_b4, mask = mask_b4)
        bkgflux_table_b4 = aperture_photometry(image2d_b4, annular_apperture_b4, mask = mask_b4)
        phot_table_b4 = hstack([rawflux_table_b4, bkgflux_table_b4], table_names = ['raw','bkg'])
        bkg_mean_b4 = phot_table_b4['aperture_sum_bkg']/annular_apperture_b4.area()
        bkg_sum_b4 = bkg_mean_b4*aperture_b4.area()
        final_sum_b4 = phot_table_b4['aperture_sum_raw'] - bkg_sum_b4
        phot_table_b4['residual_aperture_sum'] = final_sum_b4
        # Fixing Table columns
        average_time = [(time + time_b4)/2] 
        a = (phot_table_a1[0][0] + phot_table_b4[0][0])/2
        b = (phot_table_a1[0][3] + phot_table_b4[0][3])/2
        c = phot_table_a1[0][6]
        d = phot_table_b4[0][6]
        e = (phot_table_a1[0][6] + phot_table_b4[0][6])/2
        f = average_time
        a1_b4_flux.add_row([a,b,c,d,e,f])
        hdu.close()
        hdu_b4.close()
    return a1_b4_flux









def linear_bestfit(x, y, x_err, y_err, slope_guess, intercept_guess, style):
    
    """ 
    linear_bestfit prameters:
        x = x data of your plot. i.e. the time array, type = array
        y = y data of your plot. i.e. the residual flux array, type = array
        x_err = set of errors in the x direction, type = array
        y_err = set of errors in the y direction, type = array
        slope_guess = guess slope of best fit line, type = float
        intercept_guess = guess intercept of bet fit line,  type = float
        style = fmt. ie. the color and style of your curve, type = string
    About the function: 
        This function does a linear best fit to your data distribution and returns the detrended flux data in the form of an array.
    """
    norm_y = y/np.median(y)
    l_init = models.Linear1D(slope = slope_guess, intercept = intercept_guess)
    fit_l = fitting.LevMarLSQFitter()
    l = fit_l(l_init, x, norm_y)
    detrend_flux_data = norm_y/l(x)
    # Plot the data with bets fit line
    plt.subplot(1,2,1)
    plt.errorbar(x, y/np.median(y), xerr = x_err, yerr = y_err, fmt= style)
    plt.plot(x, l(x), 'k--')
    plt.xlabel('Time[sec]')
    plt.ylabel('Normalized Flux')
    plt.title('Normalized Data With Linear Best Fit')
    plt.subplot(1,2,2)
    plt.plot(x, detrend_flux_data, '.-')
    plt.xlabel('Time[sec]')
    plt.ylabel('Normalized Detrended Flux')
    plt.title('Detrended Time Series')
    return detrend_flux_data
    