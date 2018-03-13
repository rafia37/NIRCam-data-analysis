from astropy.io import fits
from astropy.table import Table, Column, hstack,vstack
from photutils import CircularAperture, RectangularAperture, CircularAnnulus, aperture_photometry
from astropy.modeling import models, fitting
import numpy as np
import everett_code as ec
import pdb, glob, itertools, warnings
import matplotlib.pyplot as plt







def fname_generator(tests, div, mmm = 'MMM'):
    """
    Takes:
    tests = list of test names as included in filenames. Type = list (of strings)
    div = the index of first test in tests that requires red files. Type = Integer
    
    Returns
    fname = Filename of NIRCam images for each test given in the tests list (for weak lens tests)
    """
    fn1 = '/data1/tso_analysis/all_tso_cv3/raw_separated_'+mmm+'/NRCN821WLP8'
    a1 = ['-*_1_481_SE_*/*.slp.fits', '-*_1_481_SE_*/*.red.fits'] # 0 for slp, 1 for red
    b4 = ['-*_1_489_SE_*/*.slp.fits', '-*_1_489_SE_*/*.red.fits']

    fname = []
    for i, t in enumerate(tests):                        
        if i < div:                                          #slp files
            globals()['%s_a1' % t] = np.sort(glob.glob(fn1 + t + a1[0]))[2:]
            globals()['%s_b4' % t] = np.sort(glob.glob(fn1 + t + b4[0]))[2:]
            
            if i==0: globals()['%s_b4' % t] =  globals()['%s_b4' % t][:-1] # for the unequal sub array issue
            
            fname.append(globals()['%s_a1' % t])
            fname.append(globals()['%s_b4' % t])
        else:                                              #red files
            globals()['%s_a1' % t] = np.sort(glob.glob(fn1 + t + a1[1]))[2:]
            globals()['%s_b4' % t] = np.sort(glob.glob(fn1 + t + b4[1]))[2:]
            fname.append(globals()['%s_a1' % t])
            fname.append(globals()['%s_b4' % t])
    
    return fname







def test_image(filename, r = False, r2 = False, f_name = False, Time = True):
    """
    Takes:
    filename = filename of a single image. Type = String
    r = red file, slope1 method. Type = Boolean
    r2 = red file, slope2 method. Type = Boolean
    f_name = name of flatfield file if it needs to be aplied. Type = String
    Time = whether you want the time or not. Type = Boolean
    
    Returns:
    image2d = 2 dimensional image array. Type = ndarray
    time = Integration time of the image. Type = Float
    header = header info of the fits file. Type = ndarray
    mask = mask that blocks out NaNs. Type = ndarray
    """
    hdu = fits.open(filename)
    image = hdu[0].data
    header = hdu[0].header
    hdu.close()
    
    if f_name != False:
        flat_file = fits.open(f_name)
        flat = flat_file[1].data
        flat_file.close()
        
    if r == False:          #.slp files
        image2d = image[0]
    elif r2 == False:       #.red file, Slope1 method
        if f_name != False:
            slope = (image[-1] - image[0])/((header['NGROUP']-1)*header['TFRAME']*header['NFRAME'])
            image2d = slope/flat
        else:
            image2d = (image[-1] - image[0])/((header['NGROUP']-1)*header['TFRAME']*header['NFRAME'])
    else:                     #.red file, Slope2 method
        if f_name != False:
            slope = image[-1]/(header['NGROUP']*header['TFRAME']*header['NFRAME'])
            image2d = slope/flat
        else:
            image2d = image[-1]/(header['NGROUP']*header['TFRAME']*header['NFRAME'])
    mask = np.isnan(image2d) == True
    if Time == True:
        time = [(header["NGROUP"] + 1) * header['TGROUP']* (header["ON_NINT"] - 1)]
    else:
        time = 0.0
    
    return image2d, time, header, mask








def photometry(image2d, cen_x, cen_y, mask, index = 0, shape = 'Circ', rad = None, r_in = None, r_out = None, ht = None, wid = None, w_in = None, w_out = None, h_out = None, ang = 0.0):
    """
    Takes:
    image2d = 2 dimensional image array. Type = ndarray
    cen_x, cen_y = x & y center position. Type = ndarray/list
    mask = mask that blocks out NaNs. Type = ndarray
    index = if cen_x and cen_Y is a list of more than 1 element, specify the desired index. Type = Integer
    shape = 'Circ':CircularAperture, 'Rect':RectangularAperture, 'CircAnn':CircularAnnulus, 'RectAnn':RectangularAnnulus
    rad, r_in, r_out, ht, wid, w_in, w_out, h_out, ang = Astropy's aperture parameters 
    
    Returns:
    flux = flux of the image extracted by the aperture described in the "shape" parameter. Type = Float
    aperture = aperture object created by astropy
    """
    if shape == 'Circ':
        aperture = CircularAperture((cen_x[index], cen_y[index]), r = rad)
    elif shape == 'Rect':
        aperture = RectangularAperture((cen_x[index], cen_y[index]), w = wid, h = ht, theta = ang)
    elif shape == 'CircAnn':
        aperture = CircularAnnulus((cen_x[index], cen_y[index]), r_in = r_in, r_out = r_out)
    elif shape == 'RectAnn':
        aperture = RectangularAnnulus((cen_x[index], cen_y[index]), w_in = w_in, w_out = w_out, h_out = h_out, theta = ang)
            
    phot_table = aperture_photometry(image2d, aperture, mask = mask)
    flux = phot_table[0][0]
    return flux, aperture









def col_col(image, mask1, cenX, cenY, r, box = 150, plots = False, col = 160):
    new_image = np.zeros_like(image)
    xcen = int(cenX)
    ycen = int(cenY)
    y, x = np.mgrid[:image.shape[0], :image.shape[1]]
    mask2 = ((x - xcen)**2 + (y - ycen)**2) <= (r**2)
    mask3 = (y > (ycen-box)) & (y < (ycen+box))
    for i in range(xcen-box, xcen+box, 1):
        final_pts = (mask1==False) & (mask2==False) & (mask3==True) & (x==i)
        flux = image[final_pts]
        y_coord = y[final_pts]
        original = image[:,i]
        a, b = ec.robust_poly(y_coord, flux, 1, sigreject=3.0, iteration=3)
        gen_y = np.arange(0, image.shape[1])
        fit = (a*gen_y) + b
        final = original - fit
        new_image[:, i] = final
        if (plots == True) & (i == col):
            plt.plot(y_coord, flux,'r.', label = "masked col %i" % i)
            plt.plot(gen_y, fit, 'k--', label = "polynomial fit")
            plt.plot(gen_y, original, 'b', label = "unmasked flux")
            plt.plot(gen_y, final, 'g', label = "detrended")
            plt.xlabel("y coordinate")
            plt.ylabel("flux")
            plt.legend(loc="best")
    return new_image







def row_row(image, mask1, cenX, cenY, r, box = 150, plots = False, row = 160):
    new_image = np.zeros_like(image)
    xcen = int(cenX)
    ycen = int(cenY)
    y, x = np.mgrid[:image.shape[0], :image.shape[1]]
    mask2 = ((x - xcen)**2 + (y - ycen)**2) <= (r**2)
    mask3 = (x > (xcen-box)) & (x < (xcen+box))
    for i in range(ycen-box, ycen+box, 1):
        final_pts = (mask1==False) & (mask2==False) & (mask3==True) & (y==i)
        flux = image[final_pts]
        x_coord = x[final_pts]
        original = image[i,:]
        a, b = ec.robust_poly(x_coord, flux, 1, sigreject=3.0, iteration=3)
        gen_x = np.arange(0, image.shape[0])
        fit = (a*gen_x) + b
        final = original - fit
        new_image[i, :] = final
        if (plots == True) & (i == row):
            plt.plot(x_coord, flux,'r.', label = "masked row %i" % i)
            plt.plot(gen_x, fit, 'k--', label = "polynomial fit")
            plt.plot(gen_x, original, 'b', label = "unmasked row %i" % i)
            plt.plot(gen_x, final, 'g', label = "detrended")
            plt.xlabel("x coordinate")
            plt.ylabel("flux")
            plt.legend(loc="best")
    return new_image
                
            

        

        





def time_series(xcenter, ycenter, filenames, r = None, r_in = None, r_out = None, rs_in = None, rs_out = None, flat_name = False, w = None, h = None, w_in = None, w_out = None, h_out = None, ws_in = None, ws_out = None, hs_out = None, red = False, red2 = False, bg_xcen = None, bg_ycen = None, mode = "astropy", src_shape = "Circ", bkg_shape = "Circ", average = "med"):

    flux_table = Table(names = ('raw_flux', 'bkg_flux', 'res_flux', 'time'))
    
    for i, hdu in enumerate(filenames):
        test_im = test_image(filename = hdu, r = red, r2 = red2, f_name = flat_name)
        image2d, time, header, mask = test_im[0], test_im[1], test_im[2], test_im[3]
        ap_phot = photometry(image2d, xcenter, ycenter, mask, index = i, shape = src_shape, rad = r, r_in = rs_in, r_out = rs_out, ht = h, wid = w, w_in = ws_in, w_out = ws_out, h_out = hs_out)
        raw_flux = ap_phot[0]
        source_ap = ap_phot[1]
       
        
        if mode == "astropy":
            if bkg_shape == "Circ":
                bkg_ap = CircularAnnulus((xcenter[i], ycenter[i]), r_in = r_in, r_out = r_out)
                bkg = aperture_photometry(image2d, bkg_ap, mask = mask)
                bkg_mean = bkg['aperture_sum']/bkg_ap.area()
            
            elif bkg_shape == "Rect":
                bkg_ap = photometry(image2d, bg_xcen, bg_ycen, mask, index = i, shape = bkg_shape, ht = h, wid = w)[1]
                bkg = aperture_photometry(image2d, bkg_ap, mask = mask)
                bkg_mean = bkg['aperture_sum']/bkg_ap.area()
            
            elif bkg_shape == "RectAnn":
                bkg_ap = RectangularAnnulus((xcenter[i], ycenter[i]), w_in = w_in, w_out = w_out, h_out = h_out, 
                                                   theta = 0.0)
                bkg = aperture_photometry(image2d, bkg_ap, mask = mask)
                bkg_mean = bkg['aperture_sum']/bkg_ap.area()
                
            else:
                warnings.warn("Not a recognized astropy shape")
            
            bkg_flux = bkg_mean*source_ap.area()
            res_flux = raw_flux - bkg_flux
                

        elif mode == "shapes":
            y, x = np.mgrid[:image2d.shape[0], :image2d.shape[1]]
            if bkg_shape == "Circ":
                bkg_pts = ((((x - xcenter[i])**2 + (y - ycenter[i])**2) > (r_in)**2) & 
                               (((x - xcenter[i])**2 + (y -ycenter[i])**2) < (r_out)**2))
            elif bkg_shape == "CIS":
                bkg_pts = ((((x - xcenter[i])**2 + (y - ycenter[i])**2) > (r_in)**2) & 
                               ((np.abs(x - xcenter[i]) < r_out) & (np.abs(y -ycenter[i]) < r_out)))
            else:
                warnings.warn("Not a recognized shape")

            if average == "med":
                bkg_med = np.nanmedian(image2d[bkg_pts])
            elif average == "avg":
                bkg_med = np.nanmean(image2d[bkg_pts])
            elif average =="mad":
                ad = np.abs(image2d[bkg_pts]-np.nanmedian(image2d[bkg_pts]))
                mad = np.nanmedian(ad)
                keep_pts = (np.abs(image2d-np.nanmedian(image2d[bkg_pts]))<(5*mad)) & bkg_pts
                bkg_med = np.nanmean(image2d[keep_pts])
            else:
                warnings.warn("Not a recognized average")

            bkg_flux = bkg_med*(np.pi*(r**2))
            res_flux = raw_flux - bkg_flux
            
        elif mode == "col_col":
            new_im = col_col(image2d, mask, xcenter[i], ycenter[i], r, box = 150)
            bkg_flux = 0
            res_flux = photometry(new_im, xcenter, ycenter, mask, index = i, shape = 'Circ', rad = r)[0]
            
        elif mode == "row_row":
            new_im = row_row(image2d, mask, xcenter[i], ycenter[i], r, box = 150)
            bkg_flux = 0
            res_flux = photometry(new_im, xcenter, ycenter, mask, index = i, shape = 'Circ', rad = r)[0]
            
        else:
            raise Warning("Not a recognized mode")
        
        flux_table.add_row([raw_flux, bkg_flux, res_flux, time])

    return flux_table











def light_curve(x, y, x_err = None, y_err = None, style = 'r.-', lbl = None):
    
    """ 
    PARAMETERS:
        x = x data of your plot. i.e. the time array; Type = Array
        y = y data of your plot. i.e. the flux array; Type = Array
        x_err = set of errors in the x direction. Default value = "None"; Type = Array
        y_err = set of errors in the y direction. Default value = "None"; Type = Array
        style = fmt. ie. the color and style of your curve. Default value = "r.-"; Type = String
        lbl = Label for the plot. Default value = "None"; Type = String
        
    RETURNS: 
        Plot = Simple light curve
    """
    
    plt.errorbar(x, y/np.median(y), xerr = x_err, yerr = y_err, fmt = style, label = lbl)
    plt.xlabel('Time[sec]')
    plt.ylabel('Normalized Flux[DN/s]')
    plt.title('Simple light curve')
    
    
    
    
    
    
    
    

def rms_vs_bin(x, y, bin_size_low, bin_size_up, bin_size_inc, num_points, style, lbl = None):
    
    """ 
    PARAMETERS:
        x = x data of your plot. i.e. the time array; Type = Array
        y = y data of your plot. i.e. the flux array; Type = Array
        bin_size_low = lower limit of bin size; Type = Int
        bin_size_up = upper limit of bin size; Type = Int
        bin_size_inc = The increment by which bin size will be increasing; Type = Int
        num_points = Number of points in the data/length of any of the data arrays; Type = Int
        style = fmt. ie. the color and style of your curve; Type = String
        lbl = Label for the plot. Default value = "None"; Type = String
        
    RETURNS: 
        Plot = rms vs. bin sive with ideal noise
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
        stdev_array.append(stdev_in_one_bin*1e6)
        time_point = x[bin_size] - x[0]
        time_array.append(time_point)

    model = stdev_array[0]/np.sqrt(bin_size_array)       
    plt.loglog(time_array,stdev_array, style, label = lbl)
    plt.loglog(time_array, model, 'k--')
    plt.xlabel('Bin size (seconds)')
    plt.ylabel('$\sigma$ (ppm)')
    
    
    
    
    
    
    
    

def norm_flux_error(flux, gain, hdu_filenames, red = False, red2 = False):
    
    """ 
    PARAMETERS:
        flux = Data i.e. the flux array, type = list/array/table
        gain = detector's gain parameter, type = float
        hdu_filenames = list of fits filenames, type = list [of strings]
        red = Whether the files are .red files or not. Default value: "False"; type = Boolean
        red2 = Whether you want to use Slope2 method or not. Default value: "False"; type = Boolean. 

    RETURNS: 
        norm_error = Normalized errors for each flux data, type = List
    """
    norm_error = []
    for i, hdus in enumerate(hdu_filenames):
        hdu = fits.open(hdus)
        header = hdu[0].header
        
        if red == False:
            errors_DNps = (np.sqrt(flux[i]*header['INTTIME']*gain))/(gain*header['INTTIME'])
        elif red2 == False:
            errors_DNps = (np.sqrt(flux[i]*((header['NGROUP']-1)*header['TGROUP'])*gain))/(((header['NGROUP']-1)*header['TGROUP'])*gain)
        else:
            errors_DNps = (np.sqrt(flux[i]*(header['NGROUP']*header['TGROUP'])*gain))/((header['NGROUP']*header['TGROUP'])*gain)
        
        errors_normalized = errors_DNps/flux[i]
        norm_error.append(errors_normalized)
        hdu.close()
    return norm_error








    

def gen_center_g2d(center_x, center_y, box_width, amp, x_std, y_std, Theta, hdu_filenames, red = False, red2 = False, flat_name = False):
    
    """
    PARAMETERS:
        center_x = x coordinate of the circular aperture; Type = float
        center_y = y coordinate of the circular aperture; Type = float
        amp = amplitude of the gaussian.  Find from the projection curve along the center; Type = float
        x_std = Standard deviation of the Gaussian in x before rotating by theta; Type = float
        y_std = Standard deviation of the Gaussian in y before rotating by theta; Type = float
        Theta = Rotation angle in radians. The rotation angle increases counterclockwise; Type = float
        hdu_filenames = list of fits filenames; Type = List [of strings]
        red = Whether the files are .red files or not. Default value: "False"; Type = Boolean
        red2 = Whether you want to use Slope2 method or not. Default value: "False"; Type = Boolean. 
    
    RETURNS:
        seperate_centers =  Center of each image; Type = Array [of tuples]
        x_values = x_value of center of each image; Type = Array
        y_values = y_value of center of each image; Type = Array
    """
    x_values = []
    y_values = []        
        
    #Fitting a gaussian model to each image in image2d list and returning center    
    for index, hdus in enumerate(hdu_filenames):
        image2d = test_image(filename = hdus, r = red, r2 = red2, f_name = flat_name)[0]
        y_pos, x_pos = np.mgrid[:image2d.shape[0],:image2d.shape[1]]
        fit_g = fitting.LevMarLSQFitter()
        gauss2D = models.Gaussian2D(amplitude = amp, x_mean = center_x, y_mean = center_y, x_stddev = x_std, y_stddev = y_std, theta = Theta)
        g = fit_g(gauss2D,x_pos[center_y-box_width:center_y+box_width,center_x-box_width:center_x+box_width],y_pos[center_y-box_width:center_y+box_width,center_x-box_width:center_x+box_width],image2d[center_y-box_width:center_y+box_width,center_x-box_width:center_x+box_width])
        g1 = fit_g(g,x_pos[center_y-box_width:center_y+box_width,center_x-box_width:center_x+box_width],y_pos[center_y-box_width:center_y+box_width,center_x-box_width:center_x+box_width],image2d[center_y-box_width:center_y+box_width,center_x-box_width:center_x+box_width])
        x_values.append(g1.x_mean[0])
        y_values.append(g1.y_mean[0])
    
    #Results
    separate_centers = zip(x_values,y_values) 
    return np.array(separate_centers), np.array(x_values), np.array(y_values)








    
    



def linear_bestfit(x, y, slope_guess, intercept_guess, show_plot = False, x_err = None, y_err = None, style = None):
    
    """ 
    PARAMETERS:
        x = x data of your plot. i.e. the time array; Type = Array
        y = y data of your plot. i.e. the residual flux array; Type = Array
        slope_guess = guess slope of best fit line; Type = Float
        intercept_guess = guess intercept of best fit line; Type = Float
        show_plot = If you want a plot of the time_series. Default value: "False"; Type = Boolean
        x_err = set of errors in the x direction. Default value: "None"; Type = Array
        y_err = set of errors in the y direction. Default value: "None"; Type = Array
        style = fmt. ie. the color and style of your curve. Default value: "None"; Type = String
        
    RETURNS: 
        detrended_flux_data = Flux data after applying linear model. Type = Array
        If show_plot = true, then:
        plot1 = Normalized Data With Linear Best Fit
        plot2 = Detrended Time Series
    """
    if show_plot == False:
        norm_y = y/np.median(y)
        l_init = models.Linear1D(slope = slope_guess, intercept = intercept_guess)
        fit_l = fitting.LevMarLSQFitter()
        l = fit_l(l_init, x, norm_y)
        detrend_flux_data = norm_y/l(x)
        return detrend_flux_data
    else:
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

    
    
    
    
    
    
    
    
    
    
    
def radius_testing(cenX1, cenY1, fnames1, cenX2, cenY2, fnames2, r_src_low, r_src_up, r_src_inc, r_in_low, r_in_up, r_in_inc, r_out_low, r_out_up, r_out_inc, Red = False, Red2 = False):
    
    r_source = np.arange(r_src_low,r_src_up,r_src_inc)
    r_inner = np.arange(r_in_low,r_in_up,r_in_inc)
    r_outer = np.arange(r_out_low,r_out_up,r_out_inc)
    rad_test = Table(names=('norm_stdev', 'r_source', 'r_in','r_out', 'rIn - r', 'rOut - rIn'))
    for R in r_source:
        for R_in in r_inner:
            for R_out in r_outer:
                if (R<R_in) and (R<R_out) and (R_in<R_out):
                    data1 = time_series(cenX1, cenY1, fnames1, r = R, r_in = R_in, r_out = R_out, red = Red, red2 = Red2)
                    data2 = time_series(cenX2, cenY2, fnames2, r = R, r_in = R_in, r_out = R_out, red = Red, red2 = Red2)
                    detrended1 = linear_bestfit(data1['time'], data1['res_flux'], 0.00002, 1)
                    detrended2 = linear_bestfit(data2['time'], data2['res_flux'], 0.00002, 1)
                    av = (detrended1 + detrended2)/2.
                    stdev = np.std(av)/np.median(av)
                    rad_test.add_row([stdev, R, R_in, R_out, R_in-R, R_out-R_in])
    
    #Finding the best combination
    a = rad_test['norm_stdev']
    min_std_dev = np.amin(a)
    best_r = rad_test['r_source'][np.argmin(a)]
    best_r_in = rad_test['r_in'][np.argmin(a)]
    best_r_out = rad_test['r_out'][np.argmin(a)]
    print "The minimum Standard deviation is %f" % min_std_dev
    print "It occurs for the radius r = %f" % best_r
    print "It occurs for the inner radius r_in = %f" % best_r_in
    print "It occurs for the outer radius r_out = %f" % best_r_out
    return rad_test

