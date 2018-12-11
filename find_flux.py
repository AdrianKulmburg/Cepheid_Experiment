from numpy import *
from general_fit import *
from astropy.io import fits
from photutils import find_peaks
from matplotlib.pyplot import *

threshold_star = 7000.0
emergency_threshold = 3000.0

initial_lower_threshold = 20000
saturation_threshold = 38000

def single_gaussian(parameters, data):
    A = parameters[0]
    
    mu_x = parameters[1]
    mu_y = parameters[2]
    
    a = parameters[3]
    b = parameters[4]
    c = parameters[5]

    D = 0.0
    
    x = data[0]
    y = data[1]
    
    return D + A / (2 * pi * sqrt(abs(a*c - b**2))) * exp(-1.0/2.0 * (c*(x - mu_x)**2 - 2*b*(x - mu_x)*(y - mu_y) + a*(y - mu_y)**2) / (a*c - b**2))

def find_peak_base_box(data, errors, peak_x, peak_y):
    """Finds the base of a given peak, where we assume
    that the threshold is zero.
    Note that for this, we assume that the peak is in
    middle of the picture, i.e. there should be no way
    to get an out of bounds error. Furthermore, we also
    assume that the source is separated enough from other
    sources for them not to impact our result."""

    box_radius = 1
    while 1:
    	first_column = data[peak_x - box_radius:peak_x + box_radius + 1, peak_y - box_radius]
    	first_row = data[peak_x - box_radius:peak_x + box_radius + 1, peak_y + box_radius]
    	last_column = data[peak_x - box_radius, peak_y - box_radius:peak_y + box_radius + 1]
    	last_row = data[peak_x + box_radius, peak_y - box_radius:peak_y + box_radius + 1]
    	
    	if (first_column == 0).all() and (first_row == 0).all() and (last_column == 0).all() and (last_row == 0).all():
    		return data[peak_x - box_radius : peak_x + box_radius + 1, peak_y - box_radius : peak_y + box_radius + 1], errors[peak_x - box_radius : peak_x + box_radius + 1, peak_y - box_radius : peak_y + box_radius + 1]
    	else:
    		box_radius += 1

def find_flux_fit(star_box, error_box, peak_data, exposure_time):
    # peak_data is a list of the shape
    # peak_data = [peak_x, peak_y, peak_height]
    # If the algorithm is used for double gaussians, then
    # peak_data = [peak_x1, peak_y1, peak_height1,
    #              peak_x2, peak_y2, peak_height2]
    
    x_length = star_box.shape[0]
    y_length = star_box.shape[1]
    
    x = linspace(-1, 1, x_length)
    y = linspace(-1, 1, y_length)
    X, Y = meshgrid(x, y)
    Z = star_box
    errors = error_box
    
    X = X.reshape(-1)
    Y = Y.reshape(-1)    
    
    Z = Z.reshape(-1)
    
    errors = errors.reshape(-1)
    
    cutout = where(Z > saturation_threshold)[0]
    
    X = delete(X, cutout)
    Y = delete(Y, cutout)
    Z = delete(Z, cutout)
    errors = delete(errors, cutout)

    # Now comes the tricky part: Since the gaussian approximation is just that, an approximation, it will eventually yield a wrong fit in some sense.
    # Therefore, we can only fit with data that is above a certain treshold (this makes sense, as can be seen from the shape of the airy pattern, which
    # is supposed to be more precise, but difficult to implement for some reason). This manifests for us another, lower treshold. However, setting
    # this treshold too high might make it difficult to find enough points to fit against. Therefore, we try to set the treshold such that there are at
    # least 4 points in x- and y-direction at the position of the estimated peak, which lie above the threshold.

    lower_threshold = initial_lower_threshold + 1000.0

    x_points = 0
    y_points = 0
    while x_points < 4 or y_points < 4:
        lower_threshold  -= 1000.0
        x_threshold, y_threshold = where((star_box > lower_threshold) & (star_box < saturation_threshold))
        x_points = len(where(x_threshold == peak_data[0])[0])
        y_points = len(where(y_threshold == peak_data[1])[0])
    

    lower_cutout = where(Z < lower_threshold)[0]

    X = delete(X, lower_cutout)
    Y = delete(Y, lower_cutout)
    Z = delete(Z, lower_cutout)
    errors = delete(errors, lower_cutout)
    
    data = row_stack((X, Y))
    
    initial_guess = [peak_data[2] * (2 * pi * 0.01), x[peak_data[0]], y[peak_data[1]], 0.01, 0.0, 0.01]
    parameter_ideal, parameter_error, p_value, SSR = general_fit(data, Z, single_gaussian, initial_guess, y_err = errors)
        
    flux = parameter_ideal[0]/exposure_time
    flux_error = parameter_error[0]/exposure_time
    return flux, flux_error, p_value, SSR, parameter_ideal, parameter_error, lower_threshold

def find_flux(star_file, error_file, create_plot = False, location = None):
    data = fits.getdata(star_file, ext = 0)
    errors = fits.getdata(error_file, ext = 0)
    
    header = fits.getheader(star_file, ext = 0)
    exposure_time = header['EXPOSURE']
    
    peaks = find_peaks(data, threshold_star)
    if len(peaks) == 0:
            print ''
            print '!!! Emergency threshold had to be used. !!!'
            peaks = find_peaks(data, emergency_threshold)

    peak_index = where( array(peaks['peak_value']) == max(array(peaks['peak_value'])))[0][0]
    peak_y = peaks['x_peak'][peak_index]
    peak_x = peaks['y_peak'][peak_index]

    box, error_box = find_peak_base_box(data, errors, peak_x, peak_y)
    
    peak_height = data[peaks['y_peak'][peak_index], peaks['x_peak'][peak_index]]

    peaks_box = find_peaks(box, threshold_star)
    if len(peaks_box) == 0:
            print ''
            print '!!! Emergency threshold had to be used. !!!'
            peaks_box = find_peaks(box, emergency_threshold)

    peak_index_box = where( array(peaks_box['peak_value']) == max(array(peaks_box['peak_value'])))[0][0]
    peak_y_box = peaks_box['x_peak'][peak_index_box]
    peak_x_box = peaks_box['y_peak'][peak_index_box]

    
    
    peak_data = [peak_x_box, peak_y_box, peak_height]
    
    flux, flux_error, p_value, SSR, parameter_ideal, parameter_error, lower_threshold_used = find_flux_fit(box, error_box, peak_data, exposure_time)
    
    if create_plot:

        x = linspace(-1, 1, box.shape[0])
        y = linspace(-1, 1, box.shape[1])
        X, Y = meshgrid(x, y)

        box_fitted = single_gaussian(parameter_ideal, [X, Y])

        i_want_3d = False

        figure(2)

        if i_want_3d:
            from mpl_toolkits import mplot3d
            fig = figure(2)
            ax = axes(projection = '3d')
        
            ax.plot_surface(X, Y, box, rstride = 1, cstride = 1, cmap = 'viridis', edgecolor = 'none')
               
            ax.scatter(X, Y, box_fitted, c = 'b')
        else:
            errorbar(X[peak_x_box, :], (box[peak_x_box, :]+1), label = 'x-axis', yerr = error_box[peak_x_box, :])
            errorbar(Y[:, peak_y_box], (box[:, peak_y_box]+1), label = 'y-axis', yerr = error_box[:, peak_y_box])
            plot(X[peak_x_box, :], (box_fitted[peak_x_box, :]+1), '--', label = 'Gaussian fit, x-axis')
            plot(Y[:, peak_y_box], (box_fitted[:, peak_y_box]+1), '--', label = 'Gaussian fit, y-axis')
            #plot(X[peak_x_box, :], (saturation_threshold) + 0.0 * X[peak_x_box, :])
            #plot(X[peak_x_box, :], (initial_lower_threshold) + 0.0 * X[peak_x_box, :])
            #plot(X[peak_x_box, :], (lower_threshold_used) + 0.0 * X[peak_x_box, :])
            legend()
            xlabel("Relative position with respect to the peak", fontsize = 20)
            ylabel("Photon count", fontsize = 20)
            title("Distribution of photon counts vertically and horizontally of the peak,\ntogether with the Gaussian fit of the curve.", fontsize = 15)
            ylim(bottom = 0.0)
        savefig(location + ".pdf")
        clf()
        close(2)
        
    return flux, flux_error
    

