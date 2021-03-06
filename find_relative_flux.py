from numpy import *
from general_fit import *
from astropy.io import fits
from photutils import find_peaks
from matplotlib.pyplot import *

threshold_star = 7000.0
emergency_threshold = 3000.0
threshold_secondary_star = 1000.0

initial_lower_threshold = 20000
saturation_threshold = 38000

maximum_padding_distance = 50.0

threshold_removal = 20.0

i_want_3d = True
i_want_2d = True

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

def find_peak_base_box(data, errors, peak_x, peak_y, side = False):
    """Finds the base of a given peak, where we assume
    that the threshold is zero.
    Note that for this, we assume that the peak is in
    middle of the picture, i.e. there should be no way
    to get an out of bounds error. Furthermore, we also
    assume that the source is separated enough from other
    sources for them not to impact our result."""

    box_radius = 1
    padding = 10
    while 1:
        first_column = data[peak_x - box_radius:peak_x + box_radius + 1, peak_y - box_radius]
        first_row = data[peak_x - box_radius:peak_x + box_radius + 1, peak_y + box_radius]
        last_column = data[peak_x - box_radius, peak_y - box_radius:peak_y + box_radius + 1]
        last_row = data[peak_x + box_radius, peak_y - box_radius:peak_y + box_radius + 1]
        
        if (first_column == 0).all() and (first_row == 0).all() and (last_column == 0).all() and (last_row == 0).all():
            mask = zeros_like(data)
            if side == False:
                mask[peak_x - box_radius:peak_x + box_radius + 1, peak_y - box_radius:peak_y + box_radius + 1] = 1
            else:
                mask[peak_x - box_radius - padding:peak_x + box_radius + 1 + padding, peak_y - box_radius - padding:peak_y + box_radius + 1 + padding] = 1
            return data[peak_x - box_radius : peak_x + box_radius + 1, peak_y - box_radius : peak_y + box_radius + 1], errors[peak_x - box_radius : peak_x + box_radius + 1, peak_y - box_radius : peak_y + box_radius + 1], mask
        else:
            box_radius += 1

def find_flux_fit(star_box, error_box, peak_data, exposure_time):
    # peak_data is a list of the shape
    # peak_data = [peak_x, peak_y, peak_height].
    
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

def find_relative_flux(star_file, error_file, defining_distance, create_plot = False, location = None):
    data = fits.getdata(star_file, ext = 0)
    errors = fits.getdata(error_file, ext = 0)
    
    header = fits.getheader(star_file, ext = 0)
    exposure_time = header['EXPOSURE'] # With this new method, the exposure time is actually meaningless; therefore, we set to 1.0 so that it doesn't impact
    # anything, while keeping the structure of the program intact, in case we need to revert back those changes for some reason.
    exposure_time = 1.0
    
    peaks = find_peaks(data, threshold_star)
    if len(peaks) == 0:
            print ''
            print '!!! Emergency threshold had to be used. !!!'
            peaks = find_peaks(data, emergency_threshold)

    peak_index = where( array(peaks['peak_value']) == max(array(peaks['peak_value'])))[0][0]
    peak_y = peaks['x_peak'][peak_index]
    peak_x = peaks['y_peak'][peak_index]

    box, error_box, box_mask = find_peak_base_box(data, errors, peak_x, peak_y)
    
    peak_height = data[peaks['y_peak'][peak_index], peaks['x_peak'][peak_index]]

    peaks_box = find_peaks(box, threshold_star)
    if len(peaks_box) == 0:
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

        figure(2)

        if i_want_3d:
            from mpl_toolkits import mplot3d
            fig = figure(2)
            ax = axes(projection = '3d')
        
            cs = ax.plot_surface(X, Y, box_fitted, rstride = 1, cstride = 1, cmap = 'viridis', edgecolor = 'none', alpha = 0.5)
               
            ax.scatter(X, Y, box, c = 'b', marker = '.', label = 'Measured\nData')
            
            ax.set_xlabel('Vertical Direction', fontsize = 15, labelpad = 8.0)
            ax.set_ylabel('Horizontal Direction', fontsize = 15, labelpad = 8.0)
            ax.set_zlabel('Photon Count', fontsize = 15, labelpad = 10.0)
            ax.set_zlim(bottom = 0.0)
            ax.set_xlim(-0.95, 0.95)
            ax.set_ylim(-0.95, 0.95)
            ax.set_title('Distribution of photon counts\ntogether with a Gaussian fit.', fontsize = 15)
            
            #viridis_proxy = Rectangle((0, 0), 1, 1, fc="viridis")
            #blue_proxy = Rectangle((0, 0), 1, 1, fc="b")
            
            #ax.legend([viridis_proxy, blue_proxy], ['Measured Data', 'Gaussian Fit'])
            ax.legend(loc = 'center left')
            
            cbar = fig.colorbar(cs)
            cbar.ax.set_ylabel('Gaussian Fit', fontsize = 15)
            
            savefig(location + "_main_3d.pdf")
            clf()
            close(2)

        if i_want_2d:
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
            savefig(location + "_main.pdf")
            clf()
            close(2)

    data = data - data * box_mask

    found = False
    starting_threshold = threshold_secondary_star

    need_to_reduce = False
    while not found:
        if starting_threshold <= 0.0:
            print ""
            print "Critical error encountered. Abandonning the computation of this star."
            raise ValueError
        side_peaks = find_peaks(data, starting_threshold)
        if len(side_peaks) == 0:
            if need_to_reduce:
                print "-> " + str(starting_threshold - threshold_removal),
            else:
                print ""
                print "Scaling back threshold to " + str(starting_threshold - threshold_removal),
                need_to_reduce = True
            starting_threshold -= threshold_removal
            continue

        defining_index = 0
        padding_distance = 0.0
        while not found and padding_distance < maximum_padding_distance:
            distances = sqrt((array(side_peaks['y_peak']) - peak_x) ** 2 + (array(side_peaks['x_peak']) - peak_y) ** 2)
            possibilities = where((distances <= defining_distance + padding_distance) & (distances >= defining_distance - padding_distance))[0]
            if len(possibilities) == 0:
                padding_distance += 0.01
            elif len(possibilities) > 1:
                print ""
                print "Oh oh, several possibilities where found as reference stars..."
                print ""
                found = True
                if need_to_reduce:
                    print "and resuming...",
                defining_index = where(possibilities == max(possibilities))[0][0]
            else:
                found = True
                if need_to_reduce:
                    print "and resuming...",
                defining_index = possibilities[0]
            if padding_distance >= maximum_padding_distance:
                
                if need_to_reduce:
                    print "-> " + str(starting_threshold - threshold_removal),
                else:
                    print ""
                    print "Scaling back threshold to " + str(starting_threshold - threshold_removal),
                    need_to_reduce = True
                starting_threshold -= threshold_removal
        

    side_peak_index = defining_index
    control = distances[side_peak_index]
    side_peak_y = side_peaks['x_peak'][side_peak_index]
    side_peak_x = side_peaks['y_peak'][side_peak_index]

    side_box, side_error_box, side_box_mask = find_peak_base_box(data, errors, side_peak_x, side_peak_y)
    
    side_peak_height = data[side_peaks['y_peak'][side_peak_index], side_peaks['x_peak'][side_peak_index]]

    side_peaks_box = find_peaks(side_box, threshold_secondary_star)
    if len(side_peaks_box) == 0:
            side_peaks_box = find_peaks(side_box, 0.0)

    side_peak_index_box = where( array(side_peaks_box['peak_value']) == max(array(side_peaks_box['peak_value'])))[0][0]
    side_peak_y_box = side_peaks_box['x_peak'][side_peak_index_box]
    side_peak_x_box = side_peaks_box['y_peak'][side_peak_index_box] 
    
    side_peak_data = [side_peak_x_box, side_peak_y_box, side_peak_height]
    
    side_flux, side_flux_error, side_p_value, side_SSR, side_parameter_ideal, side_parameter_error, side_lower_threshold_used = find_flux_fit(side_box, side_error_box, side_peak_data, exposure_time)

    if create_plot:

        x = linspace(-1, 1, side_box.shape[0])
        y = linspace(-1, 1, side_box.shape[1])
        X, Y = meshgrid(x, y)

        side_box_fitted = single_gaussian(side_parameter_ideal, [X, Y])

        figure(3)

        if i_want_3d:
            from mpl_toolkits import mplot3d
            fig = figure(3)
            ax = axes(projection = '3d')
        
            cs = ax.plot_surface(X, Y, side_box_fitted, rstride = 1, cstride = 1, cmap = 'viridis', edgecolor = 'none', alpha = 0.5)
               
            ax.scatter(X, Y, side_box, c = 'b', marker = '.', label = 'Measured\nData')
            
            ax.set_xlabel('Vertical Direction', fontsize = 15, labelpad = 8.0)
            ax.set_ylabel('Horizontal Direction', fontsize = 15, labelpad = 8.0)
            ax.set_zlabel('Photon Count', fontsize = 15, labelpad = 10.0)
            ax.set_zlim(bottom = 0.0)
            ax.set_xlim(-0.95, 0.95)
            ax.set_ylim(-0.95, 0.95)
            ax.set_title('Distribution of photon counts\ntogether with a Gaussian fit.', fontsize = 15)
            
            #viridis_proxy = Rectangle((0, 0), 1, 1, fc="viridis")
            #blue_proxy = Rectangle((0, 0), 1, 1, fc="b")
            
            #ax.legend([viridis_proxy, blue_proxy], ['Measured Data', 'Gaussian Fit'])
            ax.legend(loc = 'center left')
            
            cbar = fig.colorbar(cs)
            cbar.ax.set_ylabel('Gaussian Fit', fontsize = 15)
            
            savefig(location + "_secondary_3d.pdf")
            clf()
            close(3)
        if i_want_2d:
            errorbar(X[side_peak_x_box, :], (side_box[side_peak_x_box, :]+1), label = 'x-axis', yerr = side_error_box[side_peak_x_box, :])
            errorbar(Y[:, side_peak_y_box], (side_box[:, side_peak_y_box]+1), label = 'y-axis', yerr = side_error_box[:, side_peak_y_box])
            plot(X[side_peak_x_box, :], (side_box_fitted[side_peak_x_box, :]+1), '--', label = 'Gaussian fit, x-axis')
            plot(Y[:, side_peak_y_box], (side_box_fitted[:, side_peak_y_box]+1), '--', label = 'Gaussian fit, y-axis')
            #plot(X[side_peak_x_box, :], (saturation_threshold) + 0.0 * X[side_peak_x_box, :])
            #plot(X[side_peak_x_box, :], (initial_lower_threshold) + 0.0 * X[side_peak_x_box, :])
            #plot(X[side_peak_x_box, :], (lower_threshold_used) + 0.0 * X[side_peak_x_box, :])
            legend()
            xlabel("Relative position with respect to the peak", fontsize = 20)
            ylabel("Photon count", fontsize = 20)
            title("Distribution of photon counts vertically and horizontally of the peak,\ntogether with the Gaussian fit of the curve.", fontsize = 15)
            ylim(bottom = 0.0)
            savefig(location + "_secondary.pdf")
            clf()
            close(3)

    distance = sqrt((peak_x - side_peak_x)**2 + (peak_y - side_peak_y)**2)

    if create_plot:
        original_data = fits.getdata(star_file, ext = 0)
        original_data = original_data - original_data * side_box_mask + side_box_mask * 20000.0
        fits.writeto(location + '_secondary_detection_d_' + str(distance) + '.fit', original_data, header, overwrite = True)
        
        
    return flux, side_flux, flux_error, side_flux_error, control
    

