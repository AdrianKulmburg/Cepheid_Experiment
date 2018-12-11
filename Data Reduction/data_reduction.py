print "Charging libraries..."
from time import time
print "-# time ready."
loading_start = time()

from numpy import *
import os
from sys import argv, stdout
from scipy.stats import sem
print "-# numpy, os, sys, sem ready."
stdout.flush()

from astropy.io import fits
print "-# astropy.io.fits ready."
stdout.flush()

from astroscrappy import detect_cosmics
print "-# astroscrappy ready."
stdout.flush()

from astropy.stats import sigma_clipped_stats
from photutils import make_source_mask
print "-# photutils ready."
print "All ready."
stdout.flush()
loading_end = time()
print "Time needed for charging libraries: ", int(loading_end-loading_start), " seconds."

iterations = 1
signal_to_noise = 4.0
sigma = 10.0

def median_uncertainty(mean_uncertainty, sample_size):
    if sample_size % 2 == 0:
        return mean_uncertainty * sqrt(pi / 2.0)
    else:
        return mean_uncertainty * sqrt(pi * sample_size / (2.0 * sample_size - 1.0))

def clean_negative(data):
    f = lambda x : 0 if x < 0 else x
    return vectorize(f)(data)

def reject_masked(data, mask):
    return data * (1 - mask)

def select_masked(data, mask):
    return data * mask

def data_reduction(new_filename, prefix):
    science_directory = prefix + '/' + new_filename + '/Science'
    flats_directory = prefix + '/' + 'Flats'
    dark_for_flats_directory = prefix + '/' + 'Dark_For_Flats'
    dark_for_science_directory = prefix + '/' + new_filename + '/Dark_For_Science'
    print science_directory
    
    reading_start = time()
    problem = (not os.path.exists(science_directory))
    problem = problem or (not os.path.exists(flats_directory))
    problem = problem or (not os.path.exists(dark_for_flats_directory))
    problem = problem or (not os.path.exists(dark_for_science_directory))

    if problem:
        print 'One of the defining directories for the data reduction'
        print 'is not present. Aborting...'
        exit()

    if not os.path.exists(prefix + '/' + new_filename + "_Processed"):
        os.makedirs(prefix + '/' + new_filename + "_Processed")

    if not os.path.exists(prefix + '/' + new_filename + "_Errors"):
        os.makedirs(prefix + '/' + new_filename + "_Errors")

    if not os.path.exists(prefix + '/' + new_filename + "_With_Background"):
        os.makedirs(prefix + '/' + new_filename + "_With_Background")

    science_files = os.listdir(science_directory)
    flats_files = os.listdir(flats_directory)
    dark_for_flats_files = os.listdir(dark_for_flats_directory)
    dark_for_science_files = os.listdir(dark_for_science_directory)

    science_data = []
    science_headers = []
    flats_data = []
    flats_headers = []
    dark_for_flats_data = []
    dark_for_flats_headers = []
    dark_for_science_data = []
    dark_for_science_headers = []

    print ""
    print ""
    print "Starting to read in science files...",
    stdout.flush()
    for item in science_files:
        science_data.append(fits.getdata(science_directory + "/" + item, ext = 0))
        science_headers.append(fits.getheader(science_directory + "/" + item, ext = 0))
    print "Done."

    print "Starting to read in flats...",
    stdout.flush()
    for item in flats_files:
        flats_data.append(fits.getdata(flats_directory + "/" + item, ext = 0))
        flats_headers.append(fits.getheader(flats_directory + "/" + item, ext = 0))
    print "Done."

    print "Starting to read in darks for flats...",
    stdout.flush()
    for item in dark_for_flats_files:
        dark_for_flats_data.append(fits.getdata(dark_for_flats_directory + "/" + item, ext = 0))
        dark_for_flats_headers.append(fits.getheader(dark_for_flats_directory + "/" + item, ext = 0))
    print "Done."

    print "Starting to read in darks for science...",
    stdout.flush()
    for item in dark_for_science_files:
        dark_for_science_data.append(fits.getdata(dark_for_science_directory + "/" + item, ext = 0))
        dark_for_science_headers.append(fits.getheader(dark_for_science_directory + "/" + item, ext = 0))
    print "Done."

    print "Starting to compute master dark, master flat..."
    stdout.flush()
    master_dark = median(array(dark_for_science_data), axis = 0)
    master_flat = median(array(flats_data) - median(array(dark_for_flats_data), axis = 0) + 0.0, axis = 0)

    master_dark_error = median_uncertainty(sem(array(dark_for_science_data), axis = 0), len(dark_for_science_data))
    master_flat_error = median_uncertainty(sem(array(flats_data) - median(array(dark_for_flats_data), axis = 0), axis = 0), len(flats_data))

    mean_master_flat = master_flat.mean()
    mean_master_flat_error = sqrt((master_flat_error ** 2).sum()) / (master_flat_error.shape[0] * master_flat_error.shape[1] + 0.0)
    print "Done."
    stdout.flush()

    print "Freeing memory...",
    stdout.flush()
    del flats_data
    del flats_headers
    del dark_for_flats_data
    del dark_for_flats_headers
    del dark_for_science_data
    del dark_for_science_headers
    print "Done."

    reading_end = time()
    print "Time needed for reading: ", int(reading_end - reading_start), " seconds."

    print "Starting to compute science frames..."
    print ""

    already_done = os.listdir(prefix + '/' + new_filename + "_Processed/")
    for i, item in enumerate(science_data):

        #if science_files[i] in already_done:
        #    continue
        item_start = time()
        
        total_errors = zeros_like(item)

        print "o Computing frame ", i+1, " out of ", len(science_data)
        print "o-> Cleaning frame...",
        stdout.flush()
        cleaned_frame = (item - master_dark + 0.0) * (mean_master_flat / master_flat)
        total_errors = sqrt((master_dark_error * mean_master_flat / master_flat) ** 2 + (mean_master_flat_error * master_dark / master_flat) ** 2 + (master_dark * mean_master_flat * master_flat_error / master_flat ** 2) ** 2 )

        print "Done."
        print "o-> Eliminating cosmics...",
        stdout.flush()
        mask, no_cosmics = detect_cosmics(cleaned_frame,
                                          gain = science_headers[i]["EGAIN"],
                                          readnoise = 9.0, satlevel = 44000,
                                          objlim = 2000)
        print "Done."
        print "o-> Eliminating background noise..."
        stdout.flush()

        corrected = no_cosmics

	# This step is just to be sure that there won't be a bug with the border
	padding = 10 # 10 pixels seems to be enough
	final = zeros_like(corrected)
	final_with_background = corrected + 0.0
	corrected = corrected[padding:-padding, padding:-padding]


        for j in xrange(iterations):
            run_start = time()
            print "oo-> Starting run ", j+1, ". Mask: ",
            stdout.flush()
            mask = make_source_mask(corrected, snr = signal_to_noise, npixels = 5, sigclip_sigma = sigma, sigclip_iters = None, dilate_size = 11)
            print "Done. Correction: ",
            stdout.flush()
            mean_mask, median_mask, std_mask = sigma_clipped_stats(corrected, sigma = sigma, iters = None, mask = mask)
            corrected = clean_negative(select_masked(corrected - median_mask, mask))
            final_with_background = final_with_background - median_mask

            sample_size = corrected.shape[0] * corrected.shape[1] - mask.sum()
            std_error = median_uncertainty(std_mask / sqrt(sample_size), sample_size)
            total_errors = sqrt(total_errors ** 2 + std_error**2)

            print "Done."
            run_end = time()
            print "oo-> Time needed: ", int(run_end - run_start), " seconds."
        print "o-> Background noise elimination done."        
        print "o-> Starting to write the file...",
        stdout.flush()

	final[padding:-padding, padding:-padding] = corrected

        fits.writeto(prefix + '/' + new_filename + "_Processed/" + science_files[i], final, science_headers[i], overwrite = True)
        fits.writeto(prefix + '/' + new_filename + "_Errors/" + science_files[i], total_errors, science_headers[i], overwrite = True)
        fits.writeto(prefix + '/' + new_filename + "_With_Background/" + science_files[i], final_with_background, science_headers[i], overwrite = True)
        print "Done."
        print ""

        print "o-> Freeing memory...",
        stdout.flush()
        del cleaned_frame
        del mask
        del no_cosmics
        del mean_mask
        del median_mask
        del std_mask
	del corrected
	del final
	del final_with_background
        del total_errors
        science_data[i] = None
        print "Done."
        print ""
        stdout.flush()
        item_end = time()
        print "Time needed for one frame: ", int(item_end - item_start), " seconds."
        print ""
        print ""
    print ""
    print "All data processed."

if len(argv) < 3:
	to_treat = []
else:
	prefix = '../' + argv[1]
	to_treat = argv[2:]

absolute_total_start = time()
for new_filename in to_treat:
	print ""
	print ""
	print "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
	print "Now starting to process " + prefix + '/' + new_filename + "."
	start_time_total = time()
	data_reduction(new_filename, prefix)
	end_time_total = time()
	print "Total time needed: ", int(end_time_total - start_time_total), " seconds."
absolute_total_end = time()
print "The total time needed for everything is: ", int(absolute_total_end - absolute_total_start), " seconds."


    

    
    
