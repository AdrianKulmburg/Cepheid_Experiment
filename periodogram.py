from numpy import *
from datetime import *
from astropy.io import fits
from matplotlib.pyplot import *
from photutils import find_peaks
from find_flux import find_flux
from find_relative_flux import find_relative_flux
import os

approximate_total_files = 290

reference_year = 2018
reference_month = 10
reference_day = 1
reference_hour = 0
reference_minute = 0
reference_second = 0
reference_microsecond = 0

reference_date = datetime(year = reference_year, month = reference_month, day = reference_day, hour = reference_hour, minute = reference_minute, second = reference_second, microsecond = reference_microsecond)

plot_individual_fits = True
star_names = ['RR_Lyr', 'FF_Aql', 'V0473']
star_signs = ['x', '.', 's', 'v', '^']
star_color = ['r', 'b', 'g']
star_defining_distances = [1687.0, 1924.0, 1495.5]

complete_phases = {}
complete_fluxes = {}
complete_fluxes_errors = {}
complete_controls = {}

days_extension = os.listdir('./')

counter = 1

for i, star in enumerate(star_names):
    print 'Starting to process ' + star + '...'
    for day_number in days_extension:
        color_counter = 0
        if not ('Day' in day_number):
            continue
        day_directories = os.listdir(day_number)
        for possible_directory in day_directories:
            phases = []
            fluxes = []
            fluxes_errors = []
            side_fluxes = []
            side_fluxes_errors = []
            controls = []
            if not((star in possible_directory) and ('Processed' in possible_directory)):
                continue
            
            star_directory_name = day_number + '/' + possible_directory + '/'
            star_directory_errors = day_number + '/' + possible_directory[:-9] + 'Errors/'
            individual_fit_directory = "Individual_Fits/" + day_number + '/' + possible_directory[:-10] + '/'
            if not(os.path.exists(individual_fit_directory)) and plot_individual_fits:
                                os.makedirs(individual_fit_directory)

            events = os.listdir(star_directory_name)

            for event in events:
                print str(counter)+'/'+str(approximate_total_files),
                counter += 1
                print 'o-> Reading in file ' + star_directory_name + event + '...',
                header = fits.getheader(star_directory_name + event, ext = 0)
                data = fits.getdata(star_directory_name + event, ext = 0)
                print 'Done.'


                exposure_time = header['EXPOSURE']
                observed_time = header['DATE-OBS']
                event_year = int(observed_time[:4])
                event_month = int(observed_time[5:7])
                event_day = int(observed_time[8:10])
                event_hour = int(observed_time[11:13])
                event_minute = int(observed_time[14:16])
                event_second = int(observed_time[17:19])
                event_microsecond = int(int(observed_time[20:22]) / 100.0 * 1.0e-6)

                event_date = datetime(year = event_year, month = event_month, day = event_day, hour = event_hour, minute = event_minute, second = event_second, microsecond = event_microsecond)

                phase = event_date - reference_date
                phase = phase.days + (phase.seconds / 3600.0 + phase.microseconds / 3600.0 * 1.0e-6)/24.0
                

                if plot_individual_fits:
                    try:
                        
                        #flux, flux_error = find_flux(star_directory_name + event, star_directory_errors + event, create_plot = True, location = individual_fit_directory + event)
                        flux, side_flux, flux_error, side_flux_error, control = find_relative_flux(star_directory_name + event, star_directory_errors + event, star_defining_distances[i], create_plot = True, location = individual_fit_directory + event)
                    except ZeroDivisionError:#ValueError:
                        del data
                        del header
                        continue
                else:
                    try:
                        #flux, flux_error = find_flux(star_directory_name + event, star_directory_errors + event)
                        flux, side_flux, flux_error, side_flux_error, control = find_relative_flux(star_directory_name + event, star_directory_errors + event, star_defining_distances[i])
                    except ValueError:
                        del data
                        del header
                        continue
                    
                
                phases.append(phase)
                fluxes.append(flux)
                fluxes_errors.append(flux_error)
                side_fluxes.append(side_flux)
                side_fluxes_errors.append(side_flux_error)
                controls.append(control)

                del data
                del header

            if len(phases) == 0:
                continue

            phases = array(phases)
            fluxes = array(fluxes)
            fluxes_errors = array(fluxes_errors)
            side_fluxes = array(side_fluxes)
            side_fluxes_errors = array(side_fluxes_errors)

            side_flux = sum(side_fluxes / side_fluxes_errors ** 2) / sum(1.0 / side_fluxes_errors ** 2)
            side_flux_error = sqrt(1.0 / sum(1.0 / side_fluxes_errors ** 2))

            fluxes = fluxes / side_flux
            fluxes_errors = sqrt(fluxes_errors**2/side_flux**2 + fluxes**2*side_flux_error**2/side_flux**4)

            complete_phases[possible_directory[:-10]] = phases
            complete_fluxes[possible_directory[:-10]] = fluxes
            complete_fluxes_errors[possible_directory[:-10]] = fluxes_errors
            complete_controls[possible_directory[:-10]] = controls
            
            #figure(0)
            #errorbar(phases, fluxes, yerr = fluxes_errors, marker = star_signs[i], linestyle = '', color = star_color[color_counter], label = possible_directory)
            #draw()
            #pause(0.05)
            #color_counter += 1
#art = []
#lgd = legend(loc = 9, bbox_to_anchor=(0.5, -0.15))
#art.append(lgd)
#savefig('Periodogram.pdf', additional_artists=art, bbox_inches="tight")
#legend(loc = 'best')

import pickle
complete_data = {'phases':complete_phases, 'fluxes':complete_fluxes, 'fluxes_errors':complete_fluxes_errors}
pickle.dump(complete_data, file('Periodogram.pickle', 'w'))
pickle.dump(complete_controls, file('control.pickle', 'w'))

#show()



    
