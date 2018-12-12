from numpy import *
from matplotlib.pyplot import *
import pickle
from astropy.stats import LombScargle
from scipy.signal import find_peaks
import os
import math
from general_fit import *

number_of_days = 3

def relative_phase_basis(period, data_point):
    n = 0.0
    while data_point - n*period > 0.0:
        n += 1
    if data_point - n*period == 0.0:
        return 0.0
    else:
        return data_point - (n-1)*period
def relative_phase(period, data):
    return vectorize(lambda x : relative_phase_basis(period, x))(data)

def plot_given_events(data, events):
    for event in events:
        flux = data['fluxes'][event]
        flux_errors = data['fluxes_errors'][event]
        phase = data['phases'][event]
        errorbar(phase, flux, yerr = flux_errors, linestyle = '', marker = 'x')
        
def print_all_events(data):
    for i in data:
        for j in data[i]:
            print(j)
        print('--')

def print_all_stars_separate_days(data, star_name, displayed_star_name, secondary_magnitude, secondary_magnitude_error, ylimit, period_limits, threshold_selection, number_of_terms = 1, precision = 100000):
    T_min = period_limits[0]
    T_max = period_limits[1]
    X = array([])
    Y = array([])
    Yerr = array([])
    fig_subplots, axs = subplots(1, number_of_days, sharey=True)
    fig_subplots.subplots_adjust(wspace=0.0)
    for i in range(number_of_days):
        #figure(0)
        #subplot(1, number_of_days, i+1)
        for event in data['fluxes']:
            if ('Day' + str(i+1) in event) and (star_name in event):
                x = data['phases'][event]
                y = data['fluxes'][event]
                yerr = data['fluxes_errors'][event]
                magnitudes = secondary_magnitude - 2.5 * log10(y)
                
                magnitudes_errors = sqrt((2.5 * log10(e) * yerr / y)**2 + secondary_magnitude_error**2)
                axs[i].errorbar(x, magnitudes, yerr = magnitudes_errors, linestyle = '', marker = 'x', color = 'b')

                X = concatenate([X, x])
                Y = concatenate([Y, magnitudes])
                Yerr = concatenate([Yerr, magnitudes_errors])
        #axs[i].set_ylim(ylimit[0], ylimit[1])
##        figure(2)
##        for event in data['fluxes']:
##            if ('Day' + str(i+1) in event) and (star_name in event):
##                x = data['phases'][event]
##                y = data['fluxes'][event]
##                yerr = data['fluxes_errors'][event]
##                magnitudes = secondary_magnitude - 2.5 * log10(y)
##                
##                magnitudes_errors = 2.5 * log10(e) * yerr / y
##                errorbar(x, magnitudes, yerr = magnitudes_errors, linestyle = '', marker = 'x', color = 'r')
##        ylim(ylimit)
    order = argsort(X)
    X = X[order]
    Y = Y[order]
    Yerr = Yerr[order]
    #frequency, power = LombScargle(X, Y, Yerr, nterms = number_of_terms).autopower(nyquist_factor = 1.0)#(minimum_frequency = 1.0 / T_max, maximum_frequency = 1.0 / T_min)
    frequency = linspace(1.0/T_max, 1.0/T_min, precision)
    power = LombScargle(X, Y, Yerr, nterms = number_of_terms).power(frequency)
    fig_periodogram, ax_periodogram = subplots(1, 1)
    ax_periodogram.plot(frequency, power)
    ax_periodogram.set_title('Lomb Scargle Periodogram for ' + displayed_star_name, fontsize = 18)
    ax_periodogram.set_xlabel(r'Frequencies in $[days]^{-1}$', fontsize = 20)
    ax_periodogram.set_ylabel('Frequency power', fontsize = 20)
    ax_periodogram.tick_params(top = True, bottom = True, left = True, right = True)
    ax_periodogram.set_ylim(bottom = 0.0)

    peaks = find_peaks(power, threshold_selection)

    try:
        selected_peak = peaks[0][0]
    except IndexError:
        selected_peak = argmax(power)
    
    best_frequency = frequency[selected_peak]
    print('Selected peak:', best_frequency, 'with power', power[selected_peak])
    
    noise_level = median(power)
    peak_height = power[selected_peak] - noise_level
    half_maximum = peak_height / 2.0 + noise_level
    left = where(power <= half_maximum)[0][:selected_peak][-1]
    right = where(power <= half_maximum)[0][selected_peak:][0]
    hwhm = (frequency[right] - frequency[left])/2.0

    frequency_uncertainty = hwhm / sqrt(power[selected_peak])

    period_uncertainty = frequency_uncertainty / best_frequency ** 2
    
    p_false_alarm = LombScargle(X, Y, Yerr, nterms = number_of_terms).false_alarm_probability(best_frequency, minimum_frequency = 1.0/T_max, maximum_frequency = 1.0/T_min)
    p_level = LombScargle(X, Y, Yerr, nterms = number_of_terms).false_alarm_level(0.001)
    print("The estimated period is", 1.0 / best_frequency, "+-", period_uncertainty)

    print("The required power level for a probability of 0.001 is", p_level)
    print("The peak in question has a power of", power[selected_peak])
    print("The false alarm p-value is", p_false_alarm)


    space_for_best_fit = linspace(0.0001, 10000.0, 1000000)
    powers = LombScargle(X, Y, Yerr, nterms = number_of_terms).power(space_for_best_fit)
    best_fit_frequency = space_for_best_fit[argmax(powers)]
    space_for_mean = linspace(0.0, 1.0 / best_fit_frequency, 1000000)
    #values_for_mean = LombScargle(X, Y, Yerr, nterms = number_of_terms).model(space_for_mean, best_fit_frequency)
    #star_mean = values_for_mean.mean()
    relative_phases =  relative_phase(1.0 / best_fit_frequency, array(X))

    sinusoidal_model = lambda parameters, data: parameters[0] + parameters[1]*sin(2 * pi * best_fit_frequency * data + parameters[2])

    [star_mean, amplitude, delta], [s_star_mean, _, _], mean_p_value, mean_SSR = general_fit(relative_phases, array(Y), sinusoidal_model, [Y[0], 1.0, 0.0], y_err = array(Yerr))

    values_for_mean = sinusoidal_model([star_mean, amplitude, delta], space_for_mean)

    print("The mean apparent magnitude is", star_mean, "+-", s_star_mean)
    print("which was found for a frequency of", best_fit_frequency)
    print("and a p-value of", mean_p_value)

    fig_relative_diagram, ax_relative_diagram = subplots(1, 1)
    ax_relative_diagram.plot(space_for_mean, values_for_mean, 'k--', label = 'Best frequency')
    ax_relative_diagram.errorbar(relative_phases, array(Y), yerr = array(Yerr), marker = 'x', color = 'b', linestyle = '', label = 'Light curve measurements')
    ax_relative_diagram.set_title('Relative phases of the measurements\nwith respect to the best frequency', fontsize = 13)
    ax_relative_diagram.set_xlabel('Relative phase in [days]', fontsize = 15)
    ax_relative_diagram.set_ylabel(r'Apparent magnitude $m_V$', fontsize = 15)
    ax_relative_diagram.legend(loc = 'best')
    ax_relative_diagram.tick_params(top = True, bottom = True, left = True, right = True)


    
    for i in range(number_of_days):
        #figure(0)
        #subplot(1, number_of_days, i+1)
        x_intermediary = array([])
        y_intermediary = array([])
        for event in data['fluxes']:
            if ('Day' + str(i+1) in event) and (star_name in event):
                x_intermediary = concatenate([x_intermediary, data['phases'][event]])
                y_intermediary = concatenate([y_intermediary, data['fluxes'][event]])
        order = argsort(x_intermediary)
        x = x_intermediary[order]
        y = y_intermediary[order]

        xs = linspace(min(x)-1.0, max(x)+1.0, 1000)
        y_fit = LombScargle(X, Y, Yerr, nterms = number_of_terms).model(xs, best_frequency)
        axs[i].plot(xs, y_fit, 'k--')
        axs[i].set_xlim(min(x)-0.1/24.0, max(x)+0.1/24.0)
        axs[i].set_ylim(ylimit[0], ylimit[1])
        axs[i].tick_params(top = True, bottom = True, left = True, right = True)

        middle = (max(x) - min(x))/2.0
        tick1 = min(x) + middle/2.0
        tick2 = tick1 + middle
        axs[i].set_xticks([tick1, tick2])
        axs[i].grid(True)

    ax_main = fig_subplots.add_subplot(111, frameon=False)
    ax_main.tick_params(labelcolor='none', top = False, bottom = False, left = False, right = False)
    axs[0].set_ylabel(r'Apparent magnitude $m_V$', fontsize = 18)
    ax_main.set_xlabel('Phases in [days]', fontsize = 18)
    ax_main.set_title('Light curve of ' + displayed_star_name + ' in the visual band\nand best Lomb-Scargle-frequency.', fontsize = 12)
    axs[0].errorbar([], [], xerr = [], yerr = [], color = 'b', marker = 'x', linestyle = '', label = 'Light curve')
    axs[0].plot([], [], color = 'k', marker = '', linestyle = '--', label = 'Best frequency')
    axs[0].legend(loc = 'best', framealpha = 0.3)
##    figure(2)
##    x = linspace(min(X), max(X), 10000)
##    plot(x, LombScargle(X, Y, Yerr, nterms = number_of_terms).model(x, best_frequency), 'k:')
##    

data = pickle.load(open('Periodogram.pickle', 'rb'), encoding='latin1')

secondary_flux_RR_Lyr = 8.5013#8.80
secondary_flux_RR_Lyr_error = 0.0003
secondary_flux_FF_Aql = 8.4116#9.40
secondary_flux_FF_Aql_error = 0.0009
secondary_flux_V0473 = 6.8809#6.69
secondary_flux_V0473_error = 0.0004

magnitude_limits_RR_Lyr = (7.3, 8.5)
magnitude_limits_FF_Aql = (5.1, 6.0)
magnitude_limits_V0473 = (5.75, 6.65)

period_limits_RR_Lyr = (0.3, 1.0)
period_limits_FF_Aql = (4.0, 6.0)
period_limits_V0473 = (1.2, 10.0)

if not os.path.exists('Final_Results'):
    os.makedirs('Final_Results')



print('For RR  Lyr:')
print_all_stars_separate_days(data, 'RR_Lyr', 'RR Lyr', secondary_flux_RR_Lyr, secondary_flux_RR_Lyr_error, magnitude_limits_RR_Lyr, period_limits_RR_Lyr, 8.0, precision = 100000)
figure(2).savefig('Final_Results/RR_Lyr_periodogram.pdf')
figure(1).savefig('Final_Results/RR_Lyr_measurements.pdf')
figure(3).savefig('Final_Results/RR_Lyr_relative_diagram.pdf')
close(1)
close(2)
close(3)
print('')

print('For FF Aql:')
print_all_stars_separate_days(data, 'FF_Aql', 'FF Aql', secondary_flux_FF_Aql, secondary_flux_FF_Aql_error, magnitude_limits_FF_Aql, period_limits_FF_Aql, 600.0, precision = 100000)
figure(2).savefig('Final_Results/FF_Aql_periodogram.pdf')
figure(1).savefig('Final_Results/FF_Aql_measurements.pdf')
figure(3).savefig('Final_Results/FF_Aql_relative_diagram.pdf')
close(1)
close(2)
close(3)
print('')

print('For V 473 Lyr:')
print_all_stars_separate_days(data, 'V0473', 'V 473 Lyr', secondary_flux_V0473, secondary_flux_V0473_error, magnitude_limits_V0473, period_limits_V0473, 200.0, precision = 100000)
figure(2).savefig('Final_Results/V0473_periodogram.pdf')
figure(1).savefig('Final_Results/V0473_measurements.pdf')
figure(3).savefig('Final_Results/V0473_relative_diagram.pdf')
close(1)
close(2)
close(3)
print('')



            

