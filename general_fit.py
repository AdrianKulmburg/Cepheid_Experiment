from numpy import *
from scipy.odr import *
from scipy.stats import t, chi2, kstwobign

maxit = 1000

def cdf(dat1, dat2):
    cdfnew = []
    n = dat2.shape[0]
    for yy in dat1:
       fr = len(where(dat2 <= yy)[0])/float(n)
       cdfnew.append(fr)
    return asarray(cdfnew)

def general_fit(x, y, model_to_fit, initial_guess, x_err = None, y_err = None, return_covariance = False):
    # Function that fits a given model-function to a set of data.
    # x, y have to be numpy arrays of data to fit, x_err, y_err either equal to None in case there are no
    # errors or numpy arrays giving the erros at each point
    # model_to_fit is a function of the following form:
    #
    # def model_to_fit(parameters, data):
    #     return parameters[0] * data + parameters[1]
    #
    # where parameters is a list (and NOT a numpy array) for the parameters of the function.
    # In this case, parameters = [a, b] where a, b are the parameters for a linear function.
    # initial_guess is an initial guess for the parameters (needed as general_fit is more or less
    # just finding out a zero in the error of the fit, and thus functions similarly to fsolve for
    # example).
    # If return_covariance = True, the covariance matrix of the parameters is given as an output 
    odr_model = Model(model_to_fit)

    alpha = 0.05

    if (y_err is None) and (x_err is None):
        odr_data = RealData(x, y)
    elif x_err is None:
        odr_data = RealData(x, y, sy = y_err)
    elif y_err is None:
        odr_data = RealData(x, y, sx = x_err)
    else:
        odr_data = RealData(x, y, sx = x_err, sy = y_err)

    if len(x.shape) == 1:
        degrees_of_freedom = y.shape[0] - len(initial_guess)
    elif len(x.shape) == 2:
        degrees_of_freedom = y.shape[0]**2 - len(initial_guess)

    myodr = ODR(odr_data, odr_model, beta0 = initial_guess, maxit = maxit)

    odr_output = myodr.run()
    odr_parameter_ideal = odr_output.beta
    odr_parameter_error = odr_output.sd_beta
    
    if y_err is None: # When there is no error on the y-axis, we actually need to use a different type of test, namely the Kolmogorov-Smirnov
        # test.

        # Reference:
        # https://www.astro.rug.nl/software/kapteyn/kmpfittutorial.html

        # First we need to compute the ECDF and CDF
        N = y.shape[0]
        z = sort(y)
        cdf_data_high = arange(1.0, N + 1.0)/N
        cdf_data_lo = arange(0.0, 1.0*N)/N
        X = linspace(min(x), max(x), 300)
        fitted_values = model_to_fit(odr_parameter_ideal, X)
        fitted_values.sort()
        cdf_fitted = cdf(z, fitted_values)

        # Find biggest distance between cumulative functions
        DD = cdf_data_high - cdf_fitted
        Dip = DD.argmax()
        Dplus = DD[Dip]
        DD = cdf_fitted - cdf_data_lo
        Dim = DD.argmax()
        Dmin = DD[Dim]
        if Dplus >= Dmin:
            Di = Dip
        else:
            Di = Dim
        Dmax = max(Dplus, Dmin)

        rv = kstwobign()

        odr_p_value = 1.0 - rv.cdf(Dmax)

        odr_SSR = sum((y - model_to_fit(odr_parameter_ideal, x))**2)
        # Note that the SSR is here pretty meaningless, and only expresses what the error of the fit at each point is.

    else:
        # This is the standard case, where we use a chi squared test
        chi_squarred = sum((y - model_to_fit(odr_parameter_ideal, x))**2/(y_err**2))
        odr_p_value = 1.0 - chi2.cdf(chi_squarred, degrees_of_freedom)
        odr_SSR = sum((y-model_to_fit(odr_parameter_ideal, x))**2)

    if return_covariance:
        odr_covariance = odr_output.cov_beta
        return odr_parameter_ideal, odr_parameter_error, odr_p_value, odr_SSR, odr_covariance
    else:
        return odr_parameter_ideal, odr_parameter_error, odr_p_value, odr_SSR

    # For the output: odr_parameter_ideal gives the ideal parameters for the fit, odr_parameter_error their standard error
    # odr_p_value is the p-value of a chi squared test of the fit
    # odr_SSR is the Sum of Squared Errors of the fit
    # odr_covariance is the covariance matrix of the different parameters in odr_parameter_ideal
    # odr_parameter_ideal, odr_parameter_error and odr_covariance are given as lists, i.e. not numpy arrays
