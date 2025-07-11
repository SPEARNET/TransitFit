'''
This is a series of utility functions for TransitFit, including input validation
'''
import pandas as pd
import random
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
from .lightcurve import LightCurve
from collections.abc import Iterable
import os

def validate_variable_key(key):
    '''
    Checks that a key is valid for use with PriorInfo, and corrects when it
    is obvious what is meant. Raises KeyError if unable to correct.

    '''
    if key.lower() in ['p','period']:
        return 'P'
    if key.lower() in ['rp','r_p','radius','planet_radius', 'planet radius']:
        return 'rp'
    if key.lower() in ['a', 'semimajor axis', 'semimajor_axis']:
        return 'a'
    if key.lower() in ['inc','inclination']:
        return 'inc'
    if key.lower() in ['t0', 't_0']:
        return 't0'
    if key.lower() in ['ecc','e','eccentricity']:
        return 'ecc'
    if key.lower() in ['w', 'periastron', 'longitude_of_periastron', 'longitude of periastron']:
        return 'w'
    if key.lower() in ['q0']:
        return 'q0'
    if key.lower() in ['q1']:
        return 'q1'
    if key.lower() in ['q2']:
        return 'q2'
    if key.lower() in ['q3']:
        return 'q3'
    if key.lower() in ['p_prime']:
        return 'p_prime'
    if key.lower() in ['p_dprime']:
        return 'p_dprime'

    raise KeyError('Unable to recognise variable name {}'.format(key))


def calculate_logg(host_mass, host_radius):
    '''
    Calculates log10(g) for a host in units usable by TransitFit.
    g is in cm s^-2

    Parameters
    ----------
    host_mass : tuple
        The host mass in solar masses and the uncertainty
    host_radius : tuple
        The host radius in solar radii and the uncertainty
    '''
    m_sun = 1.989e30
    r_sun = 6.957e8
    G = 6.674e-11

    m = host_mass[0] * m_sun
    err_m = host_mass[1] * m_sun
    r = host_radius[0] * r_sun
    err_r = host_radius[1] * r_sun

    g = m * G * r ** -2 * 100 # factor 100 converts to cm s^-2
    err_g = G/(r ** 2) * np.sqrt(err_m ** 2 + (4 * m**2 * err_r ** 2)/r**2) * 100

    logg = np.log10(g)
    err_logg = err_g / (g * np.log(10))

    return logg, err_logg


def get_normalised_weights(results):
    '''
    Obtains normalised weights for the Dynesty results

    Parameters
    ----------
    results : dynesty.results.Results
        The Dynesty results object

    Returns
    -------
    weights : np.array, shape (n_iterations,)
        The normalised weights for each sample set
    '''
    return np.exp(results.logwt - results.logwt.max())/np.sum(np.exp(results.logwt - results.logwt.max()))


def get_covariance_matrix(results):
    '''
    Gets a covariance matrix from Dynesty results

    Parameters
    ----------
    results : dynesty.results.Results
        The Dynesty results object, but must also have weights as an entry

    Returns
    -------
    cov : np.array, shape (ndims, ndims)
        The covariance matrix for the results object.
    '''
    # Calculate a covariance matrix using numpy
    cov = np.cov(results.samples, rowvar=False, aweights=results.weights)

    return cov


def weighted_avg_and_std(values, weights, axis=-1, single_val=False):
    '''
    Calculates the weighted average and error on some data.

    axis defaults to the last one (-1)
    '''
    if not isinstance(values, Iterable):
        values = [values]
    if not isinstance(weights, Iterable):
        weights = [weights]

    values = np.array(values)
    weights = np.array(weights)

    if not isinstance(values[0], Iterable):
        single_val = True

    if single_val:
        # Flatten the values and weights
        try:
            flat_vals = np.array([i for epoch_vals in values for i in epoch_vals])
            flat_weights = np.array([i for epoch_vals in weights for i in epoch_vals])
        except:
            flat_vals = values
            flat_weights = weights
        try:
            average = np.average(flat_vals, weights=1/(flat_weights**2))
            uncertainty = 1 / np.sqrt(np.sum(1/(flat_weights**2)))
        except TypeError:
            average = np.average(flat_vals)
            uncertainty = 0

        return average, uncertainty

    # Make blank arrays to loop over the entries in values and weights
    average = []
    uncertainty = []

    for i in range(len(values)):
        average.append(np.average(np.array(values[i]), weights=1/(np.array(weights[i])**2)))
        uncertainty.append(1 / np.sqrt(np.sum(1/(np.array(weights[i])**2))))

    return np.array(average), np.array(uncertainty)


def AU_to_host_radii(a, R, a_err=0, R_err=0, calc_err=False):
    '''
    Converts a number in AU to a value in host radii when given the host
    radius R in Solar radii. Inverse of host_radii_to_AU.
    '''
    AU = 1.495978707e11
    R_sun = 6.957e8

    if not calc_err:
        return (a * AU) / (R * R_sun)

    err = np.sqrt( ((R_err * R_sun) * (a*AU)/((R*R_sun)**2)) ** 2 + ((a_err * AU)/(R*R_sun))**2  )

    return (a * AU) / (R * R_sun), err/(R * R_sun)


def host_radii_to_AU(a, R, a_err=0, R_err=0, calc_err=False):
    '''
    converts a separation in host radii into a separation in AU when given
    the host radius R in Solar radii. Inverse of AU_to_host_radii.
    '''
    AU = 1.495978707e11
    R_sun = 6.957e8

    if not calc_err:
        return (a * R * R_sun)/AU

    try:
        err = np.sqrt ((R * R_sun * a_err )**2 + (a * R_err * R_sun)**2)/AU
    except TypeError:
        err='-'

    return (a * R * R_sun)/AU, err


def estimate_t14(Rp, Rs, a, P):
    '''
    Estimates t14 in minutes, if P is in days
    '''
    AU = 1.495978707e11
    R_sun = 6.957e8
    R_jup = 71492000

    return (Rp * R_jup + Rs * R_sun)/(np.pi * a * AU) * P * 24 * 60

def split_lightcurve_file(path, t0, P,t14=20, cutoff=0.25, window=5,
                          new_base_fname='split_curve'):
    '''
    Split a light curve file into multiple single epoch files

    Splits a multi-epoch lightcurve data file into multiple single-epoch files
    and saves these. This is useful for dealing with multi-epoch observations
    which contain TTVs, or have long term periodic trends, since single-epoch
    observation trends can be approximated with polynomial fits. New files
    are created of the form new_base_name_X

    Parameters
    ----------
    path : str
        The path to the light curve data to be split
    t0 : float
        The centre of a transit
    P : float
        The estimated period of the planet in days
    t14 : float, optional
        The approximate transit duration in minutes. Default is 20
    cutoff : float, optional
        If there are no data within t14 * cutoff of t0, a period will be
        discarded. Default is 0.25
    window : float, optional
        Data outside of the range [t0 ± (0.5 * t14) * window] will be
        discarded. Default is 5.
    new_base_fname : str, optional
        The base name for the new files, which will have numbers appended
        depending on the epoch. This can be used to specify a relative path
        for saving. Default is `'split_curve'`.
    return_names : bool


    Returns
    -------
    paths : array_like, shape (n_curves,)
        The paths to each of the new data files
    '''
    from .io import read_data_file

    # Load the full data file as a LightCurve
    full_lightcurve = LightCurve(*read_data_file(path))

    # Split the full curve into individual epochs
    single_epoch_curves = full_lightcurve.split(t0, P, t14, cutoff, window)
    dirname = os.path.dirname(path)

    paths = []

    # Now save all the new curves
    for i, curve in enumerate(single_epoch_curves):
        fname = new_base_fname + '_{}'.format(i)
        curve.save(os.path.join(dirname, fname))
        paths.append(os.path.join(dirname, fname))

    return paths


def check_batches(allow_ttv, data_files):
    """Taken from firefly/_utils.py, with minor modifications.
    This was needed to prevent uneven batchsizes.
    Now it is not required in most of the cases, but has been kept as a failsafe.

    Args:
        allow_ttv (bool): To allow for TTV or not.
        data_files (str): The path to the data input .csv file, which contains the paths to the
        light curve data.

    Returns:
        str:  The path to modified data input .csv file, which contains the paths to the
        light curve data. It removed some lightcurves randomly to take care of batching issue in uneven batches.
    """

    print('Checking if batches can be made of equal size.')
    from natsort import natsorted
    df = pd.read_csv(data_files)
    light_curve_paths = df['Path'].to_list()

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Enforce batches be the same size
    curves = len(light_curve_paths)
    if int(curves) > 5:
        print('Still checking...')
        if allow_ttv == True:
            equal_batches = [5 + 3 * n for n in range(0, 500)]
        else:
            equal_batches = [6 + 4 * n for n in range(0, 500)]

        if int(curves) not in equal_batches:
            import bisect

            index = bisect.bisect(equal_batches, int(curves))
            new_curves = equal_batches[index - 1]
            print(
                f"\nEnforcing lightcurves to be in equal batch sizes. Discarded {curves-new_curves}."
            )
            curves = new_curves
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Sort the files into ascending order and take random sample
    if len(light_curve_paths) == int(curves):
        return data_files

    light_curve_paths = random.sample(light_curve_paths, k=int(curves))
    light_curve_paths = natsorted(light_curve_paths)
    print(
        f"\nA total of {len(light_curve_paths)}",
        "lightcurves" if len(light_curve_paths) > 1 else "lightcurve",
        "will be used.",
    )

    path_temp = os.path.dirname(os.path.abspath(data_files))
    data_path = f"{path_temp}/data_paths_modified.csv"
    df = df.loc[df['Path'].isin(light_curve_paths)]
    df.to_csv(data_path, index=False, header=True)
    print(f"New data path is {data_path}")

    return data_path
