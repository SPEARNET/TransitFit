# Inject TTV signals
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time
import os
from .ttv_fitting import get_shift_in_time_due_to_ttv,get_time_duration, taylor_series


def get_prior_value(param, priors_df):
    """    Get the prior value for a given parameter from the priors DataFrame.
    Args:
        param (str): The name of the parameter to retrieve.
        priors_df (pd.DataFrame): DataFrame containing prior information.
    Returns:
        float: The prior value for the specified parameter.
    """
    _selected_row = priors_df[priors_df['Parameter'] == param]
    if _selected_row.empty:
        raise ValueError(f"Parameter '{param}' not found in priors DataFrame.")
    if _selected_row['Distribution'].iloc[0].lower() == 'gaussian':
        value = _selected_row['Input_A'].iloc[0]
    else:
        value = np.mean(_selected_row['Input_A'].iloc[0], _selected_row['Input_B'].iloc[0])

    return value

def read_input_data(input_data):
    """
    Read input data and priors from CSV files.
    
    Args:
        input_data (str): Path to the input data CSV file.
    
    Returns:
        tuple: time, flux, flux_err.
    """
    input_data_df = pd.read_csv(input_data)

    times = []
    fluxes = []
    errors = []

    for i in range(len(input_data_df)):
        _t, _f, _e = pd.read_csv(input_data_df['Path'][i]).values.T
        times.append(_t)
        fluxes.append(_f)
        errors.append(_e)
    return times, fluxes, errors

def inject_ttv(times, P, t0, p_prime, p_dprime):
    """Inject TTV signals into the provided light curves.
    Args:
        times (list): List of time arrays for each light curve.
        P (float): Orbital period.
        t0 (float): Reference time of conjunction.
        p_prime (float): First order TTV parameter.
        p_dprime (float): Second order TTV parameter.
    Returns:
        list: List of time arrays with TTV signals injected.
    """
    times_last = np.empty(0)
    times_first = np.empty(0)

    for i in range(len(times)):
        times_last = np.append(times_last, times[i][-1])
        times_first = np.append(times_first, times[i][0])