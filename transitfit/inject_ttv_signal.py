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
    
    # We calculate t01 which is time of conjucntion for the first lightcurve. helpful when the given t0 is not the first time of conjuction.
    t0_first=np.min(times_last)-((np.min(times_last)-t0)%P)

    period_all=np.array([P])
    t0_all = np.array([t0_first])

    t_start=0
    condition=True
    initial_guess_epochs=np.array((times_last-t0_first)//P,dtype=int)
    indx=1
    _P=P+0
    for i in range(1,max(initial_guess_epochs)+1):
        tau=get_time_duration(p_prime,p_dprime,_P,t_start)

        if tau is None:
            print(None, (max(initial_guess_epochs)-i))
            break

        if tau>2*_P or tau<.5*_P:
            print(None, (max(initial_guess_epochs)-i))
            break

        # The period at the next epoch
        P_new=taylor_series(_P,p_prime,p_dprime,tau,t_start)
        _P=P_new
        t_start+=tau

        if i in initial_guess_epochs:
            if t_start+t0_first<times_first[indx] or t_start+t0_first>times_last[indx]:
                print(None, (max(initial_guess_epochs)-i))
                break
            else:
                period_all=np.append(period_all,P_new)
                t0_all=np.append(t0_all,t_start+t0_first)
                indx+=1
        
    p_list=[]
    t0_list=[]
    for i in range(len(times)):
        find_id=t0_all-times[i][-1]
        id=np.argmax(find_id[find_id<=0])

        p_list.append(period_all[id])
        t0_list.append(t0_all[id])

    times_ttv_injected=[]
    times_ttv_injected_with_shift=[]
    t0_no_ttv=t0+(initial_guess_epochs*P)