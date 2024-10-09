import numpy as np
import matplotlib.pyplot as plt
from dynesty import NestedSampler
import multiprocessing as mp
from .count_params import count_number_lcs
import batman
import pandas as pd
from .io import (
    read_input_file,
    read_priors_file,
    parse_priors_list,
    read_filter_info,
    parse_filter_list,
    print_results,
)

def get_time_duration(P_prime,P_dprime,P,t0=0):
    """Gets the time duration between two epochs

    Args:
        P_prime (float): the first derivative of the period
        P_dprime (float): second derivative of the period
        P (float): Period
        t0 (float, optional): time of conjunction. Defaults to 0.

    Returns:
        float: time duration between two epochs
    """
    a=np.float64(P_dprime/6)
    b=np.float64((P_prime/2 -1))
    c=P+0
    correction_term=np.float64(P_dprime)*t0/2
    b-=correction_term

    if a==0:
        tau=np.float64(-c/b)
    else:
        d=np.float64(b**2-4*a*c)
        if d<0:
            print("no real solution")
            return None
        else:
            x1=np.float64((-b+np.sqrt(d))/(2*a))
            x2=np.float64((-b-np.sqrt(d))/(2*a))
            if np.abs(x1-P)<np.abs(x2-P):
                tau=x1
            else:
                tau=x2
    return tau

def get_val(x,p_prime,p_dprime,P,t0):
    """value of the function at x

    Args:
        x (float): timestamp

    Returns:
        float: value of the function at x
    """
    correction_term=np.float64(p_dprime)*t0/2
    return np.float64(P*x+((p_prime/2)+correction_term)*x**2+p_dprime/6*x**3)

# Average of period from t2 to t1
def get_avg_period(t1,t2,p_prime,p_dprime,P,t0):
    """calculates average of period from t1 to t2

    Args:
        t1 (float): the first timestamp
        t2 (float): the second timestamp

    Returns:
        float: the average of period from t1 to t2
    """
    return np.float64((1/(t2-t1))*(get_val(t2,p_prime,p_dprime,P,t0)-get_val(t1,p_prime,p_dprime,P,t0)))

def get_total_shift(t1,t2,p_prime,p_dprime,P,t0):
    """calculates the total shift due to TTV between t1 and t2

    Args:
        t1 (float): the first timestamp
        t2 (float): the second timestamp

    Returns:
        float: the total shift between t1 and t2
    """
    tau=get_time_duration(p_prime,p_dprime,P,t0)
    tau_previous=get_time_duration(-p_prime,p_dprime,P,t0)

    if type(t2)==np.ndarray:
        tau_array=np.zeros(len(t2))
        #idxs=np.where(t2>t1)
        tau_array[t2>t1]=tau
        tau_array[t2<=t1]=tau_previous
        tau=tau_array
        
    return np.float64((get_avg_period(t1,t2,p_prime,p_dprime,P,t0)-P)*(t2-t1)/tau)

def get_integral_at_x(x,p_prime,p_dprime,t0):
    """value of the function at x

    Args:
        x (float): timestamp

    Returns:
        float: value of the function at x
    """
    correction_term=np.float64(p_dprime)*t0*(x**2)/2
    
    term_1=np.float64((p_prime*x**2)/2)
    term_2=np.float64((p_dprime*x**3)/6)
    correction_term=np.float64(p_dprime*t0*x**2)/2

    return term_1+term_2+correction_term

def get_shift_in_time_due_to_ttv(t1,t2,p_prime,p_dprime,P,t0):
    tau=get_time_duration(p_prime,p_dprime,P,t0)
    tau_previous=get_time_duration(-p_prime,p_dprime,P,t0)

    if type(t2)==np.ndarray:
        tau_array=np.zeros(len(t2))
        #idxs=np.where(t2>t1)
        tau_array[t2>t1]=tau
        tau_array[t2<=t1]=tau_previous
        tau=tau_array

    return (1/tau)*(get_integral_at_x(t2,p_prime,p_dprime,t0)-get_integral_at_x(t1,p_prime,p_dprime,t0))



def taylor_series(period, period_prime, period_dprime, times):
    """Calculates the period for given times using the taylor series expansion.

    Args:
        period (float): period
        period_prime (float): first derivative of period wrt time
        period_dprime (float): second derivative of period wrt time
        t0_conjunction (float): time of conjunction
        times (array): times noted in the lightcurve

    Returns:
        float: the period as a summation of Taylor's series
    """

    period_all = (
        period + (times * period_prime) + (np.power(times, 2) * period_dprime / 2)
    )

    return period_all