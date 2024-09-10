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

def get_time_duration(P_prime,P_dprime,P):
    a=P_dprime/6
    b=(P_prime/2 -1)
    c=P+0

    if a==0:
        tau=-c/b
    else:
        d=b**2-4*a*c
        if d<0:
            print("no real solution")
            return None
        else:
            x1=(-b+np.sqrt(d))/(2*a)
            x2=(-b-np.sqrt(d))/(2*a)
            if np.abs(x1-P)<np.abs(x2-P):
                tau=x1
            else:
                tau=x2
    return tau


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