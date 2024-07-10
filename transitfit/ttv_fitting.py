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