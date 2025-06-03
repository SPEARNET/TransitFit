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