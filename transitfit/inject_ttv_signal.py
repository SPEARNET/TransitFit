# Inject TTV signals
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time
import os
from .ttv_fitting import get_shift_in_time_due_to_ttv,get_time_duration, taylor_series