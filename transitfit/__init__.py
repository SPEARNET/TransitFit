'''
TransitFit package

This package is designed to fit transit light curves using BATMAN
'''
name = 'transitfit'
__version__ = '4.2.5'

from .retriever import Retriever
from .priorinfo import PriorInfo
from .lightcurve import LightCurve
from ._pipeline import run_retrieval
from ._utils import split_lightcurve_file, calculate_logg, AU_to_host_radii
from .lightcurve import LightCurve
from .error_analysis import get_quantiles_on_best_val_unweighted
