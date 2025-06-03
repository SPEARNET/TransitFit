'''
A class for parameters in TransitFit which can be retrieved. These are used
by the PriorInfo to determine dimensionality etc.

'''
from scipy.special import erf, erfinv
import numpy as np

class _Param:
    def __init__(self, value, uncertainty=None):

        self.default_value = value
        self.low_lim=None
        self.high_lim=None
        self.uncertainty=uncertainty

    def from_unit_interval(self, u):
        raise NotImplementedError


class _UniformParam(_Param):
    """Class to represent a uniform parameter in a given range."""
    def __init__(self, low_lim, high_lim, negative_allowed=True):
        if low_lim >= high_lim:
            raise ValueError('low_lim >= high_lim')

        super().__init__((high_lim + low_lim)/2)
        self.low_lim = low_lim
        self.high_lim = high_lim
        self.negative_allowed = negative_allowed

    def from_unit_interval(self, u):
        '''
        Function to convert value u in range (0,1], will convert to a value to
        be used by Batman
        '''
        if u > 1 or u < 0:
            raise ValueError('u must satisfy 0 < u < 1. ')
        val = u * (self.high_lim - self.low_lim) + self.low_lim
        if self.negative_allowed:
            return val
        return abs(val)

class _GaussianParam(_Param):
    def __init__(self, best, sigma, negative_allowed=True, clipped_gaussian=False, custom_ldcs=False):
        '''
        A GaussianParam is one which is fitted using a Gaussian prior (normal)
        distribution.
        '''

        super().__init__(best)
        self.mean = best
        self.stdev = sigma
        self.negative_allowed = negative_allowed
        self.custom_ldcs=custom_ldcs

        if clipped_gaussian:
            self.factor_for_clipping=(erf((90 - self.mean)/(self.stdev * np.sqrt(2))) + 1)/2 
        
        if custom_ldcs:
            self.ldcs_lims=(erf((0 - self.mean)/(self.stdev * np.sqrt(2))) + 1)/2 ,(erf((1 - self.mean)/(self.stdev * np.sqrt(2))) + 1)/2
        

    def from_unit_interval(self, u):
        '''
        Function to convert value u in range (0,1], will convert to a value to
        be used by Batman
        '''
        if u > 1 or u < 0:
            raise ValueError('u must satisfy 0 < u < 1')

        if self.custom_ldcs:
            u = u * (self.ldcs_lims[1] - self.ldcs_lims[0]) + self.ldcs_lims[0]

        val = self.mean + self.stdev * np.sqrt(2) * erfinv(2 * u - 1)
        if self.negative_allowed:
            return val
        return abs(val)
    
    def uniform_to_clipped_gaussian(self,u):
        '''
        Function to convert value u in range (0,1], will convert to a value to
        be used by Batman
        '''
        if u > 1 or u < 0:
            raise ValueError('u must satisfy 0 < u < 1')

        val = self.mean + self.stdev * np.sqrt(2) * erfinv(2 * u * self.factor_for_clipping - 1)
        if self.negative_allowed:
            return val
        return abs(val)