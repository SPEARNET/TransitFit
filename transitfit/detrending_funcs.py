'''
Some detrending functions for use with LightCurve
'''

import numpy as np


class NthOrderDetrendingFunction:
    '''
    Arbitrary order detrending function which conserves flux.

    Parameters
    ----------
    order : int
        The detrending function order. Must be greater than 0.
    '''
    def __init__(self, order):

        order = int(order)
        if not order > 0:
            raise ValueError('Order must be greater than 0')

        self.order = order


    def __call__(self, lightcurve, *args):
        '''
        Detrends the lightcurve and returns the flux values
        '''
        
        if not int(len(args)-2) == int(self.order):
            
            raise ValueError('Number of arguments {} supplied does not match the order {} of the detrending function.'.format(len(args)-2, self.order))

        
        t0=args[-2]
        p=args[-1]
        # We calculate t01 which is time of conjucntion for the particular lightcurve. 
        t01=np.max(lightcurve.times)-((np.max(lightcurve.times)-t0)%p)
        #times = lightcurve.times-np.mean(lightcurve.times)
        times = lightcurve.times-t01

        #times = lightcurve.times - np.floor((np.min(lightcurve.times)))
        #self.time0= time0
        
        vals = np.zeros(len(times))
        for i in range(0, self.order):
            #vals += args[i] * (times ** (i+1) - np.mean(times ** (i+1))) # Original function in JH version
            vals += args[i] * ((times)** (i+1))
            #vals += args[i] * ((times)** (i+1))

        return lightcurve.flux - vals
        """
        '''
        Detrends the lightcurve and returns the flux values
        '''
        if not len(args) == self.order:
            raise ValueError('Number of arguments {} supplied does not match the order {} of the detrending function.'.format(len(args), self.order))

        # We rescale the times because in BJD the numbers can get huge!
        times = lightcurve.times - np.floor((np.min(lightcurve.times)))

        #self.time0= time0

        vals = np.zeros(len(times))
        for i in range(0, self.order):
            #vals += args[i] * (times ** (i+1) - np.mean(times ** (i+1))) # Original function in JH version
            vals += args[i] * ((times - np.mean(times))** (i+1)) # We need to replace times with times0

        return lightcurve.flux - vals"""
