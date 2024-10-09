'''
Module to calculate the likelihood of a set of parameters
'''

import numpy as np
import batman
from copy import deepcopy
from ._paramarray import ParamArray
from .ttv_fitting import taylor_series, get_time_duration, get_total_shift, get_shift_in_time_due_to_ttv
#from collections.abc import Iterable



class LikelihoodCalculator:
    '''
    Object to quickly calculate the likelihood of a set of parameters to
    fit a given set of light curves.

    Parameters
    ----------
    lightcurves : array_like, shape (n_telescopes, n_filters, n_epochs)
        An array of LightCurves. If no data exists for a point in the array
        then the entry should be `None`.
    priors : PriorInfo
        The PriorInfo object for retrieval
    fit_ttv_taylor: bool, optional
        If True, will fit the TTVs using a Taylor expansion of the ephemeris
        equation. Default is False.
    '''
    def __init__(self, lightcurves, priors,fit_ttv_taylor=False,):
        lightcurves = deepcopy(lightcurves)
        self.lightcurves = np.array(lightcurves, dtype=object)

        self.n_telescopes = self.lightcurves.shape[0]
        self.n_filters = self.lightcurves.shape[1]
        self.n_epochs = self.lightcurves.shape[2]

        self.num_light_curves = len(np.where(self.lightcurves.flatten() != None)[0])

        self.priors = priors
        self.fit_ttv_taylor=fit_ttv_taylor

        # We need to make a separate TransitParams and TransitModels for each
        # light curve.

        # Initialise them:
        self.batman_params = ParamArray('batman_params', (self.n_telescopes, self.n_filters, self.n_epochs), True, True, True, lightcurves=self.lightcurves)
        self.batman_models = ParamArray('batman_models', (self.n_telescopes, self.n_filters, self.n_epochs), True, True, True, lightcurves=self.lightcurves)


        for i in np.ndindex(self.lightcurves.shape):
            tidx, fidx, eidx = i

            if self.lightcurves[i] is not None:
                # Set up the params
                self.batman_params.set_value(batman.TransitParams(), tidx, fidx, eidx)

                # Set up the TransitModels
                # Make some realistic parameters to setup the models with
                default_params = batman.TransitParams()
                if self.priors.allow_ttv:
                    default_params.t0 = priors.priors['t0'].default_value
                else:
                    default_params.t0 = priors.priors['t0'].default_value
                default_params.per = priors.priors['P'].default_value
                default_params.rp = priors.priors['rp'].default_value
                default_params.a = priors.priors['a'].default_value
                default_params.inc = priors.priors['inc'].default_value
                default_params.ecc = priors.priors['ecc'].default_value
                default_params.w = priors.priors['w'].default_value
                default_params.limb_dark = priors.limb_dark
                # Note that technically this is q, not u, but it doesn't matter
                # as we are only initialising here
                default_params.u = [priors.priors[qX].default_value for qX in priors.limb_dark_coeffs]

                # Now make the models
                model = batman.TransitModel(default_params, self.lightcurves[i].times)
                self.batman_models.set_value(model, tidx, fidx, eidx)
    
    # New function to handle multiple processes in the same batch.  
    def find_likelihood_parallel_processed(self, cube):
            priors = self.priors
            params = priors._interpret_param_array(cube)

            ln_likelihood = self.find_likelihood(params)
            
            if priors.fit_ld and priors.ld_fit_method not in ["independent","custom"]:
                # Pull out the q values and convert them
                u = []
                for fi in range(priors.n_filters):
                    q = [params[qX][0, fi, 0] for qX in priors.limb_dark_coeffs]
                    u.append(priors.ld_handler.convert_qtou(*q))
                return ln_likelihood + priors.ld_handler.ldtk_lnlike(
                    u, priors.limb_dark
                )
            else:
                return ln_likelihood
            
    def find_period_integration(self,params):
        times_last=np.empty(0)
        times_first=np.empty(0)
        for i in np.ndindex(self.lightcurves.shape):

            if self.lightcurves[i] is not None:
                times_last=np.append(times_last,self.lightcurves[i].times[-1])
                times_first=np.append(times_first,self.lightcurves[i].times[0])

        
        self.P=params['P'][i]
        self.p_prime=params['p_prime'][i]
        try:
            self.p_dprime=params['p_dprime'][i]
        except KeyError:
            self.p_dprime=0
        self.t0=params['t0'][i]

        # We calculate t01 which is time of conjucntion for the first lightcurve. helpful when the given t0 is not the first time of conjuction.
        self.t0_first=np.min(times_last)-((np.min(times_last)-self.t0)%self.P)
        
        self.times_last=times_last#-self.t0_first
        #self.initial_guess_epochs=np.array(self.times_last//self.P,dtype=int)
        #self.initial_guess_epochs-=self.initial_guess_epochs[0]
        #self.initial_guess_epochs+=1

        period_all=np.array([self.P])
        t0_all = np.array([self.t0_first])

        t_start=0
        condition=True
        initial_guess_epochs=np.array((times_last-self.t0_first)//self.P,dtype=int)
        indx=1
        for i in range(1,max(initial_guess_epochs)+1):
            tau=get_time_duration(self.p_prime,self.p_dprime,self.P,t_start)

            # The period at the next epoch
            P_new=taylor_series(self.P,self.p_prime,self.p_dprime,tau)
            self.P=P_new
            t_start+=tau

            if i in initial_guess_epochs:
                if t_start+self.t0_first<times_first[indx] or t_start+self.t0_first>self.times_last[indx]:
                    return None, (max(initial_guess_epochs)-i)
                else:
                    period_all=np.append(period_all,P_new)
                    t0_all=np.append(t0_all,t_start+self.t0_first)
                    indx+=1

        """while condition:
            # The time duration between the current and the next epoch
            tau=get_time_duration(self.p_prime,self.p_dprime,self.P,t_start)

            # The period at the next epoch
            P_new=taylor_series(period_all[-1],self.p_prime,self.p_dprime,tau)
            self.P=P_new

            period_all=np.append(period_all,P_new)

            #period_all=np.append(period_all,taylor_series(self.P, self.p_prime, 0, t_start))
            t_start+=tau
            t0_all=np.append(t0_all,t_start+self.t0_first)
            
            
            if t0_all[-1]>self.times_last[-1]:#+self.t0_first:
                condition=False"""
        #breakpoint()
        #if len(t0_all)<max(self.initial_guess_epochs):
        #   return -1e9

        #period_all=np.append(period_all,taylor_series(self.P, self.p_prime, 0, t_start))
        #period_all=(period_all[:-1]+period_all[1:])/2

        for i in np.ndindex(self.lightcurves.shape):

            if self.lightcurves[i] is not None:
                find_id=t0_all-self.lightcurves[i].times[-1]
                id=np.argmax(find_id[find_id<=0])

                params['P'][i]=period_all[id]
                params['t0'][i]=t0_all[id]
        
        return params
        
            


    def find_period_ttv(self,params):
        times_last=np.empty(0)
        #times_first=np.empty(0)
        for i in np.ndindex(self.lightcurves.shape):

            if self.lightcurves[i] is not None:
                times_last=np.append(times_last,self.lightcurves[i].times[-1])
                #times_first=np.append(times_first,self.lightcurves[i].times[0])

        
        self.P=params['P'][i]
        self.p_prime=params['p_prime'][i]
        #self.p_dprime=params['p_dprime'][i]
        self.t0=params['t0'][i]

        # We calculate t01 which is time of conjucntion for the first lightcurve. helpful when the given t0 is not the first time of conjuction.
        self.t0_first=np.min(times_last)-((np.min(times_last)-self.t0)%self.P)
        
        self.times_last=times_last-self.t0_first
        #self.initial_guess_epochs=np.array(self.times_last//self.P,dtype=int)
        #self.initial_guess_epochs-=self.initial_guess_epochs[0]
        #self.initial_guess_epochs+=1

        period_all=np.empty(0)
        t0_all = np.empty(0)

        t_start=0
        condition=True
        while condition:

            period_all=np.append(period_all,taylor_series(self.P, self.p_prime, 0, t_start))
            t0_all=np.append(t0_all,t_start+self.t0_first)
            t_start+=period_all[-1]
            
            if t0_all[-1]>self.times_last[-1]+self.t0_first:
                condition=False
        #breakpoint()
        #if len(t0_all)<max(self.initial_guess_epochs):
        #   return -1e9

        period_all=np.append(period_all,taylor_series(self.P, self.p_prime, 0, t_start))
        period_all=(period_all[:-1]+period_all[1:])/2

        for i in np.ndindex(self.lightcurves.shape):

            if self.lightcurves[i] is not None:
                find_id=t0_all-self.lightcurves[i].times[-1]
                id=np.argmax(find_id[find_id<=0])

                params['P'][i]=period_all[id]
                params['t0'][i]=t0_all[id]
        
        return params

    def find_likelihood(self, params):
        '''
        Finds the likelihood of a set of parameters
        '''
        all_chi2 = []
        n_data_points = 0
        total_chi2 = 0

        if self.fit_ttv_taylor:
            #params = self.find_period_ttv(params)
            response=self.find_period_integration(params)
            if type(response)==tuple:
                return -1e5*response[1]
            else:
                params=response

        for i in np.ndindex(self.lightcurves.shape):
            tidx, fidx, eidx = i

            if self.lightcurves[i] is not None:

                # GENERATE THE MODEL LIGHT CURVE

                # Convert the LDC q values to u:
                u = self.priors.ld_handler.convert_qtou(*[params[qX][i] for qX in self.priors.limb_dark_coeffs])

                # Need to update the parameters
                self.update_params(tidx, fidx, eidx,
                                   params['t0'][i],
                                   params['P'][i],
                                   params['rp'][i],
                                   params['a'][i],
                                   params['inc'][i],
                                   params['ecc'][i],
                                   params['w'][i],
                                   self.priors.limb_dark,
                                   u)

                # Now we calculate the model transits
                model = self.batman_models[i]
                # Batman considers uniform period. So we need to shift the timestamps accordingly.
                if self.fit_ttv_taylor:
                    _time=model.t
                    shift=get_shift_in_time_due_to_ttv(0,_time-params['t0'][i],params['p_prime'][i],params['p_dprime'][i],params['P'][i], params['t0'][i]-self.t0_first)                    
                    #fraction=(_time-model.t0)/self.P
                    #shift=fraction*((params['p_prime'][i]/2 - 1) * fraction*_time +params['p_dprime'][i]/6 * (fraction*_time)**2)
                    #model.t=_time-shift
                model_flux = model.light_curve(self.batman_params[i])

                # DETREND AND NORMALISE THE DATA TO COMPARE TO THE MODEL
                if self.priors.detrend:
                    # We have to work out the detrending info
                    # Get the method index
                    detrending_index = self.priors._detrend_method_index_array[i]
                    # Now get the parameter values
                    d = [params[di][i] for di in self.priors.detrending_coeffs[detrending_index]]
                    
                else:
                    d = None
                                    
                try:
                    d_new=[*d,params['t0'][i],params['P'][i]]
                except TypeError:
                    d_new=[d,params['t0'][i],params['P'][i]]

                if self.priors.normalise:
                    norm = params['norm'][i]
                else:
                    norm = 1
                detrended_flux, err = self.lightcurves[i].detrend_flux(d_new, norm)
                #detrended_flux, err = self.lightcurves[i].detrend_flux(d, norm)

                #print('Likelihood function detrended flux:')
                #print(detrended_flux)

                # Check that there is actually a transit in the model
                # otherwise we impose a large penalty to the chi2 value
                # This avoids a situation where the models try
                # to completely flatten the light curves, which is wrong!
                if np.isclose(model_flux, 1, atol=0., rtol = 1.e-7).all():
                    total_chi2 += 1e7
                    # Ensures that transit models with depths down to 1e-7 are not penalised. 
                    # This allows for the fact that folded data could be quite sensitive!
                else:
                    # Work out the chi2
                    #chi2 = sum((model_flux - detrended_flux)**2 / err**2)

                    if self.priors.error_scaling:
                        # Scaling the errorbars. https://emcee.readthedocs.io/en/stable/tutorials/line/
                        #if isinstance(params['escale'], Iterable):
                        fn = params['escale'][i]
                        sn_sq = np.power(err,2) + np.power(fn*model_flux,2)
                        total_chi2 += np.sum((np.power((model_flux - detrended_flux),2) / sn_sq) + np.log(2*np.pi*sn_sq))
                    
                    else:
                        total_chi2 += np.sum((model_flux - detrended_flux)**2 / err**2)

                #all_chi2.append(chi2)

        #print(params)
        #print(- total_chi2)
        #return - sum(all_chi2)
        return - total_chi2

    def update_params(self, telescope_idx, filter_idx, epoch_index, t0=None,
                      P=None, rp=None, a=None, inc=None, ecc=None, w=None,
                      limb_dark=None, u=None):
        '''
        Updates self.params with values given
        '''
        if t0 is not None:
            self.batman_params.array[telescope_idx, filter_idx, epoch_index].t0 = t0
        if P is not None:
            self.batman_params.array[telescope_idx, filter_idx, epoch_index].per = P
        if rp is not None:
            self.batman_params.array[telescope_idx, filter_idx, epoch_index].rp = rp
        if a is not None:
            self.batman_params.array[telescope_idx, filter_idx, epoch_index].a = a
        if inc is not None:
            self.batman_params.array[telescope_idx, filter_idx, epoch_index].inc = inc
        if ecc is not None:
            self.batman_params.array[telescope_idx, filter_idx, epoch_index].ecc = ecc
        if w is not None:
            self.batman_params.array[telescope_idx, filter_idx, epoch_index].w = w
        if limb_dark is not None:
            self.batman_params.array[telescope_idx, filter_idx, epoch_index].limb_dark = limb_dark
        if u is not None:
            self.batman_params.array[telescope_idx, filter_idx, epoch_index].u = u
