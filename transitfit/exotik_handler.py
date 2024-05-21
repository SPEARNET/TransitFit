import numpy as np
from scipy.special import erf, erfinv
from collections.abc import Iterable
import pandas as pd
from exotic_ld import StellarLimbDarkening


ld_data_path = '/home/a/Documents/GitHub/exotic'

def get_average(results):
    u0,u0_err,u1,u1_err=results#ldc_calc.get_ld()

    return np.average(u0,weights=1/u0_err**2),1/np.sqrt(np.sum(1/np.power(u0_err,2))),np.average(u1,weights=1/u1_err**2),1/np.sqrt(np.sum(1/np.power(u1_err,2)))
    


def gaussian_clipped(mean, stdev, min_value, max_value, n_samples):

    uniform_samples=np.random.uniform(0,1,n_samples)
    lims=(erf((min_value - mean)/(stdev * np.sqrt(2))) + 1)/2 ,(erf((max_value - mean)/(stdev * np.sqrt(2))) + 1)/2
    uniform_samples_clipped = uniform_samples * (lims[1] - lims[0]) + lims[0]
    gaussian_samples = mean + stdev * np.sqrt(2) * erfinv(2 * uniform_samples_clipped - 1)

    return gaussian_samples


class LDC:
    def __init__(self,teff=4400.0, logg=4.49, FeH=-0.01,model='quadratic',nsamples=10000,atmospheric_model='ATLAS9',**kwargs):

        if isinstance(teff, Iterable) and nsamples>1:
            if len(teff)==2:
    
                if atmospheric_model=='ATLAS9':
                    #ATLAS9
                    # https://exotic-ld.readthedocs.io/en/latest/views/supported_stellar_grids.html
                    min_teff,max_teff=3500,6500
                    min_logg,max_logg=4.0,5.0
                    min_FeH,max_FeH=-5,1

                
                if teff[0]<min_teff or teff[0]>max_teff:
                    raise ValueError('Teff out of bounds')
                if logg[0]<min_logg or logg[0]>max_logg:
                    raise ValueError('logg out of bounds')
                if FeH[0]<min_FeH or FeH[0]>max_FeH:
                    raise ValueError('FeH out of bounds')

                self.teff = gaussian_clipped(*teff,min_value=min_teff,max_value=max_teff,n_samples=nsamples)
                self.logg = gaussian_clipped(*logg,min_value=min_logg,max_value=max_logg,n_samples=nsamples)
                self.FeH = gaussian_clipped(*FeH,min_value=min_FeH,max_value=max_FeH,n_samples=nsamples)
                self.nsamples=nsamples
        elif isinstance(teff, Iterable):
            self.teff = teff[0]
            self.logg = logg[0]
            self.FeH = FeH[0]
            self.nsamples=nsamples
        else:
            self.teff = teff
            self.logg = logg
            self.FeH = FeH
            self.nsamples=1
        
        self.model = model

        

    def get_ld(self,wvs, throughput):
        if self.nsamples>1:
            u0s=np.empty(self.nsamples)
            u1s=np.empty(self.nsamples)
            for i in range(self.nsamples):
                sld = StellarLimbDarkening(self.FeH[i], self.teff[i], self.logg[i], ld_model = 'kurucz', ld_data_path = ld_data_path, interpolate_type='trilinear')
                u0s[i],u1s[i] = sld.compute_quadratic_ld_coeffs(wavelength_range=[min(wvs), max(wvs)],mode="custom",custom_wavelengths=wvs,custom_throughput=throughput)

            
            return u0s,u1s

        else:
            sld = StellarLimbDarkening(self.FeH, self.teff, self.logg, ld_model = 'kurucz', ld_data_path = ld_data_path, interpolate_type='trilinear')
            return sld.compute_quadratic_ld_coeffs(wavelength_range=[min(wvs), max(wvs)],mode="custom",custom_wavelengths=wvs,custom_throughput=throughput)
        

def change_priors(_filter_input,filters,host_T,host_logg, host_z, n_ld_samples,ldtk_uncertainty_multiplier,priors):
    from pathlib import Path
    _path = Path(_filter_input)
    _path=_path.parent.absolute()
    dfp=pd.read_csv(priors)

    for i,fil in enumerate(filters):
        
        ldc_calc = LDC(teff=host_T, logg=host_logg, FeH=host_z,nsamples=n_ld_samples,n_bins=1,atmospheric_model='ATLAS9')
        u0, u1=ldc_calc.get_ld(wvs=fil[0]*10, throughput=fil[1])
        u0, u0_err, u1, u1_err=np.average(u0),np.std(u0),np.average(u1),np.std(u1)

        dfp.loc[len(dfp)] = ['u0', 'gaussian', u0, u0_err, int(i)]
        dfp.loc[len(dfp)] = ['u1', 'gaussian', u1, u1_err, int(i)]

    _priors = Path(priors)
    _priors=_priors.parent.absolute()
    _priors=str(_priors)+'/_priors_with_custom_ldc.csv'
    Path(_priors).unlink(missing_ok=True)
    dfp.to_csv(_priors,index=None)
    return _priors