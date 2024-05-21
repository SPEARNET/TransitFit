from exoctk import modelgrid
from exoctk.limb_darkening import limb_darkening_fit as lf
from exoctk.throughputs import Throughput
import numpy as np
from scipy.special import erf, erfinv
from collections.abc import Iterable
import pandas as pd


# Bootstrapped
def get_average(results):
    u0,u0_err,u1,u1_err=results#ldc_calc.get_ld()

    return np.average(u0,weights=1/u0_err**2),1/np.sqrt(np.sum(1/np.power(u0_err,2))),np.average(u1,weights=1/u1_err**2),1/np.sqrt(np.sum(1/np.power(u1_err,2)))
    
def extreme_samples_clipped(mean, stdev, min_value, max_value, n_samples):
    
    samples=np.array([max(min_value,mean-stdev),min(max_value,mean+stdev)])
    
    return samples

def get_extreme(results):
    u0,u0_err,u1,u1_err=results#ldc_calc.get_ld()

    return (max(u0)+min(u0))/2,(max(u0)-min(u0))/2,(max(u1)+min(u1))/2,(max(u1)-min(u1))/2
    


def gaussian_clipped(mean, stdev, min_value, max_value, n_samples):

    uniform_samples=np.random.uniform(0,1,n_samples)
    lims=(erf((min_value - mean)/(stdev * np.sqrt(2))) + 1)/2 ,(erf((max_value - mean)/(stdev * np.sqrt(2))) + 1)/2
    uniform_samples_clipped = uniform_samples * (lims[1] - lims[0]) + lims[0]
    gaussian_samples = mean + stdev * np.sqrt(2) * erfinv(2 * uniform_samples_clipped - 1)

    return gaussian_samples


class LDC:
    def __init__(self, filter_name='MIRI.CLEAR.P750L.LRSSLIT', teff=4400.0, logg=4.49, FeH=-0.01,model='quadratic',nsamples=10000,atmospheric_model='ATLAS9',**kwargs):
        self.filter_name = filter_name
        if isinstance(teff, Iterable) and nsamples>1:
            if len(teff)==2:
    
                if atmospheric_model=='ATLAS9':
                    #ATLAS9
                    min_teff,max_teff=3500,8750
                    min_logg,max_logg=3.0,5.0
                    min_FeH,max_FeH=-0.5,0.5


                else:#=='ACES':
                    #ACES https://www.aanda.org/articles/aa/full_html/2013/05/aa19058-12/aa19058-12.html
                    # wlen: [500 A, 5.5um]
                    min_teff,max_teff=2300,7800 # 12000
                    min_logg,max_logg=0.0,6.0
                    min_FeH,max_FeH=-4.0,1.0
                
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

        if atmospheric_model=='ATLAS9':
            model_grid = modelgrid.ATLAS9(processes=8, interp=True)
            self.ld = lf.LDC(model_grid,interp=True)
        else:#=='ACES':
            model_grid = modelgrid.ACES(processes=8, interp=False)
            self.ld = lf.LDC(model_grid,interp=False)

        self.required_filter=Throughput(self.filter_name, **kwargs)

    def get_ld(self,):
        if self.nsamples>1:
            u0s=np.empty(self.nsamples)
            u0errs=np.empty(self.nsamples)
            u1s=np.empty(self.nsamples)
            u1errs=np.empty(self.nsamples)
            for i in range(self.nsamples):
                results=self.ld.calculate(self.teff[i], self.logg[i], self.FeH[i], self.model, bandpass=self.required_filter)
                u0s[i],u1s[i]=np.array(results['coeffs'])
                u0errs[i],u1errs[i]=np.array(results['errors'])
            
            return u0s,u0errs,u1s,u1errs

        else:
            results = self.ld.calculate(self.teff, self.logg, self.FeH, self.model, bandpass=self.required_filter)

            coeffs=np.array(results['coeffs'])
            errs=np.array(results['errors'])
        
            return (coeffs[0],errs[0],coeffs[1],errs[1])

class LDC_extremes:
    def __init__(self, filter_name='MIRI.CLEAR.P750L.LRSSLIT', teff=4400.0, logg=4.49, FeH=-0.01,model='quadratic',nsamples=8,atmospheric_model='ATLAS9',**kwargs):
        self.filter_name = filter_name
        if isinstance(teff, Iterable) and nsamples>1:
            if len(teff)==2:
    
                if atmospheric_model=='ATLAS9':
                    #ATLAS9
                    min_teff,max_teff=3500,8750
                    min_logg,max_logg=3.0,5.0
                    min_FeH,max_FeH=-0.5,0.5


                else:#=='ACES':
                    #ACES https://www.aanda.org/articles/aa/full_html/2013/05/aa19058-12/aa19058-12.html
                    # wlen: [500 A, 5.5um]
                    min_teff,max_teff=2300,7800 # 12000
                    min_logg,max_logg=0.0,6.0
                    min_FeH,max_FeH=-4.0,1.0
                
                if teff[0]<min_teff or teff[0]>max_teff:
                    raise ValueError('Teff out of bounds')
                if logg[0]<min_logg or logg[0]>max_logg:
                    raise ValueError('logg out of bounds')
                if FeH[0]<min_FeH or FeH[0]>max_FeH:
                    raise ValueError('FeH out of bounds')

                self.teff = extreme_samples_clipped(*teff,min_value=min_teff,max_value=max_teff,n_samples=nsamples)
                self.logg = extreme_samples_clipped(*logg,min_value=min_logg,max_value=max_logg,n_samples=nsamples)
                self.FeH = extreme_samples_clipped(*FeH,min_value=min_FeH,max_value=max_FeH,n_samples=nsamples)
                self.nsamples=8
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

        if atmospheric_model=='ATLAS9':
            model_grid = modelgrid.ATLAS9(processes=8, interp=True)
            self.ld = lf.LDC(model_grid,interp=True)
        else:#=='ACES':
            model_grid = modelgrid.ACES(processes=8, interp=False)
            self.ld = lf.LDC(model_grid,interp=False)

        self.required_filter=Throughput(self.filter_name, **kwargs)

    def get_ld(self,):
        if self.nsamples>1:
            self.nsamples=8
            u0s=np.empty(self.nsamples)
            u0errs=np.empty(self.nsamples)
            u1s=np.empty(self.nsamples)
            u1errs=np.empty(self.nsamples)
            i=0
            for teff in self.teff:
                for logg in self.logg:
                    for FeH in self.FeH:
                        results=self.ld.calculate(teff, logg, FeH, self.model, bandpass=self.required_filter)
                        u0s[i],u1s[i]=np.array(results['coeffs'])
                        u0errs[i],u1errs[i]=np.array(results['errors'])
                        i+=1
                            
            return u0s,u0errs,u1s,u1errs

        else:
            results = self.ld.calculate(self.teff, self.logg, self.FeH, self.model, bandpass=self.required_filter)

            coeffs=np.array(results['coeffs'])
            errs=np.array(results['errors'])
        
            return (coeffs[0],errs[0],coeffs[1],errs[1])
            

def change_priors(_filter_input,filters,host_T,host_logg, host_z, n_ld_samples,ldtk_uncertainty_multiplier,priors):
    from pathlib import Path
    _path = Path(_filter_input)
    _path=_path.parent.absolute()
    dfp=pd.read_csv(priors)

    for i,fil in enumerate(filters):
        df = pd.DataFrame({'Column1': fil[0]/1000, 'Column2': fil[1]})
        _tempname=str(_path)+'/_'+str(i)+'.txt'
        Path(_tempname).unlink(missing_ok=True)

        df.to_csv(_tempname, header=None, index=None, sep=' ', mode='a')

        ldc_calc = LDC(filter_name=_tempname, teff=host_T, logg=host_logg, FeH=host_z,nsamples=n_ld_samples,n_bins=1,atmospheric_model='ATLAS9')
        results=ldc_calc.get_ld()
        u0, u0_err, u1, u1_err=get_average(results)
        #print(u0, u0_err*ldtk_uncertainty_multiplier, u1, u1_err*ldtk_uncertainty_multiplier)
        Path(_tempname).unlink(missing_ok=True)

        dfp.loc[len(dfp)] = ['u0', 'gaussian', u0, u0_err, int(i)]
        dfp.loc[len(dfp)] = ['u1', 'gaussian', u1, u1_err, int(i)]

    _priors = Path(priors)
    _priors=_priors.parent.absolute()
    _priors=str(_priors)+'/_priors_with_custom_ldc.csv'
    Path(_priors).unlink(missing_ok=True)
    dfp.to_csv(_priors,index=None)
    return _priors


def change_priors_take_extremes(_filter_input,filters,host_T,host_logg, host_z, n_ld_samples,ldtk_uncertainty_multiplier,priors):
    from pathlib import Path
    _path = Path(_filter_input)
    _path=_path.parent.absolute()
    dfp=pd.read_csv(priors)
    dfp = dfp[~dfp['Parameter'].isin(['u0', 'u1', 'q0', 'q1'])]

    for i,fil in enumerate(filters):
        df = pd.DataFrame({'Column1': fil[0]/1000, 'Column2': fil[1]})
        _tempname=str(_path)+'/_'+str(i)+'.txt'
        Path(_tempname).unlink(missing_ok=True)

        df.to_csv(_tempname, header=None, index=None, sep=' ', mode='a')

        ldc_calc = LDC_extremes(filter_name=_tempname, teff=host_T, logg=host_logg, FeH=host_z,nsamples=n_ld_samples,n_bins=1,atmospheric_model='ATLAS9')
        results=ldc_calc.get_ld()
        u0, u0_err, u1, u1_err=get_extreme(results)
        #print(u0, u0_err*ldtk_uncertainty_multiplier, u1, u1_err*ldtk_uncertainty_multiplier)
        Path(_tempname).unlink(missing_ok=True)

        dfp.loc[len(dfp)] = ['u0', 'gaussian', u0, u0_err, int(i)]
        dfp.loc[len(dfp)] = ['u1', 'gaussian', u1, u1_err, int(i)]

    _priors = Path(priors)
    _priors=_priors.parent.absolute()
    _priors=str(_priors)+'/_priors_with_custom_ldc.csv'
    Path(_priors).unlink(missing_ok=True)
    dfp.to_csv(_priors,index=None)
    return _priors