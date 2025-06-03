import numpy as np
import glob
import pandas as pd
import pickle
from statsmodels.stats import weightstats

def find_avg_binned_likelihood(x, y, num_bins=100):
    """
    Find the upper envelope points without smoothing.
    
    Parameters:
    x, y: arrays of data points
    num_bins: number of bins to divide x-axis into
    
    Returns:
    bin_centers, max_points: arrays of x and y coordinates for the envelope
    """
    # Create bins along x-axis
    bins = np.linspace(min(x), max(x), num_bins + 1)
    all_li=np.zeros(len(x))

    
    # Find maximum y value in each bin
    max_points = []
    bin_centers = []
    li=[]
    likelihood_new=y-max(y)
    
    for i in range(len(bins)-1):
        mask = (x >= bins[i]) & (x < bins[i+1])
        if np.any(mask):  # Only include bins that contain points
            max_points.append(np.max(y[mask]))
            bin_centers.append((bins[i] + bins[i+1]) / 2)
            li.append(np.sum(y[mask])/len(y[mask]))
            all_li[mask]=np.sum(np.exp(likelihood_new[mask]))/len(likelihood_new[mask])
    
    return all_li

def get_error_from_binned_lkl(chosen_sample, best,logl):
    """Calculates the error from the binned likelihood of the samples.
    Args:
        chosen_sample (array): the sampled values for the parameter from dynesty
        best (float): best value among the samples
        logl (array): log likelihood of the samples
    Returns:
        tuple: the lower and upper error on the best value.
    """
    all_li=find_avg_binned_likelihood(chosen_sample, logl,num_bins=1000)

    err=np.sqrt(np.sum(np.power((chosen_sample-best),2)*all_li)/np.sum(all_li))#*err_weight
    return err, err


def make_dict(val_dict, key, vals, logl=None):
    """makes dictionary using the given values. if the key already 
    exists, then the values get appended. 

    Args:
        val_dict (dict): the dictionary where the values are to be added
        key (str): the key of the value
        vals (array): the values to be added
    """
    if key in val_dict:
        val_dict[key] = np.ndarray.flatten(np.append(val_dict[key], vals))
        if logl is not None:
            val_dict[key+'_logl'] = np.ndarray.flatten(np.append(val_dict[key+'_logl'], logl))
    else:
        val_dict[key] = np.ndarray.flatten(vals)
        if logl is not None:
            val_dict[key+'_logl'] = np.ndarray.flatten(logl)

    
    return val_dict

def check_index(index):
    """Checks if index is number or None. Returns '-' if None."""
    if index=='-':
        return index
    try:
        return str(int(index))
    except:
        return '-'
    return str(index)

def get_params(results_csv):
    """Gets the parameters from the results csv file.
    Args:
        results_csv (pd.DataFrame): the results csv file containing the parameters, telescopes, filters, epochs, best values and errors.
    Returns:
        tuple: a list of parameters to add and a list of best values.
    """
    params=results_csv['Parameter'].values
    telescope_idxs=results_csv['Telescope'].values
    filter_idxs=results_csv['Filter'].values
    epoch_idxs=results_csv['Epoch'].values
    best=results_csv['Best'].values
    error=results_csv['Error'].values

    params_to_add=[]
    params_best=[]
    for p,par in enumerate(params):
        if error[p]=='-':
            continue
        if par in ['u0','u1','u2','u3','a/AU']:
            continue
        params_to_add.append(par+'_'+check_index(telescope_idxs[p])+'_'+check_index(filter_idxs[p])+'_'+check_index(epoch_idxs[p]))
        params_best.append(best[p])
    
    return params_to_add,params_best

def get_quantiles_on_best_val_unweighted(samples, best_val):
    """Generates lower and upper limit of errors for the best values.
    Gets value of samples such that they encompass 68.27% of the samples on both sides of the best value.

    Args:
        samples (array): the sampled values for the parameter from dynesty
        best_val (float): best value among the samples

    Returns:
        tuple: the lower and upper error on the best value.
    """
    errors=-np.abs(np.percentile(samples[samples<best_val], 31.73)-best_val), np.abs(np.percentile(samples[samples>best_val], 68.27)-best_val)
    
    return errors

def get_quantiles_on_best_val(samples, weights, best_val):
    """Generates lower and upper limit of errors for the best values.
    Sorts the samples and corresponding weights. 
    Gets value of samples such that they encompass 68.27% of the samples
    by weight on both sides of the best value.

    Args:
        samples (array): the sampled values for the parameter from dynesty
        weights (array): weights of the samples
        best_val (float): best value among the samples

    Returns:
        tuple: the lower and upper error on the best value.
    """

    weights = weights/np.sum(weights)
    sorter_sample = np.argsort(samples)

    sorted_weights = weights[sorter_sample]
    sorted_samples = samples[sorter_sample]

    best_val_idx = np.argmin(np.abs(sorted_samples-float(best_val)))

    # Faster method
    best_percentile = np.sum(sorted_weights[:best_val_idx])
    model = weightstats.DescrStatsW(sorted_samples, sorted_weights)
    lower_error, upper_error = model.quantile([(1-.6827)*best_percentile, best_percentile+(.6827*(
        1-best_percentile))], return_pandas=False)-sorted_samples[best_val_idx]

    """
    # Previous method (slower)
    # Shows more detail on how we are getting the errors
    
    total_sum = 0
    for i in range(best_val_idx, len(sorted_weights)):
        total_sum += sorted_weights[i]
        if total_sum >= .6827*np.sum(sorted_weights[best_val_idx:]):
            upper_error = sorted_samples[i]-best_val
            break

    total_sum = 0
    for i in range(best_val_idx, 0, -1):
        total_sum += sorted_weights[i]
        if total_sum >= .6827*np.sum(sorted_weights[:best_val_idx]):
            lower_error = sorted_samples[i]-best_val
            break"""

    return -np.abs(lower_error), np.abs(upper_error)


def get_std_on_best_val_unweighted(samples, best_val):
    """Generates lower and upper limit of errors for the best values.
    Gets value of samples such that they encompass 68.27% of the samples on both sides of the best value.

    Args:
        samples (array): the sampled values for the parameter from dynesty
        best_val (float): best value among the samples

    Returns:
        tuple: the lower and upper error on the best value.
    """
    _e=np.power(np.sum(np.power(samples-best_val,2))/len(samples),.5)
    errors=(_e,_e)
    return errors 

def HST_detrending():
    pass

def check_files_for_samples(pathname_to_check, params_to_add, values):
    files=glob.glob(pathname_to_check)

    if len(files)>0:
        print("filter files detected")
        more_params=params_to_add
        for file in files:

            with open(file, 'rb') as handle:
                results = pickle.load(handle)
                samples = results.samples
                #weights=np.exp(results.logwt - results.logwt.max())/np.sum(np.exp(results.logwt - results.logwt.max()))
                order_of_params=results.fitting_params
            
            for o, order in enumerate(order_of_params):
                if order[0] in ['a','rp']:
                    order[0]+='/r*'
                if order[0] in ['t0']: # t0 is only epoch dependent
                    par=order[0]+'_-_-_'+check_index(order[3])
                elif 'q' in order[0] or 'rp' in order[0]: # rp and LDC are only filter dependent  
                    par=order[0]+'_-_'+check_index(order[2])+'_-'
                else:
                    par=order[0]+'_'+check_index(order[1])+'_'+check_index(order[2])+'_'+check_index(order[3])
                if par in more_params:
                    values=make_dict(values, par, samples[:, o])
    return values

def get_asymmetric_errors_updated(folder):
    if folder[-1]!='/':
        folder+='/'
    values = {}
    complete_results=glob.glob(folder+"Complete_results.csv")
    if len(complete_results)==0:
        raise RuntimeError('Please check the pathname, it may not be correct.')
    else:
        complete_results=pd.read_csv(complete_results[0])
        
    params_to_add,params_best=get_params(complete_results)

    for i, param in enumerate(params_to_add):
        values[param+'_best']=params_best[i]

    check_quicksaves=glob.glob(folder+"quicksaves/*results.pkl")

    if len(check_quicksaves)>0:
        
        for file in check_quicksaves:

            with open(file, 'rb') as handle:
                results = pickle.load(handle)
                samples = results.samples
                logl=results.logl
                weights=np.exp(results.logwt - results.logwt.max())/np.sum(np.exp(results.logwt - results.logwt.max()))
                order_of_params=results.fitting_params
            
            for o, order in enumerate(order_of_params):
                if order[0] in ['a','rp']:
                    order[0]+='/r*'
                if order[0] in ['t0']: # t0 is only epoch dependent
                    par=order[0]+'_-_-_'+check_index(order[3])
                elif 'q' in order[0] or 'rp' in order[0]: # rp and LDC are only filter dependent  
                    par=order[0]+'_-_'+check_index(order[2])+'_-'
                else:
                    par=order[0]+'_'+check_index(order[1])+'_'+check_index(order[2])+'_'+check_index(order[3])
                values=make_dict(values, par, samples[:, o], logl)

    # Check folded mode:
    folded_files=glob.glob(folder+"*_parameters/quicksaves/*results.pkl")

    if len(folded_files)>0:
        print("filter files detected")
        more_params=params_to_add-values.keys()
        for file in folded_files:

            with open(file, 'rb') as handle:
                results = pickle.load(handle)
                samples = results.samples
                logl=results.logl
                weights=np.exp(results.logwt - results.logwt.max())/np.sum(np.exp(results.logwt - results.logwt.max()))
                order_of_params=results.fitting_params
            
            for o, order in enumerate(order_of_params):
                if order[0] in ['a','rp']:
                    order[0]+='/r*'
                if order[0] in ['t0']: # t0 is only epoch dependent
                    par=order[0]+'_-_-_'+check_index(order[3])
                elif 'q' in order[0] or 'rp' in order[0]: # rp and LDC are only filter dependent  
                    par=order[0]+'_-_'+check_index(order[2])+'_-'
                else:
                    par=order[0]+'_'+check_index(order[1])+'_'+check_index(order[2])+'_'+check_index(order[3])
                if par in more_params:
                    values=make_dict(values, par, samples[:, o], logl)
            
    lower_errors=[]
    upper_errors=[]
    params_to_save=[]
    telescopes=[]
    filters=[]
    epochs=[]
    issue_with_priors=[]
    for p in params_to_add:
        samples=values[p]
        best=values[p+'_best']
        logl=values[p+'_logl']
        weights=find_avg_binned_likelihood(samples, logl,num_bins=1000)
        try:
            le,ue=get_quantiles_on_best_val(samples,weights, best)
        except:
            issue_with_priors.append(p.replace('_',', '))
            le,ue=get_error_from_binned_lkl(samples, best,logl)
        lower_errors.append(le)
        upper_errors.append(ue)
        _p=p.split('_')

        params_to_save.append("_".join(_p[0:-3]))
        telescopes.append(_p[-3])
        filters.append(_p[-2])
        epochs.append(_p[-1])

    if len(issue_with_priors)>0:
        print("The priors for following parameters might be too strict. Consider expanding the priors for them: \nParameter, Telescope, Filter, Epoch")
        for i in issue_with_priors:
            print(i)

    data = {'Parameter': params_to_save,'Telescope':telescopes,'Filter':filters,'Epoch':epochs, 'Best': params_best, 'Lower_error': lower_errors, 'Upper_error': upper_errors}
    df = pd.DataFrame(data)
    df.to_csv(folder+'results_with_asymmetric_errors.csv', index=False)


def get_asymmetric_errors(folder):
    if folder[-1]!='/':
        folder+='/'
    values = {}
    complete_results=glob.glob(folder+"Complete_results.csv")
    if len(complete_results)==0:
        raise RuntimeError('Please check the pathname, it may not be correct.')
    else:
        complete_results=pd.read_csv(complete_results[0])
        
    params_to_add,params_best=get_params(complete_results)

    for i, param in enumerate(params_to_add):
        values[param+'_best']=params_best[i]

    check_quicksaves=glob.glob(folder+"quicksaves/*output.csv")

    if len(check_quicksaves)>0:
        
        for file in check_quicksaves:
            results_file=file.replace('output.csv','results.pkl')
            priors_file=file.replace('output.csv','priors.pkl')

            with open(priors_file, 'rb') as handle:
                priors = pickle.load(handle)
                order_of_params=priors.fitting_params
                
            with open(results_file, 'rb') as handle:
                results = pickle.load(handle)
                samples = results.samples
                weights=np.exp(results.logwt - results.logwt.max())/np.sum(np.exp(results.logwt - results.logwt.max()))
            
            params_from_priors=[]
            for o, order in enumerate(order_of_params):
                if order[0] in ['a','rp']:
                    order[0]+='/r*'
                params_from_priors.append(order[0])

            param_list=[None]*len(params_from_priors)


            results_csv=pd.read_csv(file)
            #params_to_add,params_best=get_params(results_csv)


            for p, param in enumerate(set(params_from_priors)):    
                indices=[x for x, y in enumerate(params_from_priors) if y == param]

                for i,ind in enumerate(indices):
                    selected_df=results_csv[results_csv['Parameter']==param]


                    if param_list[ind] is None:
                        if param in ['t0']:
                            param_list[ind]=param+'_-_-_'+check_index(selected_df['Epoch'].values[i])
                        else:
                            param_list[ind]=param+'_'+check_index(selected_df['Telescope'].values[i])+'_'+check_index(selected_df['Filter'].values[i])+'_'+check_index(selected_df['Epoch'].values[i])

            for p, par in enumerate(param_list):
                values=make_dict(values, par, samples[:, p])

    # Check folded mode:
    folded_files=glob.glob(folder+"*_parameters/quicksaves/*output.csv")

    if len(folded_files)>0:
        print("filter files detected")
        more_params=params_to_add-values.keys()
        
        for file in folded_files:
            results_file=file.replace('output.csv','results.pkl')
            priors_file=file.replace('output.csv','priors.pkl')

            with open(priors_file, 'rb') as handle:
                priors = pickle.load(handle)
                order_of_params=priors.fitting_params
                
            with open(results_file, 'rb') as handle:
                results = pickle.load(handle)
                samples = results.samples
                weights=np.exp(results.logwt - results.logwt.max())/np.sum(np.exp(results.logwt - results.logwt.max()))
            
            params_from_priors=[]
            for o, order in enumerate(order_of_params):
                if order[0] in ['a','rp']:
                    order[0]+='/r*'
                params_from_priors.append(order[0])

            param_list=[None]*len(params_from_priors)


            results_csv=pd.read_csv(file)
            #params_to_add,params_best=get_params(results_csv)


            for p, param in enumerate(set(params_from_priors)):    
                indices=[x for x, y in enumerate(params_from_priors) if y == param]
                selected_df=results_csv[results_csv['Parameter']==param]

                for i,ind in enumerate(indices):
                    
                    if param_list[ind] is None:
                        if param in ['t0']:
                            param_list[ind]=param+'_-_-_'+check_index(selected_df['Epoch'].values[i])
                        else:
                            param_list[ind]=param+'_'+check_index(selected_df['Telescope'].values[i])+'_'+check_index(selected_df['Filter'].values[i])+'_'+check_index(selected_df['Epoch'].values[i])

            for p, par in enumerate(param_list):
                if par in more_params:
                    #print(par)
                    values=make_dict(values, par, samples[:, p])
            
    lower_errors=[]
    upper_errors=[]
    params_to_save=[]
    telescopes=[]
    filters=[]
    epochs=[]
    issue_with_priors=[]
    for p in params_to_add:
        samples=values[p]
        best=values[p+'_best']
        try:
            le,ue=get_quantiles_on_best_val_unweighted(samples, best)
        except:
            issue_with_priors.append(p.replace('_',', '))
            le,ue=get_std_on_best_val_unweighted(samples, best)
        lower_errors.append(le)
        upper_errors.append(ue)
        _p=p.split('_')

        params_to_save.append("_".join(_p[0:-3]))
        telescopes.append(_p[-3])
        filters.append(_p[-2])
        epochs.append(_p[-1])

    if len(issue_with_priors)>0:
        print("The priors for following parameters might be too strict. Consider expanding the priors for them: \nParameter, Telescope, Filter, Epoch")
        for i in issue_with_priors:
            print(i)

    data = {'Parameter': params_to_save,'Telescope':telescopes,'Filter':filters,'Epoch':epochs, 'Best': params_best, 'Lower_error': lower_errors, 'Upper_error': upper_errors}
    df = pd.DataFrame(data)
    df.to_csv(folder+'results_with_asymmetric_errors.csv', index=False)

def get_asymmetric_errors_fitting_mode_all(folder,plot_folder):
    if folder[-1]!='/':
        folder+='/'
    if plot_folder[-1]!='/':
        plot_folder+='/'
        
    df1=pd.read_csv(plot_folder+'unfolded/batch_0_samples_results_with_asymmetric_errors.csv')
    df2=pd.read_csv(folder+'Complete_results.csv')

    try:
        params=df1['Parameter'].to_list()

        if 'rp' in params:
            id=params.index('rp')
            params[id]='rp/r*'

        if 'a' in params:
            id=params.index('a')
            params[id]='a/r*'
            
        telescopes=[None]*len(params)
        filters=[None]*len(params)
        epochs=[None]*len(params)
        param_list=[None]*len(params)

        for p, param in enumerate(set(params)):    
            indices=[x for x, y in enumerate(params) if y == param]
            selected_df=df2[df2['Parameter']==param]

            for i,ind in enumerate(indices):
                telescopes[ind]=selected_df['Telescope'].values[i]
                filters[ind]=selected_df['Filter'].values[i]
                epochs[ind]=selected_df['Epoch'].values[i]
                param_list[ind]=param

        data = {'Parameter': params,'Telescope':telescopes,'Filter':filters,'Epoch':epochs, 'Best': df1['Best'].values, 'Lower_error': df1['Lower_error'].values, 'Upper_error':df1['Upper_error'].values}
        df = pd.DataFrame(data)
        df.to_csv(folder+'results_with_asymmetric_errors.csv', index=False)
    
    except:
        df1.to_csv(folder+'results_with_asymmetric_errors.csv', index=False)