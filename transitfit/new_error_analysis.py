import numpy as np
import glob
import pandas as pd
import pickle



def make_dict(val_dict, key, vals):
    """makes dictionary using the given values. if the key already 
    exists, then the values get appended. 

    Args:
        val_dict (dict): the dictionary where the values are to be added
        key (str): the key of the value
        vals (array): the values to be added
    """
    if key in val_dict:
        val_dict[key] = np.ndarray.flatten(np.append(val_dict[key], vals))
    else:
        val_dict[key] = np.ndarray.flatten(vals)
    
    return val_dict

def check_index(index):
    """Checks if index is number or None. Returns '-' if None."""
    if index=='-':
        return index
    if index.isnumeric():
        return str(index)
    if index is None or np.isnan(index):
        return '-'
    return str(index)

def get_params(results_csv):
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

    return -np.abs(np.percentile(samples[samples<best_val], 31.73)-best_val), np.abs(np.percentile(samples[samples>best_val], 68.27)-best_val)

def HST_detrending():
    pass

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

                for i,ind in enumerate(indices):
                    selected_df=results_csv[results_csv['Parameter']==param]


                    if param_list[ind] is None:
                        param_list[ind]=param+'_'+check_index(selected_df['Telescope'].values[i])+'_'+check_index(selected_df['Filter'].values[i])+'_'+check_index(selected_df['Epoch'].values[i])

            for p, par in enumerate(param_list):
                if par in more_params:
                    print(par)
                    values=make_dict(values, par, samples[:, p])
            
    lower_errors=[]
    upper_errors=[]
    params_to_save=[]
    for p in params_to_add:
        samples=values[p]
        best=values[p+'_best']
        le,ue=get_quantiles_on_best_val_unweighted(samples, best)
        lower_errors.append(le)
        upper_errors.append(ue)
        params_to_save.append(p.replace('_',','))

    data = {'Parameter,Telescope,Filter,Epoch': params_to_save, 'Best': params_best, 'Lower_error': lower_errors, 'Upper_error': upper_errors}
    df = pd.DataFrame(data)
    df.to_csv(folder+'results_with_asymmetric_errors0.csv', index=False)