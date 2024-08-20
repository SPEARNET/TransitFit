import pandas as pd
import numpy as np
from transitfit.ttv_fitting import taylor_series

def test_priors(p_prime, p_dprime,initial_guess_epochs, times_first, times_last,P_prior):
    period_all=np.empty(0)
    t0_all = np.empty(0)

    t_start=0
    for i in range(initial_guess_epochs[-1]+1):
        period_all=np.append(period_all,taylor_series(P_prior, p_prime, p_dprime, t_start))
        t0_all=np.append(t0_all,t_start)
        t_start+=period_all[-1]

    return np.logical_and(times_first<t0_all[initial_guess_epochs],times_last>t0_all[initial_guess_epochs]).all()


def get_priors(input_data, P_prior, t0_prior):

    data = pd.read_csv(input_data)

    # We prepare an initial guess for the epochs of the transits
    times=[]
    times_last=np.empty(0)
    times_first=np.empty(0)

    for d in data['Path']:
        df=pd.read_csv(d)
        times.append(df['Time'].to_numpy(dtype=np.float64))
        times_last=np.append(times_last,df['Time'].to_numpy(dtype=np.float64)[-1])
        times_first=np.append(times_first,df['Time'].to_numpy(dtype=np.float64)[0])
        
    times_last-=t0_prior   
    times_first-=t0_prior
    initial_guess_epochs=np.array(times_last//P_prior,dtype=int)


    test_priors_pprime=np.geomspace(1e-1,1e-15,100)
    p_dprime=0

    for i,p_prime in enumerate(test_priors_pprime):
        if test_priors(p_prime, p_dprime,initial_guess_epochs, times_first, times_last,P_prior):
            p_prime_max = test_priors_pprime[max(i-1,0)]
            break

    for i,p_prime in enumerate(test_priors_pprime):
        if test_priors(-p_prime, p_dprime,initial_guess_epochs, times_first, times_last,P_prior):
            p_prime_min = -test_priors_pprime[max(i-1,0)]
            break

    p_prime_limit = max(abs(p_prime_max),abs(p_prime_min))

    test_priors_pdprime=np.geomspace(p_prime_limit,1e-15,100)
    p_prime=0

    for i,p_dprime in enumerate(test_priors_pdprime):
        if test_priors(p_prime, p_dprime,initial_guess_epochs, times_first, times_last,P_prior):
            p_dprime_max = test_priors_pdprime[max(i-1,0)]
            break

    for i,p_dprime in enumerate(test_priors_pdprime):
        if test_priors(p_prime, -p_dprime,initial_guess_epochs, times_first, times_last,P_prior):
            p_dprime_min = -test_priors_pdprime[max(i-1,0)]
            break

    p_dprime_limit = max(abs(p_dprime_max),abs(p_dprime_min))
    print(f"priors for p_prime and p_dprime are: {p_prime_limit},{p_dprime_limit}")