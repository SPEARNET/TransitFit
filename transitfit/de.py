from scipy.optimize import differential_evolution, Bounds
import numpy as np
from ._utils import get_covariance_matrix, get_normalised_weights
import psutil

class ResultsDE:
    """
    Creating a Results class to store necessary information about the differential evolution results.
    """

    def __init__(self, sampler,prior_transform, callback_instance=None):
        """sampler.results are the results from the dynesty sampler
        prior_transform is the function to transform unit cube samples to parameter space.
        """
        #results = sampler#.results
        #breakpoint()
        samples, unique_idx = np.unique(np.concatenate(callback_instance.population_history),axis=0, return_index=True)
        
        
        self.samples = np.array([prior_transform(i) for i in samples])
        self.logl = - np.concatenate(callback_instance.energy_history)[unique_idx]# np.ones_like(samples[:,0]) # Placeholder since DE does not provide log-likelihoods
        self.logwt = np.ones_like(self.logl)

        # Normalise weights
        unnormalized_weights = np.exp(self.logl-np.max(self.logl))
        self.weights = unnormalized_weights / np.sum(unnormalized_weights)
        #self.weights = np.ones_like(self.logl)#get_normalised_weights(results)

        # Calculate covariance matrix and use to get uncertainties
        cov = get_covariance_matrix(self)
        diagonal = np.diag(cov)
        uncertainties = np.sqrt(diagonal)

        self.cov = cov
        self.uncertainties = uncertainties

        # Get the 16th and 84th percentiles to use as upper and lower errors
        # This is arguably better than using the covariances(???)
        median = np.median(self.samples, axis=0)
        per_16 = np.percentile(self.samples, 16, axis=0)
        per_84 = np.percentile(self.samples, 84, axis=0)

        self.median = median
        self.lower_err = abs(median - per_16)
        self.upper_err = abs(per_84 - median)
        self.per_16 = per_16
        self.per_84 = per_84

        # Save the best fit results for easy access
        self.best = np.array(prior_transform(sampler.x))#[np.argmax(self.logl)]

class CallbackDE:
    def __init__(self,):
        self.population_history = []
        self.energy_history = []

    def __call__(self, intermediate_result):
        self.population_history.append(intermediate_result['population'])
        self.energy_history.append(intermediate_result['population_energies'])

        print("Stored generations:",len(self.population_history))
        memory_info = psutil.virtual_memory()
        available_gb = memory_info.available / (1024**3)  # Convert to GB
        used_gb = memory_info.used / (1024**3)
        print(f"Available RAM: {available_gb:.2f} GB")

        if available_gb/(available_gb + used_gb) < .1:
            print(f"WARNING: Low memory! Available: {available_gb:.2f} GB, Used: {used_gb:.2f} GB")
            raise MemoryError("Insufficient memory to continue optimization.")
        
        #return False

class DifferentialEvolutionSampler:
    def __init__(self, prior_transform, log_likelihood, ndim=None):
        """
        Initialize the Differential Evolution Sampler.

        prior_transform: Function to transform unit cube samples to parameter space.
        log_likelihood: Function to compute the log likelihood of a sample.
        ndim: Number of dimensions of the parameter space (optional, inferred from prior_transform).

        """
        self.prior_transform = prior_transform
        self.log_likelihood = log_likelihood
        self.results = None
        self.sampler = None
        self.ndim = ndim

    def neg_single_run(self,i):
        """ This function takes a sample from the unit cube, transforms it, and returns the negative log likelihood."""
        unit_cube = i
        cube = self.prior_transform(unit_cube)
        likelihood = self.log_likelihood(cube)
        return -likelihood

    def run(self,nlive=1000, maxiter=None, workers=1):
        """
        Run the differential evolution sampler.
        nlive: Number of live points (population size).
        maxiter: Maximum number of iterations (optional).
        workers: Number of workers to use for multiprocessing (default is 1 to avoid pickling issues).
        """
        bounds = Bounds(np.zeros(self.ndim), np.ones(self.ndim))

        callback_instance = CallbackDE()
        sampler = differential_evolution(
                self.neg_single_run,
                bounds,
                disp=True,
                strategy='randtobest1bin',
                popsize=nlive,
                workers=workers,
                updating='deferred',
                polish=False,
                init='sobol',
                maxiter=maxiter,
                callback=callback_instance,
                #x0=x0,
            )

        results = ResultsDE(sampler, self.prior_transform, callback_instance)

        return results
