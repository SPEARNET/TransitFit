from scipy.optimize import differential_evolution, Bounds
import numpy as np
from ._utils import get_covariance_matrix


class ResultsDE:
    """
    Creating a Results class to store necessary information about the differential evolution results.
    """

    def __init__(self, sampler,prior_transform):
        """sampler.results are the results from the dynesty sampler
        prior_transform is the function to transform unit cube samples to parameter space.
        """
        results = sampler#.results

        self.logl = -np.array(results.population_energies)
        self.samples = np.array([prior_transform(i) for i in results.population])
        self.logwt = np.ones_like(self.logl)

        # Normalise weights
        self.weights = np.ones_like(self.logl)#get_normalised_weights(results)

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
                #x0=x0,
            )

        results = ResultsDE(sampler, self.prior_transform)

        return results