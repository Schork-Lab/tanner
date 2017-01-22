
from scipy.stats import norm
import pymc3 as pm
import numpy as np
from theano import shared
from scipy.stats.distributions import pareto
from scipy import optimize
import theano.tensor as T
import theano
import pandas as pd

class BayesianModel(object):
    '''
    General Bayesian Model Class for quantifying
    relationship between gene and phenotype

    Adapted from Thomas Wiecki
    https://github.com/pymc-devs/pymc3/issues/511#issuecomment-125935523

    '''

    def __init__(self, variational=True, mb=False,
                 n_chain=50000, n_trace=5000,
                 logistic=False, steps=None):
        """
        Args:
            variational (bool, optional): Use Variational Inference
            mb (bool, optional): Use minibatches
        """
        self.variational = variational
        self.cached_model = None
        self.mb = mb
        self.n_chain = n_chain
        self.n_trace = n_trace
        self.logistic = logistic
        self.steps = steps


    def cache_model(self, **inputs):
        """
        Create a cached model for the Bayesian model using
        shared theano variables for each Bayesian
        input parameter.

        Args:
            **inputs (dict): inputs for Bayesian model

        """
        self.shared_vars = self._create_shared_vars(**inputs)
        self.cached_model = self.create_model(**inputs)

    def create_model(self, **inputs):
        """
        Each instance of this class needs to define
        their PYMC3 model in here.
        """
        raise NotImplementedError('This method has to be overwritten.')

    def _create_shared_vars(self, **inputs):
        """
        For each input variable, create theano shared variable
        and set their initial values.

        Args:
            **inputs (dict): inputs for Bayesian model

        Returns:
            dict: key, value - var_name, theano.shared variable
        """
        shared_vars = {}
        for name, data in inputs.items():
            shared_vars[name] = shared(data, name=name)
        return shared_vars

    def _clean_inputs(self, inputs):
        """
        Clean the inputs, i.e. remove some
        genotype columns. Useful for some class of Bayesian models
        such as Two-Stage, where first stage involves filtering
        on certain SNPs.

        Args:
            inputs (dict): inputs for Bayesian model
        Returns:
            dict: cleaned inputs for Bayesian model
        """
        return inputs

    def run(self, **inputs):
        """
        Run cached Bayesian model using the inputs

        Args:
            **inputs (dict): inputs for Bayesian model

        Returns:
            trace: Trace of the PyMC3 inference
        """
        if self.cached_model is None:
            self.cache_model(**inputs)
        for name, data in inputs.items():
            self.shared_vars[name].set_value(data)
        if self.mb and self.variational:
            self.minibatches = zip(self._mb_generator(inputs['gwas_gen']),
                                   self._mb_generator(inputs['gwas_phen']))
        self.trace = self._inference()
        return self.trace

    def _inference(self, n_trace=None):
        """
        Perform the inference. Uses ADVI if self.variational
        is True. Also, uses minibatches is self.mb=True based
        on generators defined in self.run.

        Otherwise, uses Metropolis.

        Args:
            n_trace (int, optional): Number of steps used for trace
        Returns:
            trace: Trace of the PyMC3 inference
        """
        if n_trace is None:
            n_trace = self.n_trace

        with self.cached_model:
            if self.variational:
                if self.mb:
                    v_params = pm.variational.advi_minibatch(n=self.n_chain,
                               minibatch_tensors=self.minibatch_tensors,
                               minibatch_RVs=self.minibatch_RVs,
                               minibatches=self.minibatches,)
                else:
                    v_params = pm.variational.advi(n=self.n_chain)
                trace = pm.variational.sample_vp(v_params, draws=n_trace)
                self.v_params = v_params
            else:
                if self.steps is None:
                    self.steps = pm.Metropolis()
                start = pm.find_MAP(fmin=optimize.fmin_powell)
                trace = pm.sample(self.n_chain,
                                  step=self.steps,
                                  start=start,
                                  progressbar=True,
                                  )
                trace = trace[-n_trace:]
        self.trace = trace
        return trace

    def summary(self):
        def trace_sd(x):
             return pd.Series(np.std(x, 0), name='sd')
        def trace_quantiles(x):
            return pd.DataFrame(pm.quantiles(x, [1, 5, 25, 50, 75, 95, 99]))
        summary = pm.df_summary(self.trace, stat_funcs=[trace_sd, trace_quantiles])
        return summary




def parse_df(df, min_run):
    """
    Parse a metabolite dataframe into constituent parts.
    """
    df = df[df['run'] > min_run]  
    n_runs = len(df['run'].unique())
    run_idx = df['run'].values - (df['run'].min())

    time = df['day'].values
    time_values = np.unique(time)
    sorted_time = dict((k, i) for i, k in enumerate(time_values))
    time_idx = [sorted_time[x] for x in time]

    levels = df['metabolite'].values
    return {'run_idx': run_idx, 'time_idx': time_idx, 'time_values':time_values, 'measured_levels': levels, 'run_values': np.unique(df['run'])}

class Linear(BayesianModel):
    """
    Linear Model for Metabolite data
    """
    def __init__(self,
                 *args, **kwargs):
        """
        Args:

        """
        self.name = 'Linear'
        super(Linear, self).__init__(*args, **kwargs)

    def create_model(self, run_idx, time_idx, time_values, measured_levels, run_values):
        """
        Simple Bayesian Linear Regression

        Args:
            gwas_gen (pandas.DataFrame): GWAS genotypes
            gwas_phen (pandas.DataFrame): GWAS phenotypes

        Returns:
            pymc3.Model(): The Bayesian model
        """
        n_runs = len(np.unique(run_idx))
        n_time = len(np.unique(time_idx))
        print(n_time, n_runs, time_values)
        with pm.Model() as metabolite_model:
            intercept = pm.Normal('intercept', 0, sd=1)
            alpha = pm.Normal('alpha', mu=0, sd=1)
            sd_metabolite = pm.HalfCauchy('sd_metabolite', beta=1)
            tau = T.eye(n_time)*(1/(sd_metabolite**2))
            latent_level = pm.MvNormal('latent', mu=intercept+alpha*time_values, tau=tau, shape=n_time)
            scaling_factor = pm.HalfNormal('beta', sd=1e7, shape=n_runs)
            sd_run = pm.HalfCauchy('sd_run', beta=1, shape=n_runs)
            mu = scaling_factor[run_idx] * latent_level[time_idx]
            sd = scaling_factor[run_idx] * sd_run[run_idx]
            metabolite = pm.Normal('metabolite', mu=mu, sd=sd,
                                   observed=measured_levels)
        return metabolite_model


class Outlier(BayesianModel):
    """
    Bayesian outlier detection
    """
    def __init__(self, 
                 *args, **kwargs):
        self.name = 'Outlier'
        super(Outlier, self).__init__(*args, **kawrgs)

    def create_model(self, run_idx, time_idx, time_values, measured_levels):
        """
        Simple Bayesian Linear Regression

        Args:
            gwas_gen (pandas.DataFrame): GWAS genotypes
            gwas_phen (pandas.DataFrame): GWAS phenotypes

        Returns:
            pymc3.Model(): The Bayesian model
        """
        n_runs = len(np.unique(run_idx))
        n_time = len(np.unique(time_idx))
        print(n_time, n_runs, time_values)
        with pm.Model() as metabolite_model:
            intercept = pm.Normal('intercept', 0, sd=1)
            alpha = pm.Normal('alpha', mu=0, sd=1)
            sd_metabolite = pm.HalfCauchy('sd_metabolite', beta=1)
            mu_outlier = pm.Normal('mu_outlier', mu=0, sd=1)
            sd_outlier = pm.HalfCauchy('sd_outlier', beta=1)
            tau = t.eye(n_time)*(1/(sd_metabolite**2))
            tau_outlier = t.eye(n_time)*(1/(sd_outlier**2))
            latent_in = pm.MvNormal('latent', mu=intercept+alpha*time_values, tau=tau, shape=n_time)
            latent_out = pm.MvNormal('latent_outlier', mu=mu_outlier, tau=tau_outlier, shape=n_time)
            scaling_factor = pm.HalfNormal('beta', sd=1e7, shape=n_runs)
            sd_run = pm.HalfCauchy('sd_run', beta=1, shape=n_runs)

            frac_outliers = pm.Uniform('frac_outliers', lower=0., upper=.5)
            is_outlier = pm.Bernoulli('is_outlier', p=frac_outliers, shape=n_time) 


            metabolite_mu = scaling_factor[run_idx] * latent_in[time_idx]
            metabolite_sd = scaling_factor[run_idx] * sd_run[run_idx]
            metabolite_mu_outlier = scaling_factor[run_idx] * latent_in[time_idx]
            
            ## Extract observed y and sigma_y from dataset, encode as theano objects
            metabolite = thno.shared(np.asarray(measured_levels, dtype=theano.config.floatX), name='metabolite')

            ## Use custom likelihood using DensityDist
            likelihood = pm.DensityDist('likelihood', logp_signoise,
                                        observed={'yobs':metabolite, 'is_outlier':is_outlier,
                                          'yest_in':metabolite_mu, 'sigma_y':metabolite_sd,
                                          'yest_out':metabolite_mu_outlier})

            metabolite = pm.Normal('metabolite', mu=mu, sd=sd,
                                   observed=measured_levels)
        return metabolite_model


def logp_signoise(yobs, is_outlier, yest_in, sigma_y, yest_out):
    '''
    Equation 17 of Hogg 2010
    https://github.com/pymc-devs/pymc3/blob/master/docs/source/notebooks/GLM-robust-with-outlier-detection.ipynb


    '''   
    
    # likelihood for inliers
    pdfs_in = T.exp(-(yobs - yest_in + 1e-4)**2 / (2 * sigma_y_in**2)) 
    pdfs_in /= T.sqrt(2 * np.pi * sigma_y**2)
    logL_in = T.sum(T.log(pdfs_in) * (1 - is_outlier))

    # likelihood for outliers
    pdfs_out = T.exp(-(yobs - yest_out + 1e-4)**2 / (2 * sigma_y_in**2)) 
    pdfs_in /= T.sqrt(2 * np.pi * sigma_y**2)
    logL_out = T.sum(T.log(pdfs_out) * is_outlier)

    return logL_in + logL_out