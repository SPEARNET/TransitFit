'''
A function which will run everything when given a path to light curve and
priors!
'''

from .io import read_priors_file, read_input_file, read_filter_info, parse_data_path_list, read_data_path_array, parse_priors_list, parse_filter_list
from .retriever import Retriever
from ._utils import check_batches
from .count_params import count_number_lcs

import numpy as np
import os
import pandas as pd
from pathlib import Path


def run_retrieval(data_files, priors, filter_info=None,
                  detrending_list=[['nth order', 1]],
                  limb_darkening_model='quadratic',
                  ld_fit_method='independent', fitting_mode='auto',
                  max_batch_parameters=25, batch_overlap=2,
                  host_T=None, host_logg=None, host_z=None, host_r=None,
                  nlive=300, dlogz=None, maxiter=None, maxcall=None,
                  dynesty_sample='rslice', dynesty_bounding='multi',
                  normalise=True, detrend=True,
                  results_output_folder='./output_parameters',
                  final_lightcurve_folder='./fitted_lightcurves',
                  summary_file='summary_output.csv',
                  full_output_file='full_output.csv',
                  plot_folder='./plots', plot=True,
                  marker_color='dimgrey', line_color='black', ldtk_cache=None,
                  ldtk_samples=20000, do_ld_mc=False, data_skiprows=0,
                  allow_ttv=False, filter_delimiter=None,
                  detrending_limits=None, normalise_limits=None, bin_data=False,
                  cadence=2, binned_color='red', walks=100, slices=10, 
                  n_procs=1, check_batchsizes=False, median_normalisation=False,
                  error_scaling=False, error_scaling_limits=None, 
                  ldtk_uncertainty_multiplier=1.,
                  fit_ttv_taylor=False, 
                  use_differential_evolution=False):
    '''
    Runs a full retrieval of posteriors using nested sampling on a transit
    light curve or a set of transit light curves. For more guidance on the use
    of input files and structuring, see the TransitFit documentation.

    Parameters
    ----------
    data_files : str
        The path to the data input .csv file, which contains the paths to the
        light curve data.

    priors : str
        Path to a .csv file containing prior information on each variable to be
        fitted.

    filter_info : str, optional
        Path to a .csv file containing information on the wavelengths of the
        filters that observations were made at.

        This is required if ld_fit_method is `'single'` or `'coupled'`. If not
        None and host_T, host_logg and host_z are not None, retrieval will
        include fitting realistic limb darkening parameters for the filters.
        Default is None.

    detrending_list : array_like, shape (n_detrending_models, 2)
        A list of different detrending models. Each entry should consist
        of a method and a second parameter dependent on the method.
        Accepted methods are
            ``['nth order', order]``

            ``['custom', function, [global fit indices], [filter fit indices], [epoch fit indices]]``

            ``['off', ]``
        ``function`` here is a custom detrending function. TransitFit assumes
        that the first argument to this function is times and that all
        other arguments are single-valued - TransitFit cannot fit
        list/array variables. If 'off' is used, no detrending will be
        applied to a light curve using this model.

        If a custom function is used, and some inputs to the function
        should not be fitted individually for each light curve, but should
        instead be shared either globally, within a given filter, or within
        a given epoch, the indices of where these fall within the arguments
        of the detrending function should be given as a list. If there are
        no indices to be given, then use an empty list: []
        e.g. if the detrending function is given by::

            foo(times, a, b, c, t0, P):
                # do something

        and a should be fitted globally, then the entry in the method_list
        would be ``['custom', foo, [1], [], []]``.

    limb_darkening_model : str, optional
        The limb darkening model to use. Allowed models are

            - ``'linear'``
            - ``'quadratic'``
            - ``'squareroot'``
            - ``'power2'``
            - ``'nonlinear'``

        With the exception of the non-linear model, all models are constrained
        by the method in Kipping (2013), which can be found at
        https://arxiv.org/abs/1308.0009.
        Default is 'quadratic'.

    ld_fit_method : str, optional
        Determines the mode of fitting of limb darkening parameters. The
        available modes are:
            - ``'coupled'`` : all limb darkening parameters are fitted
              independently, but are coupled to a wavelength dependent
              model based on the host parameters through `ldkt`
            - ``'single'`` : LD parameters are still tied to a model, but only
              the first filter is actively fitted. The remaining filters
              are estimated based off the ratios given by ldtk for a host
              with the given parameters. This mode is useful for a large
              number of filters, as 'coupled' or 'independent' fitting will
              lead to much higher computation times.
            - ``'independent'`` : Each LD coefficient is fitted separately for
              each filter, with no coupling to the ldtk models.
            - ``'off'`` : no limb darkening fitting will occurr. If using this
              mode, it it strongly recommended to set values for the
              Kipping q parameters using the priors file.
            -```'custom'``` : The user can provide priors for limb darkening
            -```'exoctk'``` : LDCs are calculated using exoctk model. Which is then used as priors for limb darkening. If using this mode, it is recommended to set the number of ldtk_samples parameter to be around 10, as this step is time consuming.
            -```'exotik'``` : LDCs are calculated using exotik model. Which is then used as priors for limb darkening.
        Default is `'independent'`

    fitting_mode : {``'auto'``, ``'all'``, ``'2_stage'``, ``'folded'``, ``'batched'``}, optional
        The approach TransitFit takes towards limiting the number of parameters
        being simultaneously fitted. The available modes are:
        - ``'auto'`` : Will calculate the number of parameters required to fit
          all the data simulataneously. If this is less than max_parameters,
          will set to ``'all'`` mode, else will set to ``'folded'`` if at least
          one filter has at least 3 epochs in it. Otherwise will set to
          ``'batched'``
        - ``'all'`` : Fits all parameters simultaneously, with no folding or
          batching of curves. Should be used with caution when fitting very
          large (~< 30) numbers of parameters.
        - ``'2_stage'`` : Fits in 2 stages, first detrending the light curves
          and then fitting the detrended curves simultaneously, using the
          ``'batched'`` approach if required.
        - ``'folded'`` : Useful for fitting curves with multiple epochs for each
          filter. TransitFit will fit each filter separately and produce a
          period-folded light curve for each filter, before fitting these
          simultaneously, using the ``'batched'`` approach if required.
        - ``'batched'`` : Useful for large numbers of light curves with
          relatively few shared filters, so ``'folded'`` loses large amounts of
          multi-epoch information. This mode splits the filters into sets of
          overlapping batches, runs each batch and uses the weighted means of
          each batch to produce a final result.
        Default is ``'auto'``.

    max_batch_parameters : int, optional
        The maximum number of parameters to use in a single retrieval when
        using ``'folded'`` or ``'batched'`` fitting modes. Default is 25.

    batch_overlap : in, optional
        The number of filters or epochs to overlap adjacent batches by where
        possible. This ensures that adjacent batches share information. Default
        is 2.

    host_T : tuple or None, optional
        The effective temperature of the host star, in Kelvin. Should be given
        as a (value, uncertainty) pair. Required if ld_fit_method is
        ``'single'`` or ``'coupled'``. Default is None.

    host_logg : tuple or None, optional
        The log_10 of the surface gravity of the host star, with gravity
        measured in cm/s2. Should be given as a (value, uncertainty) pair.
        Required if ld_fit_method is ``'single'`` or ``'coupled'``. Default is
        None

    host_z : tuple or None, optional
        The metalicity of the host, given as a (value, uncertainty) pair.
        Required if ld_fit_method is ``'single'`` or ``'coupled'``. Default is
        None

    host_r : tuple or None, optional
        The host radius in Solar radii, given as a (value, uncertainty) pair.
        Required for conversion of host-planet separation from AU to host radii

    nlive : int, optional
        The number of live points to use in the nested sampling retrieval.

    normalise : bool, optional
        If True, will assume that the light curves have not been normalised and
        will fit normalisation constants within the retrieval. The range to
        fit normalisation constants c_n are automatically detected using
            ``1/f_max <= c_n <= 1/f_min``
        as the default range, where ``f_min`` and ``f_max`` are the minimum and
        maximum flux values for a given light curve. Default is ``True``.

    detrend : bool, optional
        If False, no detrending will be attempted, even if specified by
        detrending list. Default is ``True``

    dlogz : float, optional
        Retrieval iteration will stop when the estimated contribution of the
        remaining prior volume to the total evidence falls below this
        threshold. Explicitly, the stopping criterion is
        ``ln(z + z_est) - ln(z) < dlogz``, where ``z`` is the current evidence
        from all saved samples and z_est is the estimated contribution from
        the remaining volume. The default is ``1e-3 * (nlive - 1) + 0.01``.

    maxiter : int or None, optional
        The maximum number of iterations to run. If None, will
        continue until stopping criterion is reached. Default is ``None``.

    maxcall : int or None, optional
        The maximum number of likelihood calls in retrieval. If None, will
        continue until stopping criterion is reached. Default is ``None``.

    dynesty_sample : str, optional
        Method used to sample uniformly within the likelihood constraint,
        conditioned on the provided bounds. Unique methods available are:
        uniform sampling within the bounds(``'unif'``), random walks with fixed
        proposals (``'rwalk'``), random walks with variable (“staggering”)
        proposals (``'rstagger'``), multivariate slice sampling along preferred
        orientations (``'slice'``), “random” slice sampling along all
        orientations (``'rslice'``), “Hamiltonian” slices along random
        trajectories (``'hslice'``), and any callable function which follows
        the pattern of the sample methods defined in dynesty.sampling.
        'auto' selects the sampling method based on the dimensionality of
        the problem (from ``ndim``). When ndim < 10, this defaults to ``'unif'``.
        When 10 <= ``ndim`` <= 20, this defaults to ``'rwalk'``. When ndim > 20,
        this defaults to ``'hslice'`` if a gradient is provided and ``'slice'``
        otherwise. ``'rstagger'`` and ``'rslice'`` are provided as alternatives for
        ``'rwalk'`` and ``'slice'``, respectively. Default is ``'rslice'``.

    dynesty_bounding : {``'none'``, ``'single'``, ``'multi'``, ``'balls'``, ``'cubes'``}, optional
        The decomposition to use in sampling. Default is 'multi'

    results_output_folder : str, optional
        Folder to save results to. TransitFit will create subfolders within
        this if folded or batched runs are being used. Default is
        ``'./output_parameters'``

    fitted_lightcurve_folder : str, optional
        The folder to save fitted light curves to. These files contain the
        normalised and detrended light curves, as well as the best fit curve.
        Default is ``'./fitted_lightcurves'``

    summary_file : str, optional
        The file to save the summarised final parameter results to. These are
        calculated by taking weighted averages over any batched fitting.
        Default is ``'summary_output.csv'``

    full_output_file : str, optional
        The file to save the full, non-summarised results to. This file gives
        the results for each batch, without averaging over batches to get
        summaried results. Default is ``'full_output.csv'``

    plot_folder : str, optional
        Path to folder to save plots to. Default is ``'./plots'``

    plot : bool, optional
        If ``True``, will plot all fitted light curves within the fitting routine,
        including any from partial fitting (eg, single filter modes). Default
        is ``True``.

    marker_color : matplotlib color, optional
        The colour to plot data points on plots. Default is ``'dimgray'``.

    line_colour : matplotlib color, optional
        The colour to plot best fit light curves on plots. Default is ``'black'``.

    ldtk_cache : str, optional
        This is the path to cache LDTK files to. If not specified, will
        default to the LDTK default.

    ldtk_samples : int, optional
        Controls the number of samples taken by PyLDTk when calculating LDCs
        when using ``'coupled'`` or ``'single'`` modes for limb darkening
        fitting. Default is ``20000``. When using 'exoctk' mode, set this to 10. 

    do_ld_mc : bool, optional
        If ``True``, will use MCMC sampling to more accurately estimate the
        uncertainty on intial limb darkening parameters provided by PyLDTk.
        Default is ``False``.

    data_skiprows : int, optional
        The number of rows to skip when reading in light curve data from a .txt
        file. Default is ``0``.

    allow_ttv : bool, optional
        If ``True``, will fit t0 for each epoch individually. Default is
        ``False``.

    filter_delimiter : str, optional
        The delimiter in filter profile files. Default is ``None``, which will
        lead to ``pandas`` trying to auto detect the delimiter.

    detrending_limits : list, optional
        The bounds on detrending coefficients, given as (lower, upper) pair for
        each detrending method. IF not provided, will default to ±10

    bin_data : bool, optional
        If True, any folded light curves will be plotted with data binned to an
        observing cadence given by `cadence`. Default is False.

    cadence : float, optional
        The observing cadence, in minutes, to bin data to if `bin_data` is
        True. Default is 2 (mirroring TESS observations)

    binned_color : str, optional
        The color to use for binned data. Default is `'red'`.

    n_procs : int, optional
        The maximum number of processes to use when running batches. If >1,
        will run batches in parallel. Default is 1.
    
    median_normalisation : bool, optional
        Normalising lightcurves by the median of the flux value reduces the runtime in many cases.
        Default is False.

    error_scaling : bool, optional
            If True, scales the errorbars in the lightcurves following 
            https://emcee.readthedocs.io/en/stable/tutorials/line/
    
    error_scaling_limits: list, optional
            If error_scaling=True, this is the limit of the parameter.
    
    ldtk_uncertainty_multiplier: float, optional
        (From LDTK:) The uncertainty multiplier ϵ is a subjective factor
        that defines how strongly the LD profile (or the prior created
        from it) constrains the final analysis (that is, how much we
        trust the stellar atmosphere models used to create the profiles.)
    
    fit_ttv_taylor: bool, optional
        If True, will fit the TTVs using a Taylor expansion of the ephemeris
        equation. Default is False.

    use_differential_evolution: bool, optional
        If True, will use the differential evolution sampler instead of dynesty.
        Default is False.

    Returns
    -------
    results : dict
        The results returned by ``Retriever.run_dynesty()``
    '''
    print('Starting transitfit')

    # If required to check for uneven batch sizes
    if check_batchsizes:
        data_files = check_batches(allow_ttv, data_files)
    # Load in the data and work out number of telescopes, filters, and epochs
    lightcurves, detrending_index_array = read_input_file(data_files)

    n_telescopes = lightcurves.shape[0]
    n_filters = lightcurves.shape[1]
    n_epochs = lightcurves.shape[2]

    if fit_ttv_taylor:
        df_priors=pd.read_csv(priors)

        other_params=df_priors.loc[df_priors["Distribution"] != "fixed"]["Parameter"].to_list()
        if "P" in other_params:
            other_params.remove("P")
        
        number_lcs=count_number_lcs(
            inputdata=data_files,
            other_params=other_params,
            limb_darkening_model=limb_darkening_model,
            ld_fit_method=ld_fit_method,
            normalise=normalise,
            detrend=detrend,
            error_scaling=error_scaling,
            max_batch_parameters=max_batch_parameters,
            detrending_list=detrending_list,
        )
        df_inputs = pd.read_csv(data_files)
        if number_lcs<len(df_inputs):
            # Check for both "Epoch" and "Epochs" columns
            epoch_col = "Epochs" if "Epochs" in df_inputs.columns else "Epoch"
            df_inputs = df_inputs.sort_values(epoch_col, ignore_index=True)
            indices = np.unique(np.linspace(0, len(df_inputs) - 1, number_lcs, dtype=int))
            df_inputs = df_inputs.iloc[indices]
            Path(results_output_folder).mkdir(parents=True, exist_ok=True)
            
            df_inputs[epoch_col]=np.arange(0,len(df_inputs))
            
            # Fix filter mapping in priors file for the subset
            # Get the original and new filter mappings
            original_filters = df_inputs['Filter'].values
            unique_original_filters = np.unique(original_filters)
            filter_mapping = {old_filter_idx: new_idx for new_idx, old_filter_idx in enumerate(unique_original_filters)}

            df_inputs['Filter'] = df_inputs['Filter'].map(filter_mapping)
            df_inputs = df_inputs.reset_index(drop=True)
            df_inputs.to_csv(results_output_folder+'/lightcurves_for_ttv.csv',index=False)
            
            # Update priors file to match the new filter indices
            df_priors_updated = df_priors.copy()
            filter_dependent_params = ['rp', 'q0', 'q1', 'q2', 'q3', 'u0', 'u1', 'u2', 'u3']
            
            # Create new priors file with updated filter indices
            new_priors_rows = []
            for _, row in df_priors.iterrows():
                param = row['Parameter']
                if param in filter_dependent_params and not pd.isna(row.get('Filter', np.nan)):
                    old_filter_idx = int(row['Filter'])
                    if old_filter_idx in filter_mapping:
                        # This filter is still being used, update the index
                        new_row = row.copy()
                        new_row['Filter'] = filter_mapping[old_filter_idx]
                        new_priors_rows.append(new_row)
                    # If old_filter_idx not in filter_mapping, skip this row (filter not used)
                else:
                    # Non-filter dependent parameter, keep as is
                    new_priors_rows.append(row)
            
            # Save updated priors file
            df_priors_updated = pd.DataFrame(new_priors_rows)
            updated_priors_path = Path(results_output_folder) / 'priors_for_ttv.csv'
            df_priors_updated.to_csv(updated_priors_path, index=False)
            
            # Use the updated priors file
            priors = str(updated_priors_path)
        
            # Create a new filter_info file with the updated filters
            if filter_info is not None:
                df_filter_info = pd.read_csv(filter_info)
                # Select the filters that are actually used in the subset
                df_filter_info = df_filter_info[df_filter_info['Filter'].isin(unique_original_filters)]
                df_filter_info = df_filter_info.reset_index(drop=True)
                # Update the filter indices to match the new mapping
                df_filter_info['Filter'] = df_filter_info['Filter'].map(filter_mapping)
                # Save the updated filter_info file
                updated_filter_info_path = Path(results_output_folder) / 'filter_info_for_ttv.csv'
                df_filter_info.to_csv(updated_filter_info_path, index=False)
                filter_info = str(updated_filter_info_path)


            
        data_files=df_inputs
        info = data_files.values

        data_path_array, detrending_index_array = parse_data_path_list(info)

        lightcurves = read_data_path_array(data_path_array, skiprows=0)

        lightcurves, detrending_index_array = lightcurves, detrending_index_array

        # Load in the data and work out number of telescopes, filters, and epochs

        n_telescopes = lightcurves.shape[0]
        n_filters = lightcurves.shape[1]
        n_epochs = lightcurves.shape[2]

        fitting_mode = "all"
        allow_ttv=False # This mode fits t0 for each epoch individually.

    # Set up the Retriever
    retriever = Retriever(data_files, priors, n_telescopes, n_filters, n_epochs,
                          filter_info, detrending_list, limb_darkening_model,
                          host_T, host_logg, host_z, host_r, ldtk_cache,
                          ldtk_samples, do_ld_mc, data_skiprows, allow_ttv,
                          filter_delimiter, detrending_limits, normalise, 
                          normalise_limits,detrend,median_normalisation,
                          error_scaling, error_scaling_limits, 
                          ldtk_uncertainty_multiplier,ld_fit_method, fit_ttv_taylor, use_differential_evolution)

    # This part has been handled by the Retriever
    if ld_fit_method in ['exoctk','exotik']:
        ld_fit_method='custom'

    # Run the retrieval!
    results = retriever.run_retrieval(ld_fit_method, fitting_mode,
                                      max_batch_parameters, maxiter, maxcall,
                                      dynesty_sample, nlive, dlogz, plot,
                                      results_output_folder,
                                      final_lightcurve_folder, summary_file,
                                      full_output_file, plot_folder,
                                      marker_color, line_color, dynesty_bounding, normalise,
                                      detrend, batch_overlap, bin_data, cadence, binned_color,
                                      walks, slices, n_procs)

    return results
