=================
Allowing for TTVs
=================

In its default mode, ``TransitFit`` assumes that there are no TTVs and fits a single :math:`t_0` and :math: `P` that are assumed applicable to all transits. In the situation where TTVs are present, ``TransitFit`` can fit :math:`t_0` to each individual epoch, allowing then for O-C investigations.

In order to do this, set ``allow_TTV=True`` in the arguments of :meth:`~transitfit._pipeline.run_retrieval`. **Note**: ``TransitFit`` cannot automatically detect if there are TTVs present in the data. You must explicitly enable this mode.

When ``allow_TTV=True``, ``TransitFit`` cannot fit for the period value, and this must be set as a fixed value in the :ref:`data input file <Data input file>`. In order to be consistent, we recommend that any investigation into TTVs in a system should have 2 stages:

1. Run ``TransitFit`` with ``allow_TTV=False``.

2. Run ``TransitFit`` with ``allow_TTV=True``, fixing the period at the best fit value from the first step.


=================
Fitting for P_dot
=================

To fit for the rate of change of the period, :math:`\dot{P}`, set ``fit_taylor_ttv=True`` in the arguments of :meth:`~transitfit._pipeline.run_retrieval`. This method also allows fitting for t0 simultaneously. 

This method is explicitly slow when the range of priors is not well constrained. We recommend that you set the priors for :math:`\dot{P}` to be as tight as possible. An initial guess for priors can be calculated using :

:meth:`~transitfit.find_ttv_priors.get_priors(input_data, P_prior, t0_prior)`.

This assumes that all the lightcurves include the transit midpoint, and there are reasonable priors for Period and t0.

=================

While fitting for P_dot, it is necessary that all the light-curves are fitted in a single batch. To ensure this transitfit selects evenly spaced out lightcurves from the input file, to ensure that the number of parameters to be fitted is within the allowed maximum parameters in a single batch.
