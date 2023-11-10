Welcome to TransitFit!
======================

``TransitFit`` is a Python 3.X package for fitting multi-telescope, multi-filter, and multi-epoch exoplanetary transit observations. It uses the `batman <https://www.cfa.harvard.edu/~lkreidberg/batman/index.html>`_ transit model and nested sampling routines from `dynesty <https://dynesty.readthedocs.io/en/latest/index.html>`_.

``TransitFit`` can include in its likelihood calculations the effect of host star characteristics and observation filter profiles on limb-darkening coefficients (LDCs), which we refer to as :ref:`'coupling' the LDCs <Limb-darkening>`. It can also perform per-telescope detrending simultaneously with the fitting of other parameters.

See :ref:`Getting Started` for instructions on how to use ``TransitFit``. You can also find more information in the `TransitFit paper <https://ui.adsabs.harvard.edu/abs/2021arXiv210312139H/abstract>`_


Installation
============
``TransitFit`` is compatible with Python 3.X installations. To run, ``TransitFit`` requires the following packages:

* numpy
* scipy
* pandas
* matplotlib
* corner
* `batman <https://www.cfa.harvard.edu/~lkreidberg/batman/index.html>`_
* `dynesty <https://dynesty.readthedocs.io/en/latest/index.html>`_
* `ldtk <https://github.com/hpparvi/ldtk>`_

To install the most recent stable version, run::

    pip install transitfit

Alternatively, you can install it direct from the source by downloading the project from the `GitHub page <https://github.com/joshjchayes/TransitFit>`_ and running::

    pip install -e .


Citing
======
If you find ``TransitFit`` useful in your work, please cite the `accompanying paper <https://ui.adsabs.harvard.edu/abs/2021arXiv210312139H>`_. If you are using BibTeX, you can use the following citation in your .bib file::

    @article{10.1093/mnras/stad3353,
    author = {Hayes, J J C and Priyadarshi, A and Kerins, E and Awiphan, S and McDonald, I and A-thano, N and Morgan, J S and Humpage, A and Charles, S and Wright, M and Joshi, Y C and Jiang, Ing-Guey and Inyanya, T and Padjaroen, T and Munsaket, P and Chuanraksasat, P and Komonjinda, S and Kittara, P and Dhillon, V S and Marsh, T R and Reichart, D E and Poshyachinda, S and The SPEARNET Collaboration},
    title = "{TransitFit: combined multi-instrument exoplanet transit fitting for JWST, HST and ground-based transmission spectroscopy studies}",
    journal = {Monthly Notices of the Royal Astronomical Society},
    pages = {stad3353},
    year = {2023},
    month = {11},
    abstract = "{We present TransitFit†, a package designed to fit exoplanetary transit light-curves. TransitFit offers multi-epoch, multi-wavelength fitting of multi-telescope transit data. TransitFit allows per-telescope detrending to be performed simultaneously with transit parameter fitting, including custom detrending. Host limb darkening can be fitted using prior conditioning from stellar atmosphere models. We demonstrate TransitFit in a number of contexts. We model multi-telescope broadband optical data from the ground-based SPEARNET survey of the low-density hot-Neptune WASP-127 b and compare results to a previously published higher spectral resolution GTC/OSIRIS transmission spectrum. Using TransitFit, we fit 26 transit epochs by TESS to recover improved ephemeris of the hot-Jupiter WASP-91 b and a transit depth determined to a precision of 111 ppm. We use TransitFit to conduct an investigation into the contested presence of TTV signatures in WASP-126 b using 180 transits observed by TESS, concluding that there is no statistically significant evidence for such signatures from observations spanning 27 TESS sectors. We fit HST observations of WASP-43 b, demonstrating how TransitFit can use custom detrending algorithms to remove complex baseline systematics. Lastly, we present a transmission spectrum of the atmosphere of WASP-96 b constructed from simultaneous fitting of JWST NIRISS Early Release Observations and archive HST WFC3 transit data. The transmission spectrum shows generally good correspondence between spectral features present in both datasets, despite very different detrending requirements.}",
    issn = {0035-8711},
    doi = {10.1093/mnras/stad3353},
    url = {https://doi.org/10.1093/mnras/stad3353},
    eprint = {https://academic.oup.com/mnras/advance-article-pdf/doi/10.1093/mnras/stad3353/52799695/stad3353.pdf},
}


.. toctree::
    :hidden:
    :maxdepth: 2

    quickstart
    configfiles
    limb_darkening
    detrending
    manyparams
    ttvs
    faqs
    api
