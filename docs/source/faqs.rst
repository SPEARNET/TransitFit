====
FAQs
====
What prior distributions should I use for my `rp` values?
    Obviously your use case will determine a lot of this, but, as some helpful rules of thumb, we recommend that you use a wide uniform prior for your `rp` values, rather than a Gaussian based on previous measurements. 

Why do I have to use config files? Isn't that a bit outdated?
    Fair point. We've based this around config files because when being used on hundreds of light curves at once, typing all the inputs into a list becomes very difficult to read and keep track of. In future updates, we might work on streamlining the API to allow inputs to be set directly within the code by the user.

How do I cite ``TransitFit``?
    If you find ``TransitFit`` useful in your work, please cite the `accompanying paper <https://doi.org/10.1093/mnras/stad3353>`_. If you are using BibTeX, you can use the following citation in your .bib file::

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








Why have you only included nested sampling? Why can't I use MCMC?
    This comes down partly to the personal preference of the development team, but mostly because ``TransitFit`` often has to deal with high-dimensioned fitting problems. MCMC routines often struggle in this situation, especially when the posterior space can be fairly degenerate or spiky. Nested sampling can handle these situations in a more stable way.

I've found a bug - what can I do?
    Raise the issue on the `GitHub page <https://github.com/SPEARNET/TransitFit>`_. Make sure to include as much information as possible, including any traceback messages and information on your priors etc.

Can I contribute to the project?
    Absolutely! ``TransitFit`` is an open-source project (GPL-3.0 licence) - please raise a pull request for any additions, changes, or improvements want to suggest.
