
# TransitFit

[TransitFit](https://transitfit.readthedocs.io/en/latest/) has been developed as part of the ongoing effort of the Spectroscopy and Photometry of Exoplanet Atmospheres Research Network (SPEARNET). [SPEARNET](https://doi.org/10.1093/mnras/stz783) is undertaking a survey of exoplanet atmospheres using transmission spectroscopy.

**Original Author** - Joshua Hayes (University of Manchester)

**Current maintainer** - Akshay Priyadarshi (University of Manchester) [email: akshay.priyadarshi@manchester.ac.uk]

## Overview
TransitFit is designed for exoplanetary transmission spectroscopy studies and offers a flexible approach to fitting single or multiple transits of an exoplanet at different observation wavelengths.  It possesses the functionality to efficiently couple host limb-darkening parameters to a range of physical models across different wavelengths, through the use of the [Limb darkening toolkit (ldtk)](https://github.com/hpparvi/ldtk) and the [Kipping parameterisations of two-parameter limb darkening models](https://arxiv.org/abs/1308.0009). TransitFit uses [batman](https://www.cfa.harvard.edu/~lkreidberg/batman/index.html) to handle transit light curve modelling, and sampling and retrieval uses the nested sampling algorithm available through [dynesty](https://dynesty.readthedocs.io/en/latest/index.html).

<a name="installation"></a>
## Installation

Please note that TransitFit currently only runs on UNIX-based machines.

Along with an installation of Python 3 (with the standard Conda distribution packages), TransitFit requires the following packages to be installed:

- [dynesty](https://dynesty.readthedocs.io/en/latest/index.html)

- [batman](https://www.cfa.harvard.edu/~lkreidberg/batman/index.html)

- [Limb darkening toolkit (ldtk)](https://github.com/hpparvi/ldtk)


<a name="guide"></a>
## User guide
The documentation for TransitFit can be found [here](https://transitfit.readthedocs.io/en/latest/)

<a name="data"></a>
## Observational Data
The observational data of WASP-127b obtained by SPEARNET and the TRANSITFIT transmission spectrum of WASP-96b can be downloaded from [here](https://cdsarc.cds.unistra.fr/viz-bin/cat/J/MNRAS/527/4936) and [here](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/MNRAS/527/4936).

<a name="citing"></a>
## Citing TransitFit
If you have used TransitFit in your work, please cite the [accompanying paper](https://doi.org/10.1093/mnras/stad3353). If you are using BibTeX, you can add the following to your .bib file:

```
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
```
