# 2_SPHEREx_Simulate_SED

## Authors
- Yujin Yang, Woong-Seob Jeong (KASI SPHEREx Team)

This book made by Yoonsoo Bach (for HW purposes)

## Exercises
```{code-cell}
---
tags: ["hide-cell"]  # "hide-input", "hide-output"
---

%config InlineBackend.figure_format = 'retina'
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams

# We need to do it in a separate cell. See:
# https://github.com/jupyter/notebook/issues/3385
plt.style.use('default')
rcParams.update({
    'font.family': 'Times', 'font.size':12, 'mathtext.fontset':'stix',
    'axes.formatter.use_mathtext': True, 'axes.formatter.limits': (-4, 4),
    'axes.grid': True, 'grid.color': 'gray', 'grid.linewidth': 0.5,
    'xtick.top': True, 'ytick.right': True,
    'xtick.direction': 'inout', 'ytick.direction': 'inout',
    'xtick.minor.size': 4.0, 'ytick.minor.size': 4.0,  # default 2.0
    'xtick.major.size': 8.0, 'ytick.major.size': 8.0,  # default 3.5
    'xtick.minor.visible': True, 'ytick.minor.visible': True
})


import numpy as np
from astropy.table import Table
import pandas as pd
from scipy.interpolate import UnivariateSpline

def synth_phot(wlen, flux, wlen_lvf, resp_lvf, tol=1e-3):
    """
    Quick synthetic photometry routine.

    Parameters
    ----------
    wlen : `~np.ndarray`
        The wavelength of input spectrum.
    flux : `~np.ndarray`
        The flux density of input spectrum in f_nu unit
    wlen_lvf : `~np.ndarray`
        The wavelength grid of the response function
    resp_lvf : `numpy.ndarray`
        The response function at `wlen_lvf`. Assume that this is a QE.
    tol : float, optional
        Consider only wavelength range above this tolerence (peak * tol).
        The default is 1e-3.

    Returns
    -------
    Astropy.table with [wavelength, f_nu]
        wavelength is the center of the response function

    """
    idx = resp_lvf > resp_lvf.max()*tol
    mask_flux = ((wlen > wlen_lvf[idx].min())
                & (wlen < wlen_lvf[idx].max())
    wlen_all = np.concatenate([wlen[mask_flux], wlen_lvf[])

    wave_resamp = np.concatenate( (wave[index_flux], wave_lvf[index_filt]) )
    wave_resamp.sort()
    wave_resamp = np.unique(wave_resamp)
    flux_resamp = np.interp(wave_resamp, wave, flux)
    resp_resamp = np.interp(wave_resamp, wave_lvf, resp_lvf)

    return trapezoid(resp_resamp / wave_resamp * flux_resamp, wave_resamp) \
         / trapezoid(resp_resamp / wave_resamp, wave_resamp)
```

### 2.1 Simulate SPHEREx synthetic photometries of your favorite objects.
- https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs
- broadline AGNs
- strong emission line galaxies
- stars
- asteriods
- Were you able to find proper templates to simulate? If no, that's probably a good news for you!

I will select an asteroid (3200) Phaethon.

```{note}
Data used:
* Solar spectrum from [MAST HLSP SSO](https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/solsys/)
* Phaethon's reflectance spectrum from [Binzel+2019](http://smass.mit.edu/minuspubs.html) (``Data_3.0_805/a003200.visnir.txt``)
* Slope parameter $ G = 0.06 $ by [Ansdell+2014](https://ui.adsabs.harvard.edu/abs/2014ApJ...793...50A/abstract)
* Geometric albedo $ p_V $ has large uncertainty, but roughly 0.1.
* Diameter $ D $ has maximum $\sim 10\%$ uncertainty, but roughly 5 km.
```

For an asteorid at heliocentric distance $ r_h $, observer distance $ r_o $, phase angle $ \alpha $, effective size $ D $ (i.e., cross-sectional area at opposition is $ \pi D^2/4 $), geometric albedo $ p_V $ (the albedo at the opposition), reflectance spectrum $ r(\lambda) $, and phase function $ \Phi(\alpha) $, the observed spectrum is expressed as
\begin{align}
  S_o &= \frac{S_{⊙, λ}(λ)}{(r_h / \mathrm{1 au})^2} \frac{πD^2}{4} \frac{p_V r(λ)}{4π r_o^2} Φ(α)\\
  &= \frac{D^2 p_V}{16 (\mathrm{1 au})^2} S_{⊙, λ}(λ) r(λ) \frac{Φ(α)}{(r_h / \mathrm{1 au})^2 (r_o / \mathrm{1 au})^2}
   ~,
\end{align}
where the phase function is (IAU H, G model):
\begin{equation}
  Φ(α; G) = G e^{-1.87 \tan^{1.22}(α/2)} + (1-G) e^{-3.33 \tan^{0.63}(α/2)}
\end{equation}
Since $ S_{⊙, λ}(λ) $ is in W/m2/μm, I convert it to Jy by $  $