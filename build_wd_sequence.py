#!/usr/bin/env python
# coding: utf-8

# ## Building the white dwarf cooling sequence
# 
# PARSEC (CMD 3.9) does not include WD cooling sequences. Therefore, to generate a WD cooling sequence compatible with the selected $t=1$ Gyr, $Z=0.0152$ isochrone of `plot_stellar_diagrams.py` and add it to the output HR and CM diagrams, three independent ingredients are necessary:
# 
# 1. the progenitor lifetime as a function of initial mass `Mini`, for `Mini` $> 2.31$ M$_\odot$ (the turn-off mass at 1 Gyr), obtained from the same PARSEC v1.2S+COLIBRI physics used for the isochrone;
# 
# 2. initial-final mass relation (IFMR) from Kalirai et al. (2008), used to convert each initial mass into the resulting WD mass: `Mwd` $= 0.109~\cdot$ `Mini` $+~0.394$ (valid for `Mini` $\in [1,6.5]$ M$_\odot$);
# 
# 3. WD cooling tracks from Cassisi et al. (2007), downloaded from BaSTI (basti-iac.oa-abruzzo.inaf.it/wdmodels.html), and chosen to match the isochrone. They thus present the following physical features: DA hydrogen atmosphere and electron conduction opacities (standard for cooling WDs), $Z=0.020$ (the closest metallicity value to the isochrone's), and Johnson-Cousins photometry.
# 
# **Caveats**:
# - WD photometry (Johnson-Cousins, BaSTI) and the PMS-AGB sequence (Maiz-Apellaniz+Bessell, PARSEC/YBC) are not computed with the same spectral library: as a consequence, zero-point differences are small, but non-zero;
# - the Kalirai (2008) IFMR is an average relation calibrated on real open clusters. Since it is not derived from PARSEC physics, combining it with the PARSEC isochrone is methodologically inconsistent (though common practice in the literature).

# ### Imports and configuration
# 
# `BASE_DIR` is the directory containing the input files (`multiage.txt` and the `wd_models` folder with the BaSTI cooling tracks) and where the output CSV will be saved. It is set to the notebook's current working directory.

# In[ ]:


import os
import numpy as np
import pandas as pd

BASE_DIR = os.getcwd()


# ### PARSEC columns, WD masses, and the initial-mass grid
# 
# `PARSEC_COLS` names the columns of the PARSEC isochrone files, while `WD_COLS` names those of the BaSTI cooling tracks. Instead, `WD_MASSES` is the mass grid of the available cooling tracks, `TURNOFF_MINI` the turn-off mass at
# 1 Gyr, and `MINI_GRID` the grid of progenitor initial masses the WD sequence will be built on, starting just above `TURNOFF_MINI`.

# In[ ]:


PARSEC_COLS = ['Zini', 'MH', 'logAge', 'Mini', 'int_IMF', 'Mass', 'logL',
               'logTe', 'logg', 'label', 'McoreTP', 'C_O', 'period0',
               'period1', 'period2', 'period3', 'period4', 'pmode', 'Mloss',
               'tau1m', 'X', 'Y', 'Xc', 'Xn', 'Xo', 'Cexcess', 'Z',
               'mbolmag', 'Umag', 'Bmag', 'Vmag', 'Rmag', 'Imag', 'Jmag',
               'Hmag', 'Kmag']

WD_MASSES = [0.54, 0.61, 0.68, 0.77, 0.87, 1.00, 1.10]
WD_COLS = ['logt', 'Mwd', 'logTe', 'logL', 'U', 'B', 'V', 'R', 'I',
           'FUV', 'NUV', 'G2', 'BP2', 'RP2', 'G3', 'BP3', 'RP3']

TURNOFF_MINI = 2.3101215363
MINI_GRID = np.linspace(TURNOFF_MINI + 0.01, 6.5, 40)


# ### Progenitor lifetime as a function of the initial mass
# 
# Isochrones only contain stars that have not yet finished their evolution: more massive stars evolve faster, so that at any age the most massive star still present is the one that is just now completing its lifetime. This means that the total time it took for the most massive star to get there is equal to the isochrone's age.
# `lifetime_vs_mass` reads a grid of PARSEC isochrones spanning several ages and, for each one, takes this maximum mass. Sorting the resulting (mass, age) pairs by mass thus gives the lifetime as a function of the initial mass.

# In[ ]:


def lifetime_vs_mass():
    
    df = pd.read_csv(os.path.join(BASE_DIR, 'multiage.txt'), sep=r'\s+', 
                     comment='#', header=None, names=PARSEC_COLS)
    turnoff = df.groupby('logAge').Mini.max().reset_index()
    turnoff = pd.concat([turnoff, pd.DataFrame(
        {'logAge': [9.0], 'Mini': [TURNOFF_MINI]})], ignore_index=True)
    turnoff = turnoff.sort_values('Mini').reset_index(drop=True)
                                                      
    return turnoff.Mini.values, turnoff.logAge.values


# ### Reading a WD cooling track from BaSTI
# 
# Each BaSTI cooling track is stored in its own file, one per WD mass in `WD_MASSES`. The filename encodes the mass as a three-digit code equal to the mass times 100 (e.g., a WD mass of 0.54 M$_\odot$ becomes 054). `read_wd_track` derives that filename from a given mass and loads the corresponding track.

# In[ ]:


def read_wd_track(mass):
    
    tag = '%03d' % round(mass * 100)
    path = os.path.join(BASE_DIR, 'wd_models', 'COOL%sBaSTIIACZ02_DAc07opa_jcGaiaGalex' % tag)
    
    return pd.read_csv(path, sep=r'\s+', comment='#', header=None, names=WD_COLS)


# ### Interpolating a track in cooling age
# 
# BaSTI tracks tabulate temperature, luminosity, and photometry at discrete cooling ages. `interp_track` linearly interpolates between these to evaluate a track at any age in between, clamping below the track's earliest entry rather than extrapolating.

# In[ ]:


def interp_track(track, log_cool_age):
    
    log_cool_age = max(log_cool_age, track.logt.min())
    out = {}
    for col in ['logTe', 'logL', 'U', 'B', 'V', 'R', 'I']:
        out[col] = np.interp(log_cool_age, track.logt, track[col])
        
    return out


# ### Building the full WD sequence
# 
# For each initial mass on `MINI_GRID`, `build` derives a cooling age and a WD mass from the Kalirai IFMR, then interpolates between the two nearest `WD_MASSES` tracks to get the matching temperature, luminosity, and
# photometry. The results are collected into a dataframe and saved as `wd_sequence_1gyr_solar.csv`.

# In[ ]:


def build():
    
    mini_pts, logage_pts = lifetime_vs_mass()
    tracks = {m: read_wd_track(m) for m in WD_MASSES}
    rows = []
    
    for mini in MINI_GRID:
        log_lifetime = np.interp(mini, mini_pts, logage_pts)
        lifetime_yr = 10**log_lifetime
        cooling_yr = 1.0e9 - lifetime_yr
        if cooling_yr <= 0:
            continue
        mwd = 0.109 * mini + 0.394
        mwd = min(max(mwd, min(WD_MASSES)), max(WD_MASSES))
        log_cool = np.log10(cooling_yr)
        lo = max([m for m in WD_MASSES if m <= mwd], default=WD_MASSES[0])
        hi = min([m for m in WD_MASSES if m >= mwd], default=WD_MASSES[-1])
        vlo = interp_track(tracks[lo], log_cool)
        if hi == lo:
            v = vlo
        else:
            f = (mwd - lo) / (hi - lo)
            vhi = interp_track(tracks[hi], log_cool)
            v = {k: vlo[k] + f * (vhi[k] - vlo[k]) for k in vlo}
        rows.append(dict(Mini=mini, Mwd=mwd, cooling_age_yr=cooling_yr,
                          logTe=v['logTe'], logL=v['logL'],
                          Umag=v['U'], Bmag=v['B'], Vmag=v['V'],
                          Rmag=v['R'], Imag=v['I']))
        
    out = pd.DataFrame(rows).sort_values('Mini').reset_index(drop=True)
    out_path = os.path.join(BASE_DIR, 'wd_sequence_1gyr_solar.csv')
    out.to_csv(out_path, index=False)
    
    print('saved:', out_path, '(%d points)' % len(out))
    print(out.iloc[[0, len(out) // 2, -1]].to_string(index=False))
    
    return out


# ### Running the pipeline

# In[ ]:


wd_sequence = build()
wd_sequence.head()

