#!/usr/bin/env python
# coding: utf-8

# ## Hertzsprung-Russell and color-magnitude diagrams
# 
# Both diagrams are obtained from the data of the PARSEC v1.2S + COLIBRI (Bressan et al. 2012; Marigo et al. 2013; Pastorelli et al. 2019, 2020) isochrone, with age $t = 1$ Gyr and Solar metallicity $Z = 0.0152$, and a separate WD cooling sequence from PARSEC + IFMR (Kalirai et al. 2008) + cooling tracks (by Cassisi et al. 2007) from BaSTI. See `build_wd_sequence.py` for details.
# 
# The HR diagram is constructed by plotting $log(T_{eff})$ vs $log(L)$, and displayed as a solid line, representing the continuous variation of the luminosity and temperature along each evolutionary sequence. 
# Instead, the CM diagram is constructed by plotting $B-V$ vs $M_V$ (which are derived from the synthetic photometry of the isocrone, `Bmag` and `Vmag`), and displayed as a succession of discrete points to mimic single-star observational measurements.
# 
# Evolutionary sequences (PMS, MS, subgiant, RGB, CHeB, AGB, and WD) are shown in different colours.
# 
# Two figures are thus produced: the former exhibits only labelled evolutionary sequences in both the HR and the CM diagrams, whereas the latter adds the typical four regions (blue giants, red giants, red dwarfs, and white dwarfs) of the expected location of stars based on colour.

# In[ ]:


import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Latin Modern Roman', 'DejaVu Serif']
plt.rcParams['mathtext.fontset'] = 'cm'
INK = '#000000'


# ### Defining the isochrone columns

# In[1]:


PARSEC_COLS = ['Zini', 'MH', 'logAge', 'Mini', 'int_IMF', 'Mass', 
               'logL', 'logTe', 'logg', 'label', 'McoreTP', 'C_O', 
               'period0', 'period1', 'period2', 'period3', 'period4', 'pmode', 
               'Mloss', 'tau1m', 'X', 'Y', 'Xc', 'Xn', 
               'Xo', 'Cexcess', 'Z', 'mbolmag', 'Umag', 'Bmag', 
               'Vmag', 'Rmag', 'Imag', 'Jmag', 'Hmag', 'Kmag']

WD_LABEL = 20


# ### Mapping stellar temperature to display colour
# 
# `TEFF_COLOR_ANCHORS` is a lookup table of (temperature in K, hex color) pairs, giving simplified reference colours commonly used in popular astronomy, rather than colours computed from real stellar spectra. The nine anchors span the full range from red M-dwarfs to blue O-stars:
# 
# | Effective temperature (K) | Colour     | Approximate spectral type | Description         |
# |:--------------------------:|:----------:|:--------------------------:|:--------------------:|
# | 3000                        | `#ff6633`  | M                          | deep orange-red      |
# | 3700                        | `#ff9955`  | K/M                        | orange               |
# | 5000                        | `#ffd2a1`  | K                          | pale orange-yellow   |
# | 5800                        | `#fff4ea`  | G                          | warm white           |
# | 7000                        | `#fbf8ff`  | F                          | cool white           |
# | 10000                       | `#cad8ff`  | A                          | pale blue-white      |
# | 15000                       | `#a8c3ff`  | B                          | light blue           |
# | 20000                       | `#9ab6ff`  | B                          | blue                 |
# | 30000                       | `#8aabff`  | O                          | deep blue            |
# 
# `teff_to_rgb(teff_k)` turns an array of temperatures into an array of RGB colours that can be plotted directly:
# 
# 1. it splits each anchor colour into its R, G, B components;
# 2. it clips the input temperatures to the anchor range $3000-30000$ K, so that a temperature outside that range is not extrapolated, but simply takes the colour of the nearest anchor;
# 3. for each of the three channels, it linearly interpolates between the anchor temperatures to get that channel's value at the requested temperature;
# 4. it stacks the three interpolated channels back together into one (N,3)-array of RGB triples.
# 
# The result is a smooth colour gradient as a function of temperature, which is employed only for the white dwarf points (`draw_wd`) because the marker colour should reflect the star's actual temperature instead of representing its evolutionary phase, as the other colours in `PHASE_GROUPS` do.

# In[ ]:


TEFF_COLOR_ANCHORS = [
    (3000,  '#ff6633'),
    (3700,  '#ff9955'),
    (5000,  '#ffd2a1'),
    (5800,  '#fff4ea'),
    (7000,  '#fbf8ff'),
    (10000, '#cad8ff'),
    (15000, '#a8c3ff'),
    (20000, '#9ab6ff'),
    (30000, '#8aabff'),
]

def teff_to_rgb(teff_k):
    
    temps = np.array([a[0] for a in TEFF_COLOR_ANCHORS], dtype=float)
    rgbs = np.array([mcolors.to_rgb(a[1]) for a in TEFF_COLOR_ANCHORS])
    teff_k = np.clip(np.asarray(teff_k, dtype=float), temps[0], temps[-1])
    
    return np.column_stack([np.interp(teff_k, temps, rgbs[:, c]) for c in range(3)])


# ### Reading the PARSEC isochrone and the WD sequence
# 
# Three functions bring the two data sources together into one table:
# 
# 1. `read_parsec_isochrone(filename)` reads the PARSEC isochrone file and sorts it by initial mass, dropping only the single ignition point of the WD sequence (`label == 9`) and the post-AGB phase (`label == 8`);
# 
# 2. `read_wd_sequence(filename)` reads the WD sequence implemented in `build_wd_sequence.py`, and tags every row with `WD_LABEL` in order for it to be treated as an additional evolutionary phase downstream;
# 
# 3. `combine(iso, wd)` concatenates the tables produced by `read_parsec_isochrone`and `read_wd_sequence`, and re-sorts by initial mass, returning a combined table referred to as `iso`. PARSEC columns that do not exist in the WD table are left as `NaN` (this is harmless because the plotting code only reads `Mini`, `label`,`logTe`, `logL`, `Bmag`, and `Vmag`).

# In[ ]:


def read_parsec_isochrone(filename):
    
    df = pd.read_csv(filename, sep=r'\s+', comment='#', header=None, names=PARSEC_COLS)
    df = df.sort_values('Mini').reset_index(drop=True)
    df = df[~df.label.isin([8, 9])].reset_index(drop=True)
    
    return df

def read_wd_sequence(filename):
    
    df = pd.read_csv(filename)
    df['label'] = WD_LABEL
    
    return df

def combine(iso, wd):
    
    both = pd.concat([iso, wd], ignore_index=True, sort=False)
    
    return both.sort_values('Mini').reset_index(drop=True)


# ### Identifying the evolutionary phases
# 
# `PHASE_GROUPS` is a fixed, categorical colour palette assigned in evolutionary order along the isochrone, from the PMS to the WD phase. Each entry gives a short key, one or more PARSEC `label` codes associated with that phase, a description, a displaycolour, and a `source` tag distinguishing points from the continuous PARSEC physical track (`parsec`) from the separately reconstructed WD sequence (`wd`).
# 
# `phase_segments(iso)` scans the combined `iso` table and returns, for each phase present, its colour, source, and the index range where it occurs. The table is sorted by initial mass, and each evolutionary phase occupies one range of masses before the next phase begins. So all rows with a given `label` sit together in one unbroken block, and `idx.min()`/`idx.max()` mark where that block starts and ends.

# In[ ]:


PHASE_GROUPS = [
    ('PMS',  [0],         'pre-main sequence',                    '#2a78d6', 'parsec'),
    ('MS',   [1],         'main sequence',                        '#eb6834', 'parsec'),
    ('SGB',  [2],         'subgiant branch',                      '#1baf7a', 'parsec'),
    ('RGB',  [3],         'red giant branch',                     '#eda100', 'parsec'),
    ('CHeB', [4, 5],      'core helium burning',                  '#e87ba4', 'parsec'),
    ('AGB',  [6, 7],      'asymptotic giant branch (E-AGB+TP-AGB)','#008300', 'parsec'),
    ('WD',   [WD_LABEL],  'white dwarf (reconstructed sequence)', '#e34948', 'wd'),
]

def phase_segments(iso):
    
    segments = []
    for key, codes, desc, color, source in PHASE_GROUPS:
        idx = iso.index[iso.label.isin(codes)]
        if len(idx) == 0:
            continue
        segments.append((key, desc, color, idx.min(), idx.max(), source))
        
    return segments


# ### Setting the axis style for the diagrams

# In[ ]:


F_TICK, F_AX, F_TITLE = 15, 17, 17

def style(ax):
    
    ax.set_facecolor('white')
    for sp in ax.spines.values():
        sp.set_linewidth(1.1)
        sp.set_color(INK)
    ax.tick_params(which='both', direction='in', top=True, right=True,
                    colors=INK, labelsize=F_TICK, length=7, width=1.0)


# ### Drawing the WD sequence
# 
# `draw_wd` colours each point of the WD sequence by its actual effective temperature, computed from the BaSTI cooling tracks and converted to a colour via `teff_to_rgb`. A thin grey line connects the points in sequence, showing that they belong to a single ordered track even though each one is coloured individually.

# In[ ]:


def draw_wd(ax, x, y, teff, i0, i1, marker_size):
    
    ax.plot(x.iloc[i0:i1 + 1], y.iloc[i0:i1 + 1], '-', color='#b0b0b0', lw=0.9, zorder=2)
    colors = teff_to_rgb(10 ** teff.iloc[i0:i1 + 1].to_numpy())
    ax.scatter(x.iloc[i0:i1 + 1], y.iloc[i0:i1 + 1], s=marker_size,
               color=colors, edgecolor='#444444', linewidth=0.5, zorder=3)


# ### Plotting the HR diagram
# 
# `draw_theoretical` plots the isochrone as a continuous curve, since $\log(T_{eff})$ and $\log(L)$ vary smoothly with mass. Each phase is drawn as a line in its own colour, with an open circle marking its starting point.
# 
# In this way, the HR diagram is obtained.

# In[1]:


def draw_theoretical(ax, iso, segments):
    
    x, y = iso.logTe, iso.logL
    n = len(x)
    
    for i, (key, desc, color, i0, i1, source) in enumerate(segments):
        if source == 'wd':
            draw_wd(ax, x, y, iso.logTe, i0, i1, marker_size=26)
            continue
        connect = i + 1 < len(segments) and segments[i + 1][5] == source
        i1p = i1 + 1 if connect and i1 + 1 < n else i1
        ax.plot(x.iloc[i0:i1p + 1], y.iloc[i0:i1p + 1], '-', color=color, lw=2.0, solid_capstyle='round', zorder=3)
        ax.plot(x.iloc[i0], y.iloc[i0], 'o', mfc='white', mec=color, mew=1.8, ms=7, zorder=4)


# ### Plotting the CM diagram
# 
# `draw_observational` plots the same isochrone in the $B-V$ vs $M_V$ plane, using the same phase colours and `source` grouping as the HR diagram. Each phase is subsampled to roughly 18 points displayed as discrete markers, with a faint line kept in the background linking them along the isochrone.
# 
# In this way, the CM diagram is obtained.

# In[ ]:


def draw_observational(ax, iso, segments):
    
    x, y = iso.Bmag - iso.Vmag, iso.Vmag
    n = len(x)
    
    for i, (key, desc, color, i0, i1, source) in enumerate(segments):
        if source == 'wd':
            draw_wd(ax, x, y, iso.logTe, i0, i1, marker_size=20)
            continue
        connect = i + 1 < len(segments) and segments[i + 1][5] == source
        i1p = i1 + 1 if connect and i1 + 1 < n else i1
        ax.plot(x.iloc[i0:i1p + 1], y.iloc[i0:i1p + 1], '-', color=color, lw=0.8, alpha=0.55, zorder=2)
        span = i1 - i0 + 1
        step = max(1, span // 18)
        sel = list(range(i0, i1 + 1, step))
        if sel[-1] != i1:
            sel.append(i1)
        ax.scatter(x.iloc[sel], y.iloc[sel], s=16, color=color, edgecolor='white', linewidth=0.4, zorder=3)


# ### Labelling the four regions of the diagrams
# 
# `annotate_quadrants` places an italic label in each corner of both diagrams, indicating the expected positions of stars: hot stars on the left, cool stars on the right, luminous stars at the top, and faint stars at the bottom.
# 
# Therefore, the four labels are: 
# 1. hot, luminous stars (blue giants); 
# 2. cool, luminous stars (red giants); 
# 3. hot, faint stars (white dwarfs); 
# 4. cool, faint stars (red dwarfs).

# In[ ]:


def annotate_quadrants(ax):
    
    kw = dict(transform=ax.transAxes, fontsize=10.5, style='italic', color='#5a5a5a', zorder=5, linespacing=1.4,
              bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='none', alpha=0.85))
    ax.text(0.02, 0.98, 'stelle calde e luminose\n(giganti blu)', ha='left', va='top', **kw)
    ax.text(0.98, 0.98, 'stelle fredde e luminose\n(giganti rosse)', ha='right', va='top', **kw)
    ax.text(0.02, 0.02, 'stelle calde e deboli\n(nane bianche)', ha='left', va='bottom', **kw)
    ax.text(0.98, 0.02, 'stelle fredde e deboli\n(nane rosse)', ha='right', va='bottom', **kw)


# ### Building the two-panel figure
# 
# `AXIS_MARGIN` defines two sets of axis margins, `base` and `quadrants`. The `quadrants` margins are wider, leaving empty space in the four corners for `annotate_quadrants`.
# `plot_diagrams(iso, out, show_quadrants=False)` builds the full figure: it draws the HR and CM diagrams side by side, optionally adds the corner labels, and builds a shared legend below the diagrams from the entries in `PHASE_GROUPS`. If `out` is given, the figure is also saved to that path.

# In[ ]:


AXIS_MARGIN = dict(base=dict(teff=0.03, logL=0.2, bv=0.05, vmag=0.3),
                   quadrants=dict(teff=0.28, logL=1.1, bv=0.35, vmag=2.6))

def plot_diagrams(iso, out, show_quadrants=False):
    
    segments = phase_segments(iso)
    m = AXIS_MARGIN['quadrants' if show_quadrants else 'base']
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(15.5, 7.2), dpi=200)
    fig.patch.set_facecolor('white')
    fig.subplots_adjust(wspace=0.38, left=0.07, right=0.93, top=0.88, bottom=0.18)

    # HR diagram
    style(axL)
    draw_theoretical(axL, iso, segments)
    axL.set_xlim(iso.logTe.max() + m['teff'], iso.logTe.min() - m['teff'])
    axL.set_ylim(iso.logL.min() - m['logL'], iso.logL.max() + m['logL'])
    axL.set_xlabel(r'$\log\,T_{\rm eff}\ \mathrm{[K]}$', fontsize=F_AX)
    axL.set_ylabel(r'$\log\,(L/L_\odot)$', fontsize=F_AX)
    axL.set_title('Diagramma di Hertzsprung-Russell', fontsize=F_TITLE, pad=12)

    # CM diagram
    bv = iso.Bmag - iso.Vmag
    style(axR)
    draw_observational(axR, iso, segments)
    axR.set_xlim(bv.min() - m['bv'], bv.max() + m['bv'])
    axR.set_ylim(iso.Vmag.max() + m['vmag'], iso.Vmag.min() - m['vmag'])
    axR.set_xlabel(r'$B-V$', fontsize=F_AX)
    axR.set_ylabel(r'$M_V$', fontsize=F_AX)
    axR.set_title('Diagramma colore-magnitudine', fontsize=F_TITLE, pad=12)

    if show_quadrants:
        annotate_quadrants(axL)
        annotate_quadrants(axR)

    handles = []
    for key, _, color, _, _, source in segments:
        if source == 'wd':
            wd_color = teff_to_rgb(np.array([15000.0]))[0]
            handles.append(Line2D([0], [0], color='#b0b0b0', lw=1.2, marker='o',
                                   mfc=wd_color, mec='#444444', mew=0.7, ms=9, label=key))
        else:
            handles.append(Line2D([0], [0], color=color, lw=3, solid_capstyle='round', label=key))
    fig.legend(handles=handles, loc='lower center', ncol=len(handles), frameon=False,fontsize=13, 
               bbox_to_anchor=(0.5, 0.02), handlelength=1.6, columnspacing=1.8, handletextpad=0.6)

    if out is not None:
        fig.savefig(out, facecolor='white', dpi=200)
        
    return fig, axL, axR


# ### Running the pipeline

# In[ ]:


here = os.getcwd()
iso = read_parsec_isochrone(os.path.join(here, 'isochrone_parsec_1Gyr_solar_full.txt'))
wd = read_wd_sequence(os.path.join(here, 'wd_sequence_1gyr_solar.csv'))
iso = combine(iso, wd)

out_path = os.path.join(here, 'hr_cm_diagrams.png')
out_path_regions = os.path.join(here, 'hr_cm_diagrams_regions.png')

fig1, axL1, axR1 = plot_diagrams(iso, out_path, show_quadrants=False)
fig1


# In[ ]:


fig2, axL2, axR2 = plot_diagrams(iso, out_path_regions, show_quadrants=True)
fig2

