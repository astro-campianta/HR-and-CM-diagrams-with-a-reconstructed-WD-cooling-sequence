# HR and CM diagrams with a reconstructed white dwarf cooling sequence

This project provides the two Python scripts used to generate the stellar-evolution diagrams featured in *Nozioni di base dell'astrofisica stellare* (*Basic notions of stellar astrophysics*), the episode following the introduction of *L'atlante stellare* (*The stellar atlas*), a personal video-podcast series on stellar astrophysics developed by Camilla Pianta. 
The scripts build theoretical Hertzsprung-Russell and colour-magnitude diagrams from a PARSEC isochrone, extended with a white dwarf (WD) cooling sequence reconstructed from PARSEC progenitor lifetimes, the Kalirai (2008) initial-final mass relation, and BaSTI cooling tracks (since PARSEC alone does not include the WD phase).

---

## Pipeline

The two scripts must be run in order, since `plot_stellar_diagrams.py` reads the CSV produced by `build_wd_sequence.py`:

```
python build_wd_sequence.py
python plot_stellar_diagrams.py
```

---

## `build_wd_sequence.py`

Reconstructs a WD cooling sequence compatible with the 1 Gyr, Z = 0.0152 PARSEC isochrone used in `plot_stellar_diagrams.py`, and saves it as `wd_sequence_1gyr_solar.csv`. Three ingredients are combined: the progenitor lifetime as a function of initial mass (from PARSEC), the initial-final mass relation of Kalirai et al. (2008), and the BaSTI cooling tracks of Cassisi et al. (2007).

### Functions

- **lifetime_vs_mass**: reads a grid of PARSEC isochrones spanning several ages and, for each age, takes the highest surviving initial mass (the turn-off mass), returning the resulting lifetime-vs-initial-mass relation;
- **read_wd_track**: derives the BaSTI filename corresponding to a given WD mass and loads the corresponding cooling track;
- **interp_track**: linearly interpolates a BaSTI cooling track (temperature, luminosity, photometry) at a given cooling age, clamping below the track's earliest tabulated entry rather than extrapolating;
- **build**: combines the previous three functions over a grid of progenitor initial masses, deriving for each one a cooling age and a WD mass (via the Kalirai IFMR), interpolating between the two nearest cooling tracks, and saving the resulting sequence.

### Parameters

- **BASE_DIR**: directory containing the input files (`multiage.txt` and the `wd_models` folder) and where the output CSV is saved (default: current working directory);
- **PARSEC_COLS**: column names of the PARSEC isochrone files;
- **WD_COLS**: column names of the BaSTI cooling track files;
- **WD_MASSES**: mass grid (in M☉) of the available BaSTI cooling tracks;
- **TURNOFF_MINI**: turn-off initial mass at 1 Gyr;
- **MINI_GRID**: grid of progenitor initial masses on which the WD sequence is built, starting just above `TURNOFF_MINI`;
- **mass**: initial mass used to select a BaSTI track (`read_wd_track`);
- **track, log_cool_age**: a BaSTI cooling track and the cooling age (log₁₀ yr) at which to evaluate it (`interp_track`).

### Caveats

WD photometry (BaSTI) and the PMS-AGB sequence (PARSEC/YBC) are not computed with the same spectral library, so zero-point differences, though small, are non-zero. The Kalirai (2008) IFMR is an empirical relation calibrated on open clusters rather than derived from PARSEC physics, so combining it with the PARSEC isochrone is methodologically approximate, though common practice in the literature.

---

## `plot_stellar_diagrams.py`

Builds the HR diagram and the CM diagram from the combined PARSEC isochrone and WD sequence, with each evolutionary phase (PMS, MS, SGB, RGB, CHeB, AGB, WD) shown in its own colour. Two figures are produced: one with only the coloured evolutionary sequences, and one that additionally labels the four expected regions of the diagrams (blue giants, red giants, red dwarfs, white dwarfs).

### Functions

- **teff_to_rgb**: converts an array of effective temperatures into RGB colours by linear interpolation between a lookup table of reference colours, used to colour the WD points by their actual temperature;
- **read_parsec_isochrone**: reads a PARSEC isochrone file, sorted by initial mass, dropping the WD ignition point and the post-AGB phase;
- **read_wd_sequence**: reads the WD sequence produced by `build_wd_sequence.py` and tags it as an additional evolutionary phase;
- **combine**: concatenates the isochrone and WD sequence tables into a single table sorted by initial mass;
- **phase_segments**: identifies, for each evolutionary phase, the colour and the index range it occupies along the combined table;
- **style**: applies consistent axis styling (tick direction, spine width and colour) to diagrams;
- **draw_wd**: draws the WD sequence, colouring each point by its actual effective temperature rather than by evolutionary phase;
- **draw_theoretical**: draws the HR diagram as a continuous curve, coloured by evolutionary phase, with an open circle marking the start of each phase;
- **draw_observational**: draws the CM diagram as a subsampled sequence of discrete points per phase, mimicking single-star observational measurements;
- **annotate_quadrants**: adds the four descriptive corner labels (blue giants, red giants, red dwarfs, white dwarfs) to the diagrams;
- **plot_diagrams**: builds the full two-panel figure (HR and CM diagrams side by side) with a shared legend, optionally with the corner labels, and saves it if an output path is given.

### Parameters

- **filename**: path to the PARSEC isochrone file or to the WD sequence CSV (`read_parsec_isochrone`, `read_wd_sequence`);
- **iso**: the combined isochrone and WD sequence table;
- **wd**: the WD sequence table;
- **segments**: the list of evolutionary-phase colour/index ranges returned by `phase_segments`;
- **ax**: the matplotlib axis to draw on;
- **x, y, teff, i0, i1, marker_size**: coordinates, temperature array, index range, and marker size of the WD points to draw (`draw_wd`);
- **out**: output path for the saved figure, or `None` to skip saving;
- **show_quadrants**: whether to add the four corner labels to the diagrams (default: `False`);
- **PARSEC_COLS**, **WD_LABEL**, **TEFF_COLOR_ANCHORS**, **PHASE_GROUPS**, **AXIS_MARGIN**: module-level constants defining, respectively, the isochrone column names, the label code assigned to WD points, the temperature-to-colour lookup table, the evolutionary-phase colour scheme, and the axis margins used with and without corner labels.

---

## Installation

Requirements:

- Python 3.x
- [NumPy](https://numpy.org/)
- [pandas](https://pandas.pydata.org/)
- [Matplotlib](https://matplotlib.org/)

```
pip install numpy pandas matplotlib
```

---

## Data

Three input files are required, none of which are derived from the code itself:

- **multiage.txt**: a grid of PARSEC isochrones (CMD 3.9, v1.2S+COLIBRI physics) spanning several ages, used by `build_wd_sequence.py` to derive the lifetime-vs-initial-mass relation;
- **isochrone_parsec_1Gyr_solar_full.txt**: a single PARSEC isochrone at t = 1 Gyr, Z = 0.0152, used by `plot_stellar_diagrams.py`;
- **wd_models**: seven BaSTI white dwarf cooling tracks (Cassisi et al. 2007), with DA hydrogen atmosphere, Z = 0.020, and Johnson-Cousins photometry, one file per WD mass (0.54 to 1.10 M☉).

Both PARSEC isochrones were obtained from the PARSEC CMD web interface (stev.oapd.inaf.it/cgi-bin/cmd), and the cooling tracks from the BaSTI web interface (basti-iac.oa-abruzzo.inaf.it/wdmodels.html).

Outputs: `wd_sequence_1gyr_solar.csv` (from `build_wd_sequence.py`), and `hr_cm_diagrams.png` / `hr_cm_diagrams_regions.png` (from `plot_stellar_diagrams.py`).

