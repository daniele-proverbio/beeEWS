# CSD-based EWS for honeybee colony fragility
Code and data to reproduce the results in paper on EWS and honeybee colony dynamics.

## Codes
### Model simulations and computational EWS
`WC_model_clean.ipynb`:
- Simulates the model from Khoury et al (2011) described in the paper, in its stochastic version using _sdeint_.
- Save its results in a database for subsequent analysis.
- Extracts and estimates EWS related to AC(1), Variance and Shannon Entropy, for different sampling resolutions.
- Generates corresponding figures.

The code is in Python, tested on v. 3.13. It requires the _scipy_ and _sdeint_ libraries. 

### Data analysis
In folder **MATLABL_src**.  
The code is in Matlab, tested on v. 2025b. The functions require the _Symbolic Toolbox_ and the _Signal Processing Toolbox_. 

`analysisTemperature_low_freq.m`: Analyse temperature data for the low-frequency dataset.  
- 1: Loads the respective datasets about Infested and Uninfested colonies, as well as External temperature (provided as a single xlsx files), as well as colony strength (given as manual input inside the code
- 2: Visualizes the data (Temperature and colony strength, over time)
- 3a: a failed test of naive detrending by external temperature (negative result to check that such detrending is not meaningful.
- 3b: Change point analysis and Granger causality test, to statistically observe if T-ext influences directly T_int
- 4: Look for law between T_ext and T_int, using a polynomial fit (calls _FitTemp_ function).

`analysisTemperature_high_freq.m`: For points 1 and 2 as in _analysisTemperature_low_freq.m_, but for the high-frequency dataset.

`EWS_extract_low_freq`: Extract and analyse statistical indicators and EWS for the low-frequency dataset, for infested and uninfested hives.
- 1: Load data
- 2: Calculate Var, AC(1) and CV after LOESS-detrending on an expanding window
- 3: plot (incl. chage points detected in _analysisTemperature_low_freq.m_)
- 4: Estimate true and false positive signals, for the ROC, by testing various thresholds of kendall's $\tau$ from the trends of the considered indicators
- 5: ROC plot, AUC calculation, estimation of best threshold for $\tau$

`EWS_extract_high_freq`: Extract and analyse statistical indicators and EWS for the high-frequency dataset, for infested and uninfested hives
- 1: Load data. Also remove the daily oscillations via differencing
- 2: Calculate Var, AC(1), CV and spectral slope (here, uses the _estimate_noise_color.m_ function) after LOESS-detrending on an expanding window
- 3: Plot indicators, on the same timeframe as for the low-frequency dataset
- 4: Estimate true and false positive signals, for the ROC, by testing various thresholds of kendall's $\tau$ from the trends of the considered indicators
- 5: ROC plot, AUC calculation, estimation of best threshold for $\tau$
- 6: noise properties (distribution and frequency properties) before and after the changepoint, after Savitzky-Golay global smoothing. Uses  _estimate_noise_color_ and _plot_psd_fit_ functions.

`phase_plane_plotter_khoury.m`: Estimates nullclines, vector field, equilibria and their stability, for the model by Khoury et al (2011) analysed in the paper. Plots the results, for different values of $m$; change index of mu_values at line 47 to get the tested values (0.2, 0.4 or $m_c$).

`other_figures.m`: Creates other companion plots for the article (examples of noisy dynamics and bifurcation plots of the colony models).

### Utils
Companion functions, used by the main codes:
- _estimate_noise_color_: obtains statistics about the detrended temperature noise
- _fitTemp_: best fit of T_int by T_ext
- _plot_psd_fit_: plots power spectral density from noise color fit

## Data
The folder contains published (shared with permission) and unpublished data about hive **temperature** (internal and external), sampled at different time intervals, for infested and uninfested hives.
- _Dataset_internal_infested_ and _Dataset_internal_uninfested_: one-day averaged internal temperature data, from the publication "Frizzera, D., et al. _Varroa destructor exacerbates the negative effect of cold contributing to honey bee mortality_, J. Insect Physiol., 151, 104571 (2023)". Cite the original publication if you reuse this dataset.
- _Dataset_external_: average daily external temperature data were derived from the regional meteorological observatory (ARPA FVG – OSMER and GRN, https://www.meteo.fvg.it/)
- _all_temperature_highfreq_: internal temperatures at 20-minutes sampling rate, for each hive.
Data about **colony strength** are already embedded into the MATLAB files, along with their dates.

## Reproducibility
### Requirements
Python library versions:
- scipy==1.17.0
- sdeint==0.3.0
- matplotlib==3.10.9

MATLAB: v. 2025b as stable; requires _Symbolic Toolbox_ and the _Signal Processing Toolbox_.

### Installation
The code files (ipynb and **Data analysis** files) are standalone. Download and open them one by one, with your preferred editor. Keep the folder structure to link the code directly to the corresponding data.

### Output
Once run, the code files generate all figures and results included in the paper. The output is typically retrned within a few seconds or minutes on a standard laptop.

## Credits
- **Code** developed by Daniele Proverbio in the context of the ERC INSPIRE grant (project number 101076926). If you reuse the code in full or in part, please cite the original manuscript: "Early waring signals based on critical slowing down have little practical applicability to anticipate honeybee health collapse". Code released under Apache 2.0 licence, see LICENCE for details.
- **Data** collected by Davide Frizzera and Desiderato Annoscia in the context of the Horizon 2020 PoshBee project (no. 773921) and Italian PRIN 2017 - UNICO (2017954WNT). If you reuse the datasets, please cite the first manuscript they appeared in ("Frizzera, D., et al. _Varroa destructor exacerbates the negative effect of cold contributing to honey bee mortality_, J. Insect Physiol., 151, 104571, 2023") as well as the original manuscript "Early waring signals based on critical slowing down have little practical applicability to anticipate honeybee health collapse"
