# beeEWS
Code for paper on EWS and bee fitness

## Content

### analysisTemperature_low_freq: 
Analyse temperature data for the low-frequency dataset.  
-  1: Loads the respective datasets about Infested and Uninfested colonies, as well as External temperature (provided as a single csv files), as well as colony strength (already input manually, data from previous publications as described in the Main Text of the article).  
- 2: Visualizes the data (Temperature and colony strength, over time)
- 3a: a failed test of naive detrending by external temperature (negative result to check that such detrending is not meaningful.  
3b: Change point analysis and Granger causality test, to statistically observe if T-ext influences directly T_int
- 4: Look for law between T-ext and T_int, using a polynomial fit (calls **FitTemp.m** function).

### analysisTemperature_high_freq: 
Similar as for the previous one, but for the high-frequency dataset. Covers only points 1 and 2.

### EWS_extract_low_freq
Extract and analyse statistical indicators and EWS for the low-frequency dataset, for infested and uninfested hives
- 1: Load data
- 2: Calculate Var, AC(1) and CV after LOESS-detrending on an expanding window
- 3: plot (incl. chage points detected in **analysisTemperature_low_freq**
- 4: Estimate true and false positive signals, for the ROC, by testing various thresholds of kendall's \tau from the trends of the considered indicators
- 5: ROC plot, AUC calculation, estimation of best threshold for \tau

### EWS_extract_high_freq
Extract and analyse statistical indicators and EWS for the high-frequency dataset, for infested and uninfested hives
- 1: Load data. Also remove the daily oscillations via differencing
- 2: Calculate Var, AC(1), CV and spectral slope (here, uses the **estimate_noise_color** function) after LOESS-detrending on an expanding window
- 3: Plot indicators, on the same timeframe as for the low-frequency dataset
- 4: Estimate true and false positive signals, for the ROC, by testing various thresholds of kendall's \tau from the trends of the considered indicators
- 5: ROC plot, AUC calculation, estimation of best threshold for \tau
- 6: noise properties (distribution and frequency properties) before and after the changepoint, after Savitzky-Golay global smoothing. Uses  **estimate_noise_color** and **plot_psd_fit** functions.

### phase_plane_plotter_khoury
Estimates nullclines, vector field, equilibria and their stability, for the model by Khoury et al (2011) analysed in the paper. Plots the results, for different values of $m$; change index of mu_values at line 47 to get the tested values (0.2, 0.4 or $m_c$).

## Utils
The code is in MATLAB. The functions require the symbolic toolbox and the signal processing toolbox. Tested on MATLAB 2025b.

## Data
Data avialability as described in the paper.

## Credits
Code developed by Daniele Proverbio. If you reuse the code in full or in part, please cite the original manuscript: "Early waring signals based on critical slowing down have little
practical applicability to anticipate honeybee health collapse"

