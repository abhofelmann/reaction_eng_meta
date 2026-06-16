# Electrified reaction engineering with metamaterial reactors 
This repository contains data and code used to model and optimize the deisgn of an inductively heated metamaterial reactor. 

## Overview of this repository 
```
reaction_eng_meta
|-code
|-data
```
**code**: contains MATLAB scripts used to fit 1D pseudo-homogeneous model and predict temperature profiles 
- `Correct_sweep_MaxT.mlx`: MATLAB live script to rank reactor configurations by defined FOM, then re-run top 30 structures to maintin max temperature, and re-rank.
- `Cr_infty_wradial_combine_runs.mxl`: MATLAB live script used to combine and rank the predicted temperature profiles, conversion, and power data generated using `Temp_Crinf_wRadial.m`. Files called by this program can be found in `data/model_predict_Crinfty`. 
- `Predict_Temp.m`: MATLAB script containing 1D pseudo-homogeneous model which predicts temperature profile for given flow, max. temperature, and reactor configuraiton
- `Running_TestPoints.mlx`: MATLAB live script which verifies model training by predicting temeprature profile of given test points
- `Temp_Crinf_wRadial.m`: MATLAB script which constains modified version of `Predict_Temp.m` for predicting the temperature profile of reactor configurations with infinite AC resistance ratios, accounting for radial loss
- `compile_data.mxl`: MATLAB live script used to clean and compile 100 training data points. Produces file `Compileddata.mat` in `\experimental_reaction`. 
- `sweep1.m`: MATLAB script used to predict temperature profiles of 2^14 possible reactor configurations.

**data**: contains experimental and simulation data associated with fitting and utilizing the scripts in `code/`
- `experimental_reaction/`: raw and cleaned data which makes up the 100 data points used to train and test 1D pesudo-homogeneous model
  - `Compileddata.mat`: .mat file with temperature, power, conversion, structure, and flow data for each of the 100 data points
  - `RWGS_`: pre-fix for 5 .mat files containing cleaned up expeirmental data for the 5 different reactor configurations used to collect training data. File name is of the form `RWGS_SSSS_TXXX_F_F_F_F_Fslpm.mat` where SSSS is the structure configuration, XXX is the list of maximum temepratures in Celsius, and F is the flow rate in slpm.   
- `model_predict_Cr2.3/`: .mat files containing predicted temperature profiles (using `Predict_Temp.m`) for the six experimentally verified operating conditions found in Figure 4 of the main text. The file names include the flow rate in slpm (fX), the maximum temperature (500 C or 525 C).
-   - `_w2deg_` denotes the first run of the 2^14 structures with maximum temperatures within 2 degrees Celsius of the set temperature
    - `_top30_` denotes the top 30 structures of a given operating condition, used to ensure the maximum temeprature profile is within 0.5 degrees Celsius of the set temperature.
- `model_predict_Crinfty/`: .mat files containing predicted temepratures (using `Temp_Crinf_Radial.m`) for 2^14 reactor configurations with an AC resistance ratio of infinity at a maximum temperature of 500 degrees Celsius and a flow rate of 2500 h^-1. Each of the files were created using the same operating conditions. The numbers indicate the index of the reactor configuration. 



