Atqasuk Lake winter precipitation calibration workflow: 

1. Open USGS_SM_Atqasuk_data_preprocess-v0.ipynb and run preprocess_Atqasuk function. It will produce the data/Atqasuk.dat. Note, that to run pre-processing script successfully, make sure data/South Meade.txt and data/atqasuk-atq-2014-2015-meteorology-timeseries-calon.csv are in right folders (see preproccess_Atqasuk_met_data.py) 
2. Run the model. Note that setup and inputs files can be found here https://data.ess-dive.lbl.gov/view/doi:10.15485/1808368. 
3. Use post-proces_Atqasuk_LAKE_results_v0 to plot the results
4. Go back to step 1, change the `wprecip` and `p_reduce` (see `preprocess_Atqasuk` function), do steps 1-3. In step 3, check difference between previous and current results.