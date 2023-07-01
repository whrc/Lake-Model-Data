[![DOI](https://zenodo.org/badge/420248298.svg)](https://zenodo.org/badge/latestdoi/420248298)

# LAKE2.0_processing scripts
 Processing scripts for LAKE2.0 model input and output
 Lake2.0 model outputs and lake temperature validation datasets.

scripts\
	
	Pre-processing
		Toolik
			scripts\formatToolikhrly_final.R
				Used to format toolik met data for model
			scripts\formatToolikLakeInterpolate_final.R
				Used to combine measured lake temp data and interpolate to regular 1m grid
			scripts\ToolikDischarge_final.R
				Used to process raw discharge measurements to inflow and outflow files for model
		FoxDen
			scripts\formatFoxDenDay_final.R
				Used to format foxden met data for model
		Atqasuk
			scripts\Atqasuk\USGS_SM_Atqasuk_data_preprocess-v0.ipynb
				Used to format and calibrate met data for model
		All
			scripts\standardizeLakeTempData.R
				Used to standardize measured lake water temperature data for use with post processing scripts.
	
	Post-processing
		All
			scripts\plotLaketemp_Final.R
				Plot LAKE model temp data compared to obs data and met data. Figs 1-4 in paper. 
			scripts\plotLakeTempScenariosCM_final.R
				Calculate errors for scenarios, Combine errors for color matrix plot. Fig 5
			scripts\plotLakeTempScenariosCM_season_final.R	
				Calculate errors for scenarios with seasons, Combine errors for seasonal barplot. Fig 6

model output\
	
	Toolik
		model output\Toolik\layers 1 1.dat
		model output\Toolik\soil_temp 1 1.dat
		model output\Toolik\water_temp 1 1.dat
	
	FoxDen
		model output\FoxDen\layers 1 1.dat
		model output\FoxDen\soil_temp 1 1.dat
		model output\FoxDen\water_temp 1 1.dat
	
	Atqasuk
		model output\Atqasuk\layers 1 1.dat
		model output\Atqasuk\soil_temp 1 1.dat
		model output\Atqasuk\water_temp 1 1.dat

validation data\
	
	validation data\Toolik_LakeTemp.csv
	validation data\FoxDen_LakeTemp.csv
	validation data\Atqasuk_LakeTemp.csv
