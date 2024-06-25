import pandas as pd

def parse_T_fluxes(filepath):
    '''
    Read in LAKE model methane file 'methane_series 1 1.dat'

    Args:
    filepath (string): full filepath to methane file

    Output:
    methane (pd.DataFrame)
    '''
    t_flux=pd.read_csv(filepath, delimiter=r"\s+", skiprows=35, index_col=None, header=None)
    t_flux.columns=['year','month','day','hour','integration_time',
                    'surface temperature, C',
                    'water skin temperature, C',
                    'water surface temperature, C',
                    'mean temperature of water coloumn, C',
                    'maximal temperature in the water coloumn, C',
                    'zero-dimensional model temperature, C',
                    'upper ice surface temperature, C',
                    'upper snow surface temperature, C',
                    'sensible heat flux,    W/m**2',
                    'latent heat flux,      W/m**2',
                    'downward heat flux at the upper lake surface, W/m**2',
                    'downward heat flux at the lake bottom, W/m**2',
                    'friction velocity at the surface (waterside), m/s',
                    'friction velocity at the bottom, m/s',
                    'wind work at the water surface, W/m**2',
                    'albedo of the lake-atmosphere interface, n/d',
                    'shortwave radiation penetrated below surface, W/m**2',
                    'significant wave height, m',
                    'bottom ice salinity, kg/kg',
                    'discharge in x direction, m**3/s',
                    'discharge in y direction, m**3/s',
                    'toptsoil_sc1', 'toptsoil_sc2', 'toptsoil_sc3', 'toptsoil_sc4', 'toptsoil_sc5'] #mg/(m**2*day)
    
    t_flux['Date'] = pd.to_datetime({'Year': t_flux['year'], 'Month': t_flux['month'], 'Day':t_flux['day']})

    return t_flux

def parse_methane_series(filepath):
    '''
    Read in LAKE model methane file 'methane_series 1 1.dat'

    Args:
    filepath (string): full filepath to methane file

    Output:
    methane (pd.DataFrame)
    '''
    methane=pd.read_csv(filepath, delimiter=r"\s+", skiprows=35, index_col=None, header=None)
    methane.columns=['year', 'month', 'day', 'hour', 'integration_time',  'talik depth, m',
                    'lake surface methane concentration, mol/m**3',
                    'lake bottom methane concentration, mol/m**3',
                    'soil bottom methane concentration, mol/m**3',
                    'lake surface oxygen concentration, mol/m**3',
                    'lake bottom oxygen concentration, mol/m**3',
                    'total methane production due to young C decomposition, mol/(m**2*s)',
                    'total methane production due to old C decomposition, mol/(m**2*s)',
                    'methane ebullition flux averaged over the lake bottom, upwards, mol/(m**2*s)',
                    'methane ebullition flux at the surface, mol/(m**2*s)',
                    'methane plant-mediated flux at the lake bottom, mol/(m**2*s)',
                    'methane diffusion flux averaged over the lake bottom, upwards, mol/(m**2*s)',
                    'methane turbulent flux at the lake surface, upwards, mol/(m**2*s)',
                    'methane ebullition flux averaged over the lake bottom, upwards, mg/(m**2*day)',
                    'methane ebullition flux at the surface, mg/(m**2*day)',
                    'methane plant-mediated flux at the lake bottom, mg/(m**2*day)',
                    'methane diffusion flux averaged over lake bottom, upwards, mg/(m**2*day)',
                    'methane turbulent flux at the lake surface, upwards, mg/(m**2*day)',
                    'methane turbulent flux at the bottom of mixed layer normalized by  surface area, upwards, mg/(m**2*day)',
                    'methane flux from sediments in the mixed layer normalized by  surface area, upwards, mg/(m**2*day)',
                    'methane bubble flux at the bottom of the mixed layer normalized by  surface area, upwards, mg/(m**2*day)',
                    'total methane oxidation in water normalized by  surface area, mg/(m**2*day)',
                    'methane oxidation in mixed layer normalized by  surface area, mg/(m**2*day)',
                    'co2 turbulent flux at the lake surface, upwards, mol/(m**2*s)',
                    'co2 ebullition flux at the surface, mol/(m**2*s)',
                    'oxygen turbulent flux at the lake surface, upwards, mol/(m**2*s)',
                    'oxygen ebullition flux at the surface, mol/(m**2*s)',
                    'methane flux by inlets, normalized by surface area, mol/(m**2*s)',
                    'methane flux through outlet, normalized by surface area, mol/(m**2*s)',
                    'methane_ebul_sc1', 'methane_ebul_sc2', 'methane_ebul_sc3', 'methane_ebul_sc4', 'methane_ebul_sc5'] #mg/(m**2*day)
    
    methane['Date'] = pd.to_datetime({'Year': methane['year'], 'Month': methane['month'], 'Day':methane['day']})
    methane = methane[['Date','year', 'month', 'day', 'hour',
                       'talik depth, m',
                       'lake surface methane concentration, mol/m**3',
                       'lake bottom methane concentration, mol/m**3',
                       'total methane production due to young C decomposition, mol/(m**2*s)',
                       'total methane production due to old C decomposition, mol/(m**2*s)', 
                       'methane diffusion flux averaged over the lake bottom, upwards, mol/(m**2*s)',
                       'methane ebullition flux at the surface, mol/(m**2*s)',
                       'methane ebullition flux at the surface, mg/(m**2*day)',
                       'methane turbulent flux at the lake surface, upwards, mol/(m**2*s)',
                       'methane turbulent flux at the lake surface, upwards, mg/(m**2*day)',
                       'total methane oxidation in water normalized by  surface area, mg/(m**2*day)',
                       'methane plant-mediated flux at the lake bottom, mg/(m**2*day)',
                       'co2 turbulent flux at the lake surface, upwards, mol/(m**2*s)',
                       'co2 ebullition flux at the surface, mol/(m**2*s)',
                       'oxygen turbulent flux at the lake surface, upwards, mol/(m**2*s)',
                       'oxygen ebullition flux at the surface, mol/(m**2*s)',
                       'methane_ebul_sc1', 'methane_ebul_sc2', 'methane_ebul_sc3', 'methane_ebul_sc4', 'methane_ebul_sc5']]
    
    methane.columns = ['Date', 'year', 'month', 'day', 'hour',
                       'talik depth, m',
                       'lake surface methane concentration, mol/m**3',
                       'lake bottom methane concentration, mol/m**3',
                       'methane_prod_young_mol/(m**2*s)', 
                       'methane_prod_old_mol/(m**2*s)', 
                       'methane_diffuse_bot_mol/(m**2*s)',
                       'methane_ebul_mol/(m**2*s)', 
                       'methane_ebul_mg/(m**2*day)', 
                       'methane_turb_flux_mol/(m**2*s)', 
                       'methane_turb_flux_mg/(m**2*day)',
                       'methane_oxid_mg/(m**2*day)', 
                       'methane_plant_med_flux_mg/(m**2*day)',
                       'co2_turb_flux_mol/(m**2*s)', 
                       'co2_ebul_mol/(m**2*s)', 
                       'ox_turb_flux_mol/(m**2*s)', 
                       'ox_ebul_mol/(m**2*s)',
                      'methane_ebul_sc1', 'methane_ebul_sc2', 'methane_ebul_sc3', 'methane_ebul_sc4', 'methane_ebul_sc5']
    
    methane['run_name'] = filepath.split('/')[-2]

    return methane
    
def parse_layer_file(filepath):
    '''
    Read in LAKE model layers file 'layers 1 1.dat'

    Args:
    filepath (string): full filepath to layers file

    Output:
    layers (pd.DataFrame)
    '''
    layers=pd.read_csv(filepath, delimiter=r"\s+", skiprows=19, index_col=None, header=None)
    layers.columns=['year', 'month', 'day', 'hour', 'integration_time', 'water layer thickness, m', 
                           'W mixed layer thickness, m', 'E mixed layer thickness, m', 'S mixed layer thickness, m', 'N mixed layer thickness, m',
                           'W lower layer thickness, m', 'E lower layer thickness, m', 'S lower layer thickness, m', 'N lower layer thickness, m',
                           'ice layer thickness,   m', 'snow layer thickness,  m', 'bottom ice thickness,  m', 'reservoir volume,  m**3', 'volume deficit (accumulated),  m**3']
    layers['Date'] = pd.to_datetime({'Year': layers['year'], 'Month': layers['month'], 'Day':layers['day']})
    layers['mean_mixed_layer_thickness'] = layers[['W mixed layer thickness, m', 'E mixed layer thickness, m', 'S mixed layer thickness, m', 'N mixed layer thickness, m']].mean(axis=1)
    layers['mean_lower_layer_thickness'] = layers[['W lower layer thickness, m', 'E lower layer thickness, m', 'S lower layer thickness, m', 'N lower layer thickness, m']].mean(axis=1)
    layers = layers.drop(columns = ['W mixed layer thickness, m', 'E mixed layer thickness, m', 'S mixed layer thickness, m', 'N mixed layer thickness, m',
                                    'W lower layer thickness, m', 'E lower layer thickness, m', 'S lower layer thickness, m', 'N lower layer thickness, m'])
    layers['run_name'] = filepath.split('/')[-2]

    return layers

def parse_univariate_file(filepath, variable_name):
    '''
    Read in LAKE model univariate timeseries file in second format '{variable} 1 1f2.dat'. Outputs variable time series to pandas DataFrame.
    Expects first five columns are 'year', 'month', 'day', 'hour', and 'integration_time'
    --------------------------------------------------------------------------------------
    Tested for
    water temperature: filename = 'water_temp  1  1f2.dat'
    co2 conventation: filename = 'co2_water  1  1f2.dat'
    ch4 concentration: filename = 'methane_water  1  1f2.dat'
    o2 concentration: filename = 'oxygen_water  1  1f2.dat'
    doc concentration: filename = 'DOC  1  1f2.dat'
    --------------------------------------------------------------------------------------

    Args:
    filepath (string): full filepath to time series file
    vaeriable_name (string): name for time series variable column

    Output:
    layers (pd.DataFrame)
    '''
    
    file=pd.read_csv(filepath, delimiter=r"\s+", skiprows=7, index_col=None, header=None)
    file.columns=['year', 'month', 'day', 'hour', 'integration_time', 'depth', variable_name]
    file['Date'] = pd.to_datetime({'Year': file['year'], 'Month': file['month'], 'Day':file['day']})

    return file