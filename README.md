# V2Karst_model

V2Karst is a large-scale conceptual semi-distributed model that simulates sub-daily (from V1.1) and daily potential groundwater recharge in karst regions. The model includes an explicit representation of both karst and land cover properties, which allows assessment of the combined impact of future changes in climate and land cover on karst groundwater recharge. The model builds on its previous version (Varkarst), which had a simple formulation of evapotranspiration and did not represent explicitly land cover properties.

A complete description of VarKarst is provided in [Hartmann et al. (2015)](https://doi.org/10.5194/gmd-8-1729-2015) and a description of V2Karst in [Sarrazin et al. (2018)](https://gmd.copernicus.org/articles/11/4933/2018/). If you make use of the V2Karst, please acknowledge these two articles.

V2Karst is provided under the terms of the GNU General Public License version 3.0 and therefore **without warranty of any kind**.

The V2Karst Model is written in matlab and is composed of 6 functions. This division of the model in separate modules facilitates future model improvements and the testing of alternative model equations by users.

An overview of the structure of the V2Karst model is provided below, while further information on the usage of each function is provided in the function 'help'. Additionally, in-code comments explain the model equations.

**- V2Karst.m**<br />
This function simulates the V2Karst model and calls the other 5 functions (interception_routine.m, soil_epikarst_routine.m, Penman_Monteith.m, aerodynamic_resistance.m and LAI_seasonality.m)
**This is the function that users need to call to simulate V2Karst.**

**- interception_routine.m**<br />
The function evaluates daily/sub-daily actual evaporation from canopy interception.

**- soil_epikarst_routine.m**
This function evaluates the daily/sub-daily soil and epikarst water balance.

**- Penman_Monteith.m**<br />
This function evaluates the daily/sub-daily potential evapotranspiration using the Penman Monteith equation.

**- aerodynamic_resistance.m**<br />
This function evaluates the daily/sub-daily aerodynamic resistance.

**- LAI_seasonality.m**<br />
This function evaluates the daily/sub-daily value of LAI.

**- snow_routine.m**<br />
This function evaluates the daily/sub-daily effective precipitation (Pre_eff) from the variation of the snow pack and set potential transpiration and potential soil evaporation to 0 in presence of snow pack.

**- net_radiation.m**<br />
This function function evaluated net radiation from downward longwave and shortwave radiation, surface albedo and air temperature as a proxy for surface temperature.
The assumption that air temperature is a proxy for surface temperature has only been tested for daily time step. **When net radiation has to be calculated from downward radiations, we recommend to use a daily time step (and NOT sub-daily) to simulate V2Karst**.

Additional functions are provided to process the input data, as reported below.

**- RH_WS_blending_height.m**<br />
This function computes relative daily/sub-daily humidity and wind speed at a blending height. This is useful when wind speed and humidity measurements are provided below canopy level.
**This is the function that users need to call to transform humidity and wind speed data.**

**- RH_blending_height.m**<br />
This function computes relative humidity at a blending height.

**- WS_blending_height.m**<br />
This function computes wind speed at a blending height.

DOI for latest version (V1.1):
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1484282.svg)](https://doi.org/10.5281/zenodo.1484282)
