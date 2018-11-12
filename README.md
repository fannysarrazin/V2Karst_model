# V2Karst_model

V2Karst is a large-scale conceptual semi-distributed model that simulates sub-daily (from V1.1) and daily potential groundwater recharge in karst regions. The model includes an explicit representation of both karst and land cover properties, which allows assessment of the combined impact of future changes in climate and land cover on karst groundwater recharge. The model builds on its previous version (Varkarst), which had a simple formulation of evapotranspiration and did not represent explicitly land cover properties.

A complete description of VarKarst is provided in [Hartmann et al. (2015)](https://doi.org/10.5194/gmd-8-1729-2015) and a description of V2Karst in [Sarrazin et al. (in review)](https://www.geosci-model-dev-discuss.net/gmd-2017-315/). If you make use of the V2Karst, please acknowledge these two articles.

V2Karst is provided under the terms of the GNU General Public License version 3.0 and therefore **without warranty of any kind**.

The V2Karst model is written in matlab and is composed of 6 functions. This division of the model in separate modules facilitates future model improvements and the testing of alternative model equations by users.

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

[![DOI](https://zenodo.org/badge/114650321.svg)](https://zenodo.org/badge/latestdoi/114650321)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1484282.svg)](https://doi.org/10.5281/zenodo.1484282)
