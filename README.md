# SpaFHy v1.0

This is repository for SpaFHy_v1 source code (under MIT copyright licence). Contains working examples of catchment scale and point-scale models (with example data).

***Reference:***

Launiainen, S., Guan, M., Salmivaara, A., and Kieloaho, A.-J.: Modeling boreal forest evapotranspiration and water balance at stand and catchment scales: a spatial approach, Hydrol. Earth Syst. Sci., https://www.hydrol-earth-syst-sci.net/23/3457/2019/ 

***Abstract:***

Vegetation is known to have strong influence on
evapotranspiration (ET), a major component of terrestrial
water balance. Yet hydrological models often describe ET
by methods unable to include the variability of vegetation
characteristics in their predictions. To take advantage of the
increasing availability of high-resolution open GIS data on
land use, vegetation and soil characteristics in the boreal
zone, a modular, spatially distributed model for predicting
ET and other hydrological processes from grid cell to catchment
level is presented and validated. An improved approach
to upscale stomatal conductance to canopy scale using information
on plant type (conifer/deciduous) and stand leafarea
index (LAI) is proposed by coupling a common leafscale
stomatal conductance model with a simple canopy radiation
transfer scheme. Further, a generic parametrization for
vegetation-related hydrological processes for Nordic boreal
forests is derived based on literature and data from a boreal
FluxNet site. With the generic parametrization, the model
was shown to reproduce daily ET measured using an eddycovariance
technique well at 10 conifer-dominated Nordic
forests whose LAI ranged from 0.2 to 6.8m2 m-2. Topography,
soil and vegetation properties at 21 small boreal headwater
catchments in Finland were derived from open GIS data
at 16m grid size to upscale water balance from stand
to catchment level. The predictions of annual ET and specific
discharge were successful in all catchments, located from 60
to 68N, and daily discharge was also reasonably well predicted
by calibrating only one parameter against discharge
measurements. The role of vegetation heterogeneity in soil
moisture and partitioning of ET was demonstrated. The proposed
framework can support, for example, forest trafficability
forecasting and predicting impacts of climate change
and forest management on stand and catchment water balance.
With appropriate parametrization it can be generalized
outside the boreal coniferous forests.

***Drained peatland version (https://github.com/LukeEcomod/SpaFHy-Peat):***

Leppä, K., Sarkkola, S., Peltoniemi, M., Hökkä, H., Saarinen, M., Lehtonen, A., Laiho, R., Mäkipää, R., Launiainen, S., Nieminen, M. 2020. Selection cuttings as a tool to control water table level in boreal drained peatland forests. Front. Earth Sci., 09 October 2020, https://doi.org/10.3389/feart.2020.576510


***Applications:***

Tyystjärvi V., Kemppinen J., Luoto M., Aalto T., Markkanen T., Launiainen S., Kieloaho A-J, Aalto J. 2021. Modelling spatio-temporal soil moisture dynamics in mountain tundra. Hydrological Processes.2022;36:e14450https://doi.org/10.1002/hyp.14450

Salmivaara, A., Launiainen S., et al. 2020. Towards dynamic Forest Trafficability Prediction using Open Spatial Data, Hydrological Modelling and Sensor Technology, Forestry 2020; 1–13, doi:10.1093/forestry/cpaa010

***Minimal User-guide to run the model:***

1) clone repository
2) create folder '/results'
3) change file paths accordingly in spafhy_parameters

4) To run catchment-scale model, see demo_spafhy_C3. append parent folder to sys.path if necessary
5) To run point-scale model, see demo_spafhy_point. append parent folder to sys.path if necessary

Folders:

/data - contains necessary example data

/misc - old and obsolete stuff; may be re-cyclable

/FigsC3 - figures etc.

/results - (create manually!) to store results.

Tested in Py3.6, should also work with 2.7.
Required packages: os, sys, numpy, pandas, NetCDF4, pickle, scipy, matplotlib, (seaborn for plotting)
Pandas is constantly evolving so depending on version some fixes may be necessary.


***Contact:*** samuli.launiainen@luke.fi
