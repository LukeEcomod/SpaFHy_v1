%HydeDaily2000-2010.txt'file info

% 
% Samuli Launiainen, Luke 7.12.2016
%
% FOR CANOPYGRID / SPATHY CALIBRATION IN PYTHON
% 
% Contains daily sums or daily averages of Hyytiälä SMEAR II-measurements (Kolari et al., 2009 Global Change Biol.; Launiainen 2010 Biogeosciences)
% 
% Missing data = NaN
%
% For data use permissions see: https://avaa.tdata.fi/web/smart/smear/terms-of-use
% To obtain data from SMEAR -stations run by Institute of Atmospheric and Earth System Research (INAR) of University of Helsinki, visit
% https://avaa.tdata.fi/web/smart/smear/download
%
% Gap-filling & Reco separation as in Kolari et al. (2009) Bor. Env. Res.
% 
% ET and sensible heat flux Snow.Data(:,1:7) series are gap-filled as in Reichstein et al. 2005 using combination of look-up table methods 
% and mean diurnal variability.
% 
% Gap-filling of meteorological data (Tair, PAR, Rnet, CO2) by linear interpolation.
% 
% Columns
% 
% 1 yyyy, 2 mm, 3 dd
% 4 doy
% 5 NEE (g C m-2 d-1)
% 6 GPP (g C m-2 d-1)
% 7 TER (g C m-2 d-1)
% 8 ET  (mm d-1)
% 9 H	  (W d-1) daily ave
% 10 NEEflag (-); fraction of gapfilled data. 0 = all are measured, 1 = all are gap-filled
% 11 ETflag
% 12 Hflag
% 13 PAR (umol m-2 s-1) daily ave
% 14 Rnet (W m-2) daily ave
% 15 Tair (degC) daily ave
% 16 VPD (kPa) daily ave
% 17 CO2 (ppm) daily ave, from SMEAR II
% 18 Precip (mm d-1) from SMEAR II
% 19 Precip (mm d-1) from Hyytiälä FMI station
% 20 Wind speed (m s-1) ave
% 21 Pressure (kPa) from SMEAR II
% 22 Snow water equivalent (mm); total amount of water stored in snowpack at the ground
% 23 SWC (soil water content, m3 m-3) Humus
% 24 SWC (soil water content, m3 m-3) A-horizon
% 25 SWC (soil water content, m3 m-3) B-horizon
% 26 SWC (soil water content, m3 m-3) C-horizon
% 27 Tsoil (deg C) Humus
% 28 Tsoil (deg C)  A-horizon
% 29 Tsoil (deg C)  B-horizon
% 30 Tsoil (deg C)  C-horizon
% 31 Rnetflag (-); fraction of gapfilled data. 0 = all are measured, 1 = all are gap-filled
% 32 Trfall (mm/d)
% 33 Snowdepth (cm)
% 34 Snowdepth std (cm)
% 35 SWE (mm)
% 36 SWE std (mm)
% 37 Roff1 (mm/d)
% 38 Roff2 (mm/d)

%HydeCage.... .txt readme

Columns
1 yyyy, 2 mm, 3 dd
4 doy
5 NEE (g C m-2 d-1)
6 GPP (g C m-2 d-1)
7 TER (g C m-2 d-1)
8 ET  (mm d-1)
9 H	  (W d-1) daily ave
10 NEEflag (-); fraction of gapfilled data. 0 = all are measured, 1 = all are gap-filled
11 ETflag
12 Hflag
13 PAR (umol m-2 s-1) daily ave
14 Rnet (W m-2) daily ave
15 Tair (degC) daily ave
16 VPD (kPa) daily ave
17 CO2 (ppm) daily ave, from SMEAR II
18 Soil water cont. (m3 m-3) at 0-10 cm, daily ave
19 Precip (mm d-1) from SMEAR II
20 Precip (mm d-1) from Hyytiälä FMI station
21 Wind speed (m s-1) ave
22 Pressure (kPa) from SMEAR II
