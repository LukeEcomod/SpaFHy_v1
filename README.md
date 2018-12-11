# SpaFHy v1.0
Spatial Forest Hydrology (Launiainen et al. 2019)

This is repository for SpaFHy_v1 to be submitted into Geosci. Model Dev.

Contains working examples of catchment scale and point-scale models (with example data).

To run the model:

1) clone repository
2) create folder '/results'
3) change file paths accordingly in spafhy_parameters

4) To run catchment-scale model, see demo_spafhy_C3. append parent folder to sys.path if necessary
5) To run point-scale model, see demo_spafhy_point. append parent folder to sys.path if necessary

Folders:

/data - contains necessary example data

/misc - old and obsolete stuff; may be re-cyclable

/Figs_C3 - figures etc.

/results - (create!) to store results.

Tested in Py 3.6, should also work with 2.7.
Required packages: numpy, pandas, NetCDF4, pickle, scipy, matplotlib, (seaborn for plotting)
Pandas is constantly evolving so depending on version some fixes may be necessary.

Antti: can you check through the codes and required versions of packages. Also, test using 2.7

Contact: samuli.launiainen@luke.fi
