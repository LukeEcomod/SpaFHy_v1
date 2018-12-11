#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 13:48:49 2017

@author: MG
"""
import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from spathy_calibration import sve_calibrations
# import matplotlib.lines as mlines
# import matplotlib.gridspec as gridspec
# from src import analysis as at
# from scipy.stats import lognorm
# from datetime import datetime, timedelta

# from iotools import read_setup, create_catchment, read_SVE_runoff,read_FMI_weather,read_AsciiGrid, read_climate_prj

# spathy_path = os.path.join(os.path.expanduser('~'),'projects','spathy') #path to spathy folder

eps = np.finfo(float).eps

spathy_path = os.path.join('c:', 'c:\datat\spathydata')
results_path = os.path.join(spathy_path, 'Results', 'Cal')
setupfile = os.path.join(r'c:\repositories\spathy\ini', 'spathy_default.ini')
Nreps = 50

chm=[['1', '2013-01-01', '2015-12-31', '2013-12-31'],   # lompolojanganoja 514 ha
     ['2', '2006-01-01', '2009-12-31', '2006-12-31'],   # liuhapuro 170 ha
     ['3', '2008-01-01', '2015-12-31', '2008-12-31'],   # porkkavaara 72 ha
     ['10', '2011-01-01', '2013-12-31', '2011-12-31'],  # kelopuro 74 ha. 2014 gappy, 2015 runoff is low
     ['11', '2014-01-01', '2015-12-31', '2014-12-31'],  # hauklammenoja 137 ha
     ['13', '2014-01-01', '2015-12-31', '2014-12-31'],  # rudbacken 436 ha
     ['14', '2011-01-01', '2015-12-31', '2011-12-31'],  # paunulanpuro 154 ha
     ['16', '2011-01-01', '2015-12-31', '2011-12-31'],  # huhtisuonoja 500 ha. very flat, large fraction is drained peatlands
     ['17', '2006-01-01', '2009-12-31', '2006-12-31'],  # kesselinpuro 2100 ha
#     ['18','2011-01-01', '2015-12-31', '2011-12-31'],  # korpijoki, area 12200 ha so not suitable
     ['19', '2011-01-01', '2015-12-31', '2011-12-31'],  # pahkaoja 2344 ha
     ['20', '2011-01-01', '2015-12-31', '2011-12-31'],  # vaarajoki 1900 ha
     ['21', '2011-01-01', '2015-12-31', '2011-12-31'],  # myllypuro 1053 ha
     ['22', '2011-01-01', '2015-12-31', '2011-12-31'],  # vaha-askanjoki 1600 ha
#    [ '23','2011-01-01', '2015-12-31', '2011-12-31'],  # ylijoki 5600 ha, very large and slow
     ['24', '2011-01-01', '2015-12-31', '2011-12-31'],  # kotioja 1800 ha
     ['25', '2011-01-01', '2015-12-31', '2011-12-31'],  # kohisevanpuro 1070 ha
     ['26', '2011-01-01', '2015-12-31', '2011-12-31'],  # iittovuoma 1160 ha
     ['27', '2011-01-01', '2015-12-31', '2011-12-31'],  # laanioja 1362 ha
     ['28', '2013-01-01', '2015-12-31', '2013-12-31'],  # kroopinsuo 179 ha
     ['29', '2012-01-01', '2015-12-31', '2012-12-31'],  # surnui 71 ha, poor data quality
     ['30', '2011-01-01', '2015-12-31', '2011-12-31'],  # pakopirtti 795 ha, uncertain catchment boundaries
     ['31', '2011-01-01', '2015-12-31', '2011-12-31'],  # ojakorpi 33 ha
     ['32', '2011-01-01', '2015-12-31', '2011-12-31'],  # rantainrahka 38 ha
     ['33', '2011-01-01', '2012-12-31', '2011-12-31'],  # kivipuro 54 ha
     ]
#subset = [0,1,2,3,4,5,6,7,8,11,13,14,15,16,17,20,21,22]#
subset = [9, 10, 12, 18, 19] # these were missing!
# marker = ['o','.',',','v','>','*','h','s','D','p','o','.',',','v','>','*','h','s','D','p']
# ids = np.empty(len(subset))
# subset = [1]# , 2, 3, 4, 6]
#subset = [0, 1, 2, 3, 4, 5, 6, 7, 17, 18, 20, 21, 22]

#subset = range(0, len(chm))
for k in subset:
    print chm[k]
    cid = chm[k][0]
    start = chm[k][1]
    end = chm[k][2]
    spinup_end = chm[k][3]

    # run full calibration
    spot, res = sve_calibrations(setupfile, cid, start, end, spinup_end, reps=Nreps)

#    # run only topmodel calibration
#   _, _ = sve_topmodel_calibration(fn, cid, start, end, spinup_end, reps=Nreps)

