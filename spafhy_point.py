# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 13:13:45 2017

@author: slauniai

spafhy_point: 
    combines canopygrid and bucketgid for solving point-scale water balance
v. 260618 / Samuli
"""
import numpy as np
import pandas as pd
from canopygrid import CanopyGrid
from bucketgrid import BucketGrid

eps = np.finfo(float).eps  # machine epsilon

 
class SpaFHy_point():
    """
    SpaFHy for point-scale simulation. Couples Canopygrid and Bucketdrid -classes
    """
    def __init__(self, pgen, pcpy, pbu, FORC, cmask=np.ones(1),
                 cpy_outputs=True, bu_outputs=True):
        """
        creates SpaFHy_point -object
        Args:
            pgen - parameter dict
            pcpy - canopy parameter dict
            pbu - bucket model parameter dict
            FORC - forcing data (pd.DataFrame)
            cmask - catchment mask; in case multiple cells are simulated. np.array
            cpy_outputs - saves canopy outputs within self.cpy.results
            bu_outputs - saves bucket outputs within self.cpy.results
        Returns:
            object
        """
        self.dt = pgen['dt']  # s
        self.id = pgen['catchment_id']
        self.spinup_end = pgen['spinup_end']
        self.pgen = pgen
    
        self.FORC = FORC 
        self.Nsteps = len(self.FORC)

        """--- initialize CanopyGrid and BucketGrid ---"""
        
        # sub-models require states as np.array; so multiply with 'cmask'
        cstate = pcpy['state'].copy()

        for key in cstate.keys():
            cstate[key] *= cmask

        self.cpy = CanopyGrid(pcpy, cstate, outputs=cpy_outputs)

        for key in pbu.keys():
            pbu[key] *= cmask
            
        self.bu = BucketGrid(pbu, outputs=bu_outputs)

                          
    def _run(self, fstep, Nsteps, soil_feedbacks=True, results=False):
        """ 
        Runs SpaFHy_point
        IN:
            fstep - index of starting point [int]
            Nsteps - number of timesteps [int]
            soil_feedbacks - False sets REW and REE = 1 and ignores feedback 
                        from soil state to Transpi and Efloor
            results - True returns results dict
        OUT:
            updated state,  optionally saves results within self.cpy.results and
            self.bu.results and/or returns at each timestep as dict 'res'
        """
        dt = self.dt

        for k in range(fstep, fstep + Nsteps):
            print('k=' + str(k))

            # forcing
            doy = self.FORC['doy'].iloc[k]; ta = self.FORC['T'].iloc[k]
            vpd = self.FORC['VPD'].iloc[k]; rg = self.FORC['Rg'].iloc[k]
            par = self.FORC['Par'].iloc[k]; prec = self.FORC['Prec'].iloc[k]
            co2 = self.FORC['CO2'].iloc[k]; u = self.FORC['U'].iloc[k]
            if not np.isfinite(u):
                u = 2.0            

            if soil_feedbacks:
                beta0 = self.bu.Ree # affects surface evaporation
                rew0 = self.bu.Rew # affects transpiration
            else:
                beta0 = 1.0
                rew0 = 1.0

            # run CanopyGrid
            potinf, trfall, interc, evap, et, transpi, efloor, mbe = \
                self.cpy.run_timestep(doy, dt, ta, prec, rg, par, vpd, U=u, CO2=co2,
                                      beta=beta0, Rew=rew0, P=101300.0)

            # run BucketGrid water balance
            infi, infi_ex, drain, tr, eva, mbes = self.bu.watbal(dt=dt, rr=1e-3*potinf, tr=1e-3*transpi,
                                                           evap=1e-3*efloor, retflow=0.0)

            if results == True:
                # returns fluxes and and state in dictionary
                res = {'Wliq': self.bu.Wliq, 'Wliq_top': self.bu.Wliq_top,
                       'Evap': evap, 'Transpi': transpi, 'Efloor': efloor,
                       'Interc': interc, 'Drain': 1e3*drain, 'Infil': 1e3*infi, 
                       'Qs': 1e3*infi_ex
                       }
                return res
            
