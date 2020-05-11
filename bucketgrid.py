# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 10:52:25 2017
Version 02.10.18/SL

@author: slauniai
"""
import numpy as np
eps = np.finfo(float).eps

class BucketGrid(object):
    """
    Two-layer soil water bucket model for gridded use in SpaFHy.
    """
    def __init__(self, spara, outputs=False):
        """
        Initializes BucketGrid:
        Args:
            REQUIRED:
            spara - dictionary of soil properties. keys - values np.arrays
                depth [m]
                poros [m3m-3]
                fc [m3m-3]
                wp [m3m-3]
                ksat [ms-1]
                beta [-]
                maxpond [-]
                org_depth [m]
                org_poros [m3m-3]
                org_fw [m3m-3]
                org_rw [m3m-3]
            
                pond_sto - initial pond storage [m]
                org_sat - initial saturation or organic layer [-]
                rootzone_sat - initial saturation of root zone [-]
            OPTIONAL:  
            outputs - True appends output grids to dict stored within object 

        CHANGES:
            05.05.2020 removed typo in watbal mbe computation and added outputs
        """
        
        """ set object properties. All will be 1d or 2d arrays of same shape """
        # above-ground pond storage [m]
        self.MaxPond = spara['maxpond']
        
        # top layer is interception storage, which capacity depends on its depth [m]
        # and field capacity
        self.D_top = spara['org_depth']     # depth, m3 m-3
        self.poros_top = spara['org_poros'] # porosity, m3 m-3
        self.Fc_top = spara['org_fc']       # field capacity m3 m-3
        self.rw_top = spara['org_rw']       # ree parameter m3 m-3
        self.MaxStoTop = self.Fc_top * self.D_top # maximum storage m 

        # root-zone layer properties
        self.D = spara['depth']             # depth, m
        self.poros = spara['poros']         # porosity, m3 m-3     
        self.Fc = spara['fc']               # field capacity, m3 m-3 
        self.Wp = spara['wp']               # wilting point, m3 m-3
        self.Ksat = spara['ksat']           # sat. hydr. cond., m s-1
        self.beta = spara['beta']           # hyd. cond. exponent, -
        # self.soilcode = spara['soilcode']   # soil type integer code
        
        self.MaxWatSto = self.D*self.poros  # maximum soil water storage, m

        """
        set buckets initial state: given as arrays
        """
        self.PondSto = np.minimum(spara['pond_sto'], self.MaxPond)
        
        # toplayer storage and relative conductance for evaporation
        self.WatStoTop = self.MaxStoTop * spara['org_sat']
        self.Wliq_top = self.poros_top *self.WatStoTop / self.MaxStoTop
        self.Ree = np.minimum(self.Wliq_top / self.rw_top, 1.0) # relative ecaporation rate (-)
        
        # root zone storage and relative extractable water
        self.WatSto = np.minimum(spara['rootzone_sat']*self.D*self.poros, self.D*self.poros)
        
        self.Wliq = self.poros*self.WatSto / self.MaxWatSto
        self.Wair = self.poros - self.Wliq
        self.Sat = self.Wliq/self.poros
        self.Rew = np.minimum((self.Wliq - self.Wp) / (self.Fc - self.Wp + eps), 1.0)
        
        # grid total drainage to ground water [m]
        self._drainage_to_gw = 0.0
        
        # create dictionary of empty lists for saving results
        if outputs:
            self.results = {'Infil': [], 'Retflow': [], 'Drain': [], 'Roff': [], 'ET': [],
            'Mbe': [], 'Wliq': [], 'PondSto': [], 'Wliq_top': [], 'Ree': []}

    def watbal(self, dt=1.0, rr=0.0, tr=0.0, evap=0.0, retflow=0.0):
        """
        Computes 2-layer bucket model water balance for one timestep dt
        Top layer is interception storage and contributes only to evap.
        Lower layer is rootzone and contributes only tr and creates drainage.
        Capillary interaction between layers is neglected and connection from bottom up
        is only in case of excess returnflow.
        Pond storage can exist above top layer.
        
        IN:
            dt [s]
            rr = potential infiltration [m]
            tr = transpiration from root zone [m]
            evap = evaporation from top layer [m]
            retflow = return flow from ground water [m]
        OUT: dict with 
            inflow [m] - total inflow to root zone
            roff [m] - surface runoff
            drain [m] - drainage from root zone
            tr [m] - transpiration from root zone
            mbe [m] - mass balance error

        """
        gridshape = np.shape(self.Wliq)  # rows, cols
    
        if np.shape(retflow) != gridshape:
            retflow = retflow * np.ones(gridshape)
        if np.shape(rr) != gridshape:
            rr = rr * np.ones(gridshape)
        
        rr0 = rr.copy()
       
        # add current Pond storage to rr & update storage
        PondSto0 = self.PondSto.copy()
        rr += self.PondSto
        self.PondSto = np.zeros(gridshape)
        
        WatSto0 = self.WatSto.copy()
        WatStoTop0 = self.WatStoTop.copy()
        
        
        #top layer interception & water balance
        interc = np.maximum(0.0, (self.MaxStoTop - self.WatStoTop))\
                    * (1.0 - np.exp(-(rr / self.MaxStoTop)))
        
        self.WatStoTop = np.maximum(0.0, self.WatStoTop + interc)  
        evap = np.minimum(evap, self.WatStoTop)
        self.WatStoTop -= evap
      
        # infiltration to rootzone
        rr = rr - interc
                
        # ********* compute bottom layer (root zone) water balance ***********

        # transpiration removes water from rootzone
        tr = np.minimum(tr, self.WatSto - eps)
        self.WatSto -= tr
        
        # drainage: at gridcells where retflow > 0, set drain to zero.
        # This delays drying of cells which receive water from topmodel storage
        # ... and removes oscillation of water content at those cells.
        drain = np.minimum(self.hydrCond() * dt, np.maximum(0.0, (self.Wliq - self.Fc))*self.D)
        drain[retflow > 0.0] = 0.0
        
        # inflow to root zone: restricted by potential inflow or available pore space
        Qin = (retflow + rr)  # m, pot. inflow
        inflow = np.minimum(Qin, self.MaxWatSto - self.WatSto + drain)
        
        dSto = (inflow - drain)
        self.WatSto = np.minimum(self.MaxWatSto, np.maximum(self.WatSto + dSto, eps))
                
        # if inflow excess after filling rootzone, update first top layer storage
        exfil = Qin - inflow
        to_top_layer = np.minimum(exfil, self.MaxStoTop - self.WatStoTop - eps)
        # self.WatStoTop = self.WatStoTop + to_top_layer
        self.WatStoTop += to_top_layer
        
        # ... and then pond storage ...
        to_pond = np.minimum(exfil - to_top_layer, self.MaxPond - self.PondSto - eps)
        self.PondSto += to_pond
 
        # ... and route remaining to surface runoff
        roff = exfil - to_top_layer - to_pond
        
        # compute diagnostic state variables at root zone:
        self.setState()
        
        # update grid total drainage to ground water [m]
        self._drainage_to_gw = np.nansum(drain)
        
        # mass balance error [m]
        mbe = (self.WatSto - WatSto0)  + (self.WatStoTop - WatStoTop0) + (self.PondSto - PondSto0) \
            - (rr0 + retflow - tr - evap - drain - roff)
        
        # append results to lists; use only for testing small grids!
        if hasattr(self, 'results'):
            self.results['Infil'].append(inflow - retflow)   # infiltration through top boundary
            self.results['Retflow'].append(retflow)         # return flow from below 
            self.results['Roff'].append(roff)       # surface runoff
            self.results['Drain'].append(drain)     # drainage
            self.results['ET'].append(tr + evap)            
            self.results['Mbe'].append(mbe)
            self.results['Wliq'].append(self.Wliq)
            self.results['PondSto'].append(self.PondSto)
            self.results['Wliq_top'].append(self.Wliq_top)
            self.results['Ree'].append(self.Ree)
        
        return inflow, roff, drain, tr, evap, mbe
    
    def setState(self):
        """ updates state variables"""
        # root zone
        self.Wliq = self.poros*self.WatSto / self.MaxWatSto
        self.Wair = self.poros - self.Wliq
        self.Sat = self.Wliq / self.poros
        self.Rew = np.maximum(0.0,
              np.minimum((self.Wliq - self.Wp) / (self.Fc - self.Wp + eps), 1.0))
        
        # organic top layer; maximum that can be hold is Fc
        self.Wliq_top = self.Fc_top * self.WatStoTop / self.MaxStoTop
        self.Ree = self.relative_evaporation()
        
    def hydrCond(self):
        """
        returns hydraulic conductivity [ms-1] based on Campbell -formulation
        """
        k = self.Ksat*self.Sat**(2*self.beta + 3.0)
        return k

    def relative_evaporation(self):
        """
        returns relative evaporation rate from the organic top layer; loosely
        based on Launiainen et al. 2015 Ecol. Mod. Moss-module
        Returns:
            f - [-], array or grid of 
        """
        f = np.maximum(0.0, np.minimum(0.98*self.Wliq_top / self.rw_top, 1.0))
        return f
