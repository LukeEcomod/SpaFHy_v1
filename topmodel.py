# -*- coding: utf-8 -*-
"""
Created on Fri Jun 02 15:14:11 2017

@author: slauniai

******************************************************************************
TopModel (Beven & Kirkby) -implementation for SpatHy -integration

Topmodel() allows spatially varying soil depths and transmissivity
Topmodel_Homogenous() assumes constant properties and hydrologic similarity \n
retermined from TWI = log (a / tan(b))

(C) Samuli Launiainen, 2016-
Last edit: 7.2.2018 / Samuli Launiainen
******************************************************************************
"""

import numpy as np
# import matplotlib.pyplot as plt
eps = np.finfo(float).eps  # machine epsilon

class Topmodel_Homogenous():
    def __init__(self, pp, cellarea, cmask, flowacc, slope, S_initial=None,
                 outputs=False):
        """
        sets up Topmodel for the catchment assuming homogenous
        effective soil depth 'm' and sat. hydr. conductivity 'ko'.
        This is the 'classic' version of Topmodel where hydrologic similarity\
        index is TWI = log(a / tan(b)).
        
        Args:
            pp - parameter dict with keys:
                dt - timestep [s]
                ko - soil transmissivity at saturation [m/s]
                m -  effective soil depth (m), i.e. decay factor of Ksat with depth
                twi_cutoff - max allowed twi -index
                so - initial catchment average saturation deficit (m)
            cmask - catchment mask, 1 = catchment_cell
            cellarea - gridcell area [m2]
            flowacc - flow accumulation per unit contour length (m)
            slope - local slope (deg)
            S_initial - initial storage deficit, overrides that in 'pp'
            outputs - True stores outputs after each timestep into dictionary
        """
        if not S_initial:
            S_initial = pp['so']

        self.dt = float(pp['dt'])
        self.cmask = cmask
        self.CellArea = cellarea
        dx = cellarea**0.5
        self.CatchmentArea = np.size(cmask[cmask == 1])*self.CellArea

        # topography
        self.a = flowacc*cmask  # flow accumulation grid
        self.slope = slope*cmask  # slope (deg) grid

        # effective soil depth [m]
        self.M = pp['m']
        # lat. hydr. conductivity at surface [m2/timestep]
        # self.To = pp['ko']*pp['m']*self.dt
        self.To = pp['ko']*self.dt
        
        """ 
        local and catchment average hydrologic similarity indices (xi, X).
        Set xi > twi_cutoff equal to cutoff value to remove tail of twi-distribution.
        This concerns mainly the stream network cells. 'Outliers' in twi-distribution are
        problem for streamflow prediction
        """
        slope_rad = np.radians(self.slope)  # deg to rad

        xi = np.log(self.a / dx / (np.tan(slope_rad) + eps))
        # apply cutoff
        clim = np.percentile(xi[xi > 0], pp['twi_cutoff'])
        xi[xi > clim] = clim
        self.xi = xi
  
        self.X = 1.0 / self.CatchmentArea*np.nansum(self.xi*self.CellArea)

        # baseflow rate when catchment Smean=0.0
        self.Qo = self.To*np.exp(-self.X)

        # catchment average saturation deficit S [m] is the only state variable
        s = self.local_s(S_initial)
        s[s < 0] = 0.0
        self.S = np.nanmean(s)

        # create dictionary of empty lists for saving results
        if outputs:
            self.results = {'S': [], 'Qb': [], 'Qr': [], 'Qt': [], 'qr': [],
                            'fsat': [], 'Mbe': [], 'R': []
                           }

    def local_s(self, Smean):
        """
        computes local storage deficit s [m] from catchment average
        """
        s = Smean + self.M*(self.X - self.xi)
        return s

    def subsurfaceflow(self):
        """subsurface flow to stream network (per unit catchment area)"""
        Qb = self.Qo*np.exp(-self.S / (self.M + eps))
        return Qb

    def run_timestep(self, R):
        """
        runs a timestep, updates saturation deficit and returns fluxes
        Args:
            R - recharge [m per unit catchment area] during timestep
        OUT:
            Qb - baseflow [m per unit area]
            Qr - returnflow [m per unit area]
            qr - distributed returnflow [m]
            fsat - saturated area fraction [-]
        Note: 
            R is the mean drainage [m] from bucketgrid.
        """

        # initial conditions
        So = self.S
        s = self.local_s(So)

        # subsurface flow, based on initial state
        Qb = self.subsurfaceflow()

        # update storage deficit and check where we have returnflow 
        S = So + Qb - R
        s = self.local_s(S)

        # returnflow grid
        qr = -s
        qr[qr < 0] = 0.0  # returnflow grid, m
        
        # average returnflow per unit area
        Qr = np.nansum(qr)*self.CellArea / self.CatchmentArea

        # now all saturation excess is in Qr so update s and S. 
        # Deficit increases when Qr is removed 
        S = S + Qr
        self.S = S

        # saturated area fraction
        ix = np.where(s <= 0)
        # fsat = np.max(np.shape(ix))*self.CellArea / self.CatchmentArea
        fsat = len(ix[0])*self.CellArea / self.CatchmentArea
        del ix
        
        # check mass balance
        dS = (So - self.S)
        dF = R - Qb - Qr
        mbe = dS - dF

        # append results
        if hasattr(self, 'results'):
            self.results['R'].append(R)
            self.results['S'].append(self.S)
            self.results['Qb'].append(Qb)
            self.results['Qr'].append(Qr)
            self.results['qr'].append(qr)
            self.results['fsat'].append(fsat)
            self.results['Mbe'].append(mbe)

        return Qb, Qr, qr, fsat


#class Topmodel():
#    def __init__(self, pp, ksat=None, cellarea=None, cmask=None, flowacc=None,
#                 slope=None, S_initial=None, outputs=False):
#        """
#        Topmodel, allows for spatially varying 'effective soil depth m' and
#        transmissivity 'to'.
#        Based on Saulnier et al. 1997. J.Hydrol. 202, 158-172.
#        Args:
#            pp - parameter dict with keys {'dt','to','m', 'So'}
#            ksat - saturated hydr. conductivity (ms-1), array or grid
#            cmask - catchment mask, 1 = catchment_cell
#            cellarea - area of each cell (or index class)
#            flowacc - flow accumulation per unit contour length (m)
#            slope - local slope (rad)
#            outputs - True stores outputs after each timestep
#
#            in pp: keys
#            dt - timestep [s]
#            to - soil transmissivity at saturation [m/s]
#            m -  effective soil depth (m), i.e. decay factor of Ksat with depth
#            so - initial catchment average saturation deficit (m)
#        """
#        if not S_initial:
#            S_initial = pp['so']
#
#        self.dt = pp['dt']
#        self.cmask = cmask
#        self.CellArea = cellarea
#        self.CatchmentArea = np.size(cmask[cmask == 1])*self.CellArea
#
#        # topography
#        self.a = flowacc*cmask  # flow accumulation grid
#        self.slope = slope*cmask  # slope (deg) grid
#
#        # local m and average M effective soil depth [m]
#        self.m = pp['m']*cmask
#        self.M = 1.0 / self.CatchmentArea*np.nansum(self.m*self.CellArea)
#
#        self.To = pp['to']*self.dt
#
#        # local xi and catchment average X hydrologic similarity indices
#        rad = np.radians(self.slope)  # deg to rad
#
#        self.xi = self.m*np.log(self.a / (self.To*(np.tan(rad) + eps)))
#        self.X = 1.0 / self.CatchmentArea*np.nansum(self.xi*self.CellArea)
#        # print('X', self.X)
#
#        # baseflow rate when catchment Smean=0.0
#        self.Qo = np.exp(-self.X / self.M)
#
#        # state vars.: local s and catchment average S saturation deficit [m]
#        s = self.local_s(S_initial)
#        s[s < 0] = 0.0
#        self.S = np.nanmean(s)
#
#        # create dictionary of empty lists for saving results
#        if outputs:
#            self.results = {'S': [], 'Qt': [], 'Qb': [], 'Qf': [], 'Roff': []}
#
#    def local_s(self, Smean):
#        """computes local storage deficit s [m] """
#        s = Smean + self.M*(self.X - self.xi)
#        return s
#
#    def subsurfaceflow(self):
#        """subsurface flow to stream network (per unit catchment area)"""
#        Qb = self.Qo*np.exp(-self.S / (self.M + eps))
#        return Qb
#
#    def run_timestep(self, R):
#        """
#        runs a timestep, updates saturation deficit and returns fluxes
#        IN: R - recharge [m per unit catchment area] during timestep
#        OUT:
#            Qb - baseflow [m per unit area]
#            Qf - returnflow [m per unit area]
#            Roff - runoff [m per unit area], recharge to initially sat. area
#            fsat - saturated area fraction [-]
#            qf - returnflow [m per unit area] in matrix
#        """
#
#        # initial conditions
#        So = self.S
#        s = self.local_s(So)
#
#        # subsurface flow, based on initial state
#        Qb = self.subsurfaceflow()
#
#        # for outputs & use of outcommented lines
#        Rec = R
#        Roff = R - Rec
#
##       use this snippet if input R includes drainage also from sat. cells.
##       input to ground water storage is from unsaturated cells only
##        ix = np.where(s > 0)  # unsaturated cells
##        uf = np.max(np.shape(ix))*self.CellArea/self.CatchmentArea
##        Rec = R*uf
##        del ix
##
##        #and remaining forms runoff
##        Roff = R - Rec
##        roff = np.zeros(np.shape(s))
##        roff[np.where(s <= 0)] = R
#
#        # update catchment water balance & comp. returnflow to surface
#        S = So + Qb - Rec
#        s = self.local_s(S)
#
#        # return flow (saturation excess) occurs from 'new' saturated cells
#        ix = np.where(s <= 0)
#        fsat = np.max(np.shape(ix))*self.CellArea / self.CatchmentArea
#        del ix
#
#        # returnflow grid
#        qf = -s
#        qf[qf < 0] = 0.0  # returnflow grid, m
#        # average returnflow per unit area
#        Qf = np.nansum(qf)*self.CellArea/self.CatchmentArea
#
#        # now all saturation excess is in Qf so update s and S
#        S = S + Qf
#        self.S = S
#
#        # check mass balance
#        # dS = (self.S - So)
#        # dF = R - Roff - Qb - Qf
#        # print ('topmbe [mm]', 1000*dS + 1000*dF)
#
#        # append results
#        if hasattr(self, 'results'):
#            self.results['S'].append(self.S)
#            self.results['Qt'].append(Qb + Qf + Roff)
#            self.results['Qb'].append(Qb)
#            self.results['Qf'].append(Qf)
#            self.results['Roff'].append(Roff)
#            self.results['fsat'].append(fsat)
#        return Qb, Qf, Roff, fsat, qf
    