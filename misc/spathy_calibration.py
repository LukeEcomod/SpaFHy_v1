# -*- coding: utf-8 -*-
"""
Created on Thu Dec 08 18:13:21 2016

CODE FOR SPATHY HYDROLOGIC MODULE CALIBRATION USING SPOTPY -MODULE.

sve_calibrations: runs full model, optimizes topmodel parameters
sve_topmodel_calibration: runs spathy once, saves potential drainage to 'topmodel'
and runs topmodel and optimizes its parameters

@author: slauniai
"""
import os
import spotpy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import graphics
from spathy_sve import SpatHy
from canopygrid import CanopyGrid
# from bucketgrid import BucketGrid
from topmodel import Topmodel_Homogenous as Topmodel

from iotools import read_setup, create_catchment, read_SVE_runoff

eps = np.finfo(float).eps  # machine epsilon

# spathy_path = os.path.join(os.path.expanduser('~'),'projects','spathy') #path to spathy folder
spathy_path = os.path.join('c:', 'c:\datat\spathydata')


def sve_calibrations(setupfile, catchment, start_date, end_date, spinup_end, reps=None):
    """
    main function for parameter optimization using spotpy.
    Runs full Spathy model 'nreps' and samples parameters from distributions given in
    spotpy_setup_spathy(). Evaluates model goodness using selected objective function
    """
    cal_path = os.path.join(spathy_path, 'results', 'cal')
    # fig_path = os.path.join(spathy_path, 'figs')
    print cal_path
    """read parameter file into dicts"""
    pgen, pcpy, pbu, ptop = read_setup(setupfile)
    pgen['catchment_id'] = catchment
    pgen['start_date'] = start_date
    pgen['end_date'] = end_date

    # create spathy run for spotpy
    spot = spotpy_setup_spathy(pgen, pcpy, pbu, ptop, spinup_end)

    # create spotpy sampler & sample
    db = os.path.join(cal_path, 'dbc' + catchment)
    print db

    # run monte-carlo or linear hypercube sampler
    # sampler = spotpy.algorithms.mc(spot, dbname=db, dbformat='csv')
    sampler = spotpy.algorithms.lhs(spot, dbname=db, dbformat='csv')
    sampler.sample(reps)

    # get results and run with best parameterset
    results = sampler.getdata()
    #spotpy.analyser.plot_parametertrace(results) 
    #spotpy.analyser.plot_objectivefunction(results,spot.evaldata)
    #spotpy.analyser.plot_regression(results,spot.evaldata)
    #spotpy.analyser.plot_allmodelruns(results,spot.evaldata)
    #parname=spotpy.analyser.get_parameternames(results)
    #par=spotpy.analyser.get_parameters(results)
    data = spotpy.analyser.load_csv_results(db)
    #print data['like1'], data['parm']
    ########## create a file to store m, ko, and obj############
    par_obj = {'m': np.empty([reps-1,1]),'ko': np.empty([reps-1,1]),'obj': np.empty([reps-1,1])}
    par_obj['m'] = data['parm']
    # par_obj['ko'] = data['parko']
    par_obj['ko'] = np.ones(np.shape(data['parm']))*ptop['ko'] # fix ko
    par_obj['obj'] = data['like1']

    #print par_obj
    graphics.plot_par_objective(catchment, par_obj)
    #graphics.mesh_par_objective(catchment,par_obj)
    ############################################################
    p = list(spotpy.analyser.get_best_parameterset(results)[0])
    print p
    # p=[spot.ptop['m'], spot.ptop['ko']]
    res = spot.simulation(p)
    ####################mingfu##########
    maxNSC=spotpy.objectivefunctions.nashsutcliffe(spot.evaldata, res)
    maxNSE=spotpy.objectivefunctions.lognashsutcliffe(spot.evaldata, res)
    print maxNSC, maxNSE
    
    #######################################
    # plot and save figure
    plt.figure(num=catchment)
    plt.clf()
    plt.tick_params(labelsize=10)
    
    plt.plot(spot.evaldates, spot.evaldata, 'k-', spot.evaldates, res, 'r-')
    Prec = spot.spathy.FORC['Prec'].ix[spot.evaldates].values.tolist()
    P = np.nansum(Prec)*spot.spathy.dt
    R = np.nansum(res)
    Qm = sum(spot.evaldata)

    #txt = 'P=' + str(P) + '  Qm=' + str(Qm) + ' R=' + str(R) + '  Qm/P=' + str(Qm/P) + '  R/P=' + str(R/P) + '  NSC=' + str(maxNSC)+', '+ str(maxNSE) 
    txt = 'Q/P=%.2f Q$_{mod}$/P=%.2f RI=%.2f NSE=%.2f'  %(Qm/P, R/P, max(par_obj['obj']), maxNSE)
    plt.xlabel(txt, fontsize=8)
    head = 'C' + str(catchment) + ': m=%.3f' %(p[0]) # + ', k_o=' + str(p[1]) 
    # +',gs='+str(p[3]) # +',w=' +str(p[3]) +',ws='+str(p[4]) 
    plt.title(head, fontsize=8)
    plt.ylabel('R mm/d')
#    manager = plt.get_current_fig_manager()
#    manager.window.showMaximized()  
    plt.savefig(db + '.png')
    plt.close('all')
#        posterior=spotpy.analyser.get_posterior(results, threshold=0.9);
#        spotpy.analyser.plot_parameterInteraction(posterior)
#        plt.title(catchment); plt.savefig(db+'_post_0.9.png')
#        plt.close('all')
    return spot, res  # spot.evaldates, spot.evaldata#, res, Prec


class spotpy_setup_spathy(object):
    """
    spotpy parameter optimization module setup script for SpatHy
    """
    def __init__(self, pgen, pcpy, pbu, ptop, spinup_end, p=None, flatten=True):
        """
        spinup_end - 'yyyy-mm-dd'; end of spinup-period
        p - list of variables & bounds
            p=[ ['varname', low,high,stepsize,optquess,minbound,maxbound], ...]
        catchment_id - to override catchment_id in setupfile
        """
        # from spathy import SpatHy
        self.pgen, self.pcpy, self.pbu, self.ptop = pgen, pcpy, pbu, ptop

        gisdata = create_catchment(self.pgen['catchment_id'],
                                        fpath=self.pgen['gis_folder'],
                                        plotgrids=False, plotdistr=False)

        self.gisdata = gisdata

        """ 
        To speed up calculations, omit cells outside catchment and flatten arrays into 1d
        NOTE: spatial location is lost and current netcdf format does not work with this.
        Re-mapping to original grid is easy, though.
        """
        if flatten:                                
            cmask = gisdata['cmask'].copy()
            ix = np.where(np.isfinite(cmask))
            keys = {'cmask', 'LAI_conif', 'LAI_decid', 'hc', 'cf', 'soil', 'flowacc', 'slope'}
            flatgis = {key: gisdata[key][ix] for key in keys}
            
            flatgis['cellsize'] = gisdata['cellsize']
            self.gisdata = flatgis
            del flatgis
            
        self.cmask = self.gisdata['cmask']
        
        """ create spathy model object """
        self.spathy = SpatHy(self.pgen, self.pcpy, self.pbu, self.ptop, self.gisdata, flatten=True) 
        self.Nsteps = self.spathy.Nsteps
        print('INIT M',self.spathy.top.M, 'init to', self.spathy.top.To)
        
        """ **** load evaluation runoff data **** """
        
        roff = read_SVE_runoff(self.pgen['catchment_id'], self.pgen['start_date'],
                               self.pgen['end_date'], pgen['runoff_file'])
        roff = roff.dropna()
        roff[roff < eps] = eps
        # plt.plot(roff)

        self.evaldates = roff[roff.index > spinup_end].index
        # print self.evaldates[0]; print self.evaldates[-1]
        self.evaldata = roff[roff.index > spinup_end].values.tolist()
        del roff
        
        """
        # print len(self.evaldates); print len(self.evaldata)
        ##################################################################################################
        #sourcefile = 'D:\LUKE_work\SVE_catchment\measured\\'+self.pgen['catchment_id']+'\\tmp\\roff.dat'
        sourcefile = pgen['runoff_file']
        print sourcefile
        roff = pd.read_csv(sourcefile)
        utctime=pd.to_datetime(roff['utctime'],format='%Y-%m-%d')
        roff.index=utctime
        roff=roff.drop('utctime',1)
        roff = roff[(roff.index >= self.pgen['start_date']) & (roff.index <= self.pgen['end_date'])]
        roff.columns = ['Qm']
        self.evaldates = roff[roff.index > spinup_end].index
        self.evaldata = roff['Qm'][roff.index > spinup_end].values.tolist()
        ##################################################################################################

        """
        """
        create parameter distributions: MODIFY HERE ranges and parameters to be
        calibrated
        """
        # spotpy.parameter.distric('paramname',low,high,stepsize,optquess,minbound,maxbound) 
        self.params = [spotpy.parameter.Uniform('m', 0.002, 0.05, 0.005, 0.01)]
                       # spotpy.parameter.Uniform('ko', 1e-4, 1e-3, 0.01, 5e-4)] # fix ko
                       # spotpy.parameter.Uniform('gsref',1.8e-3,2.7e-3,1e-4,2.3e-3),
                       # spotpy.parameter.Uniform('f',0.6,0.75,0.01,0.68),
                       # spotpy.parameter.Uniform('wmax',1.5,2.5,0.1,2.0),
                       # spotpy.parameter.Uniform('wmaxsnow',3.75,6.25,0.1,5.0),
                       # spotpy.parameter.Uniform('kmelt',2e-5,4e-5,1e-6,2.8934e-05)

        self.simcount=0

    def parameters(self):
        # samples parameters from distribution
        return spotpy.parameter.generate(self.params)

    def simulation(self, vector):
        """
        run simulation, vary parameters. vector is in order:
        vector[0]=m, [1]=ko
        """

        """ adjust parameter values and update spa components """
        self.ptop['m'] = vector[0]
        # self.ptop['ko'] = vector[1] #fix ko

#        """--- CanopyGrid and BucketGrid ---: not calibrated now """
#
#        if self.pgen['spatial_cpy'] == 'no': #spatially constant stand properties
#            self.spathy.cpy=CanopyGrid(self.pcpy,cmask=self.GisData['cmask']); 
#        else:
#            self.spathy.cpy=CanopyGrid(self.pcpy,lai=self.GisData['LAI_pine']+self.GisData['LAI_spruce']+self.GisData['LAI_decid'], ba=self.GisData['ba'],\
#                hc=self.GisData['hc']); 
#            
#        if self.pgen['spatial_soil'] == 'no': #spatially constant soil properties
#            self.spathy.bu=BucketGrid(self.pbu, cmask=self.GisData['cmask'])
#        else:
#            self.spathy.bu=BucketGrid(self.pbu, soiltypefile=self.pgen['soil_paramfile'],sdata=self.GisData['soil'])

        """ --- initialize homogenous topmodel --- """
        self.spathy.top = Topmodel(self.ptop, self.gisdata['cellsize']**2, self.cmask,
                            self.gisdata['flowacc'], slope=self.gisdata['slope'], outputs=False)

        """ ---- run Spathy and return modeled streamflow (mm/d)---"""
        self.simcount += 1
        # print 'Spathy run : ' +str(self.simcount)

        res = self.spathy._run(0, self.Nsteps, calibr=True)
        simulations = res['Qm'].ix[self.evaldates].values.tolist()

        # print('simlen',len(simulations))
        return simulations  # returns Qmod

    def evaluation(self, evaldates=False):
        if evaldates:
            return self.evaldates
        else:
            # print('evallen',len(self.evaldata))
            return self.evaldata

    def objectivefunction(self, simulation, evaluation):
        
        evaluation = evaluation + eps
        simulation = simulation + eps
        
        def modified_nashsutchliffe(evaluation, simulation):
            """ modified nash-suchliffe -criteria"""
            if len(evaluation) == len(simulation):
                s, e = np.array(simulation), np.array(evaluation)
                # s, e = simulation, evaluation
                mean_observed = np.mean(e)
                # compute numerator and denominator
                numerator = sum(abs(e - s))
                denominator = sum(abs(e - mean_observed))
                # compute coefficient
                obj = 1 - (numerator / (denominator + eps))
                print obj
                return 1 - (numerator / (denominator + eps))

            else:
                print("Error: evaluation and simulation lists does not have the same length.")
                return np.nan

        def modified_agreementindex(evaluation, simulation):
            """ modified willmot's agreement index"""
            if len(evaluation) == len(simulation):
                s, e = np.array(simulation), np.array(evaluation)
                # s, e = simulation, evaluation
                
                # compute numerator and denominator
                numerator = sum(abs(e - s))
                denominator = sum(abs(s - np.mean(e)) + abs(e - np.mean(e)))
                obj = 1 - (numerator / (denominator + eps))
                print obj    
                return 1 - (numerator / (denominator + eps))

            else:
                print("Error: evaluation and simulation lists does not have the same length.")
                return np.nan

        def relative_agreementindex(evaluation, simulation):
            """ relative willmot's agreement index"""
            if len(evaluation) == len(simulation):
                s, e = np.array(simulation), np.array(evaluation)
                # s, e = simulation, evaluation
                
                # compute numerator and denominator
                numerator = sum(((e - s) / e)**2)
                denominator = sum(((abs(s - np.mean(e)) + abs(e - np.mean(e))) / np.mean(e))**2)
                obj = 1 - (numerator / (denominator + eps))
                print obj, numerator, denominator
                return 1 - (numerator / (denominator + eps))

            else:
                print("Error: evaluation and simulation lists does not have the same length.")
                return np.nan
            

        # # normalize to zero mean, unit variance to remove bias
#        evaluation=evaluation - np.mean(evaluation)
#        evaluation=list(evaluation / np.std(evaluation))
#        simulation=simulation-np.mean(simulation)
#        simulation=list(simulation / np.std(simulation))

        # select objective function to be used
        # objectivefunction = _MNS(evaluation, simulation)
        # objectivefunction = -spotpy.objectivefunctions.rmse(evaluation,simulation)
        # print('simu', len(simulation), 'eva', len(evaluation))
        #objectivefunction = spotpy.objectivefunctions.nashsutcliffe(evaluation, simulation)

        #objectivefunction = spotpy.objectivefunctions.lognashsutcliffe(evaluation, simulation)
        #objectivefunction = spotpy.objectivefunctions.agreementindex(evaluation, simulation)
        #objectivefunction = relative_agreementindex(evaluation, simulation)
        objectivefunction = modified_agreementindex(evaluation, simulation)

        return objectivefunction
    
    """ ********** calibration script for stand-alone topmodel ************* """
def sve_topmodel_calibration(setupfile, catchment, start_date, end_date, spinup_end, reps=None):

    """
    main function for Topmodel parameter optimization using spotpy.
    Runs full Spathy model once, saves drainage to list and repeats topmodel 
    run 'nreps' times and samples parameters from distributions given in
    spotpy_setup_topmodely(). Evaluates model goodness using selected objective function.

    saves outputs to folder spathy_path/results, csv -format
    """
    from spathy import initialize_netCDF
    from netCDF4 import Dataset
    
    cal_path = os.path.join(spathy_path, 'results')

    """read parameter file into dicts"""
    pgen, pcpy, pbu, ptop = read_setup(setupfile)
    pgen['catchment_id'] = catchment
    pgen['start_date'] = start_date
    pgen['end_date'] = end_date

    gisdata = create_catchment(pgen['catchment_id'], fpath=pgen['gis_folder'],
                               plotgrids=False, plotdistr=False)
              
    """ greate SpatHy object, run and save results to netcdf"""
    spa = SpatHy(pgen, pcpy, pbu, ptop, gisdata, bu_outputs=False, flatten=False)
    Nsteps = spa.Nsteps

    # evaluation data == measured runoff
    roff = spa.Qmeas
    roff = roff.dropna()
    roff = roff[roff.index > spinup_end]
    forcingdates = spa.FORC.index

    """ create temporary netcdf storage """
    ncf, ofi = initialize_netCDF(spa.id, spa.GisData, spa.FORC, fpath=spa.pgen['output_folder'], fname='tmp.nc')
                                      
    """ ----- Run spathy, get recharge to topmodel from bu.results['Drain']  ----- """
    spa._run(0, Nsteps, calibr=False, ncf=ncf)
    
    # select only elements inside catchment    
    cmask = gisdata['cmask'].copy()
    ix = np.isfinite(cmask)
    cmask = cmask[ix].copy()
    flowacc = gisdata['flowacc'][ix].copy()
    slope = gisdata['slope'][ix].copy()
    
#    recharge = spa.bu.results['Drain'];
#    del spa.bu.results
    """ get drainage data from netcdf and remove file"""
    tmp = Dataset(ofi, 'r')
    aa = tmp['bu']['Drain'][:].copy()

    recharge = aa[:, ix].copy()
    del aa
    tmp.close()
    os.remove(ofi)

    # create spotpy setup
    spot = spotpy_setup_topmodel(ptop, spa.GisData['cellsize'], cmask, flowacc,
                                 slope, recharge, roff, forcingdates, spinup_end)

    # create spotpy sampler & sample
    db = os.path.join(cal_path, 'top_dbc' + catchment)

    # run monte-carlo or linear hypercube sampler
    sampler = spotpy.algorithms.mc(spot, dbname=db, dbformat='csv')
    sampler = spotpy.algorithms.lhs(spot, dbname=db, dbformat='csv')
    sampler.sample(reps)

    # get results and run with best parameterset
    results = sampler.getdata()
    p = list(spotpy.analyser.get_best_parameterset(results)[0])
    print p
    # p=[spot.ptop['m'], spot.ptop['ko']]
    res = spot.simulation(p)

    # plot and save figure
    plt.figure(num=catchment)
    plt.clf()
    plt.tick_params(labelsize=10)

    plt.plot(spot.evaldates, spot.evaldata, 'k-', spot.evaldates, res, 'r-')
    Prec = spa.FORC['Prec'].loc[spot.evaldates].values.tolist()
    P = np.nansum(Prec)*spa.dt
    R = np.nansum(res)
    Qm = sum(spot.evaldata)

    txt = 'P=' + str(P) + '  Qm=' + str(Qm) + ' R=' + str(R) + '  Qm/P=' + str(Qm/P) + '  R/P=' + str(R/P)
    plt.xlabel(txt, fontsize=10)
    head = 'C' + str(catchment) + ': m=' + str(p[0]) + ', k_o=' + str(p[1]) 
    # +',gs='+str(p[3]) # +',w=' +str(p[3]) +',ws='+str(p[4]) 
    plt.title(head, fontsize=10)
    plt.ylabel('R mm/d')
    plt.savefig(db+'.png')
    
#        posterior=spotpy.analyser.get_posterior(results, threshold=0.9);
#        spotpy.analyser.plot_parameterInteraction(posterior)
#        plt.title(catchment); plt.savefig(db+'_post_0.9.png')
#        plt.close('all')
    return spot, res  # spot.evaldates, spot.evaldata#, res, Prec


class spotpy_setup_topmodel(object):
    """
    spotpy parameter optimization module setup script for topmodel.
    NOTE: test calibrating both m and ko, then fix ko and re-do calibration only for m    
    """
    def __init__(self,ptop, cellsize, cmask, flowacc, slope, recharge, roff, forcingdates, spinup_end):
        """
        ptop - topmodel parameter dict
        cellsize - cellsize (float)
        cmask - catchment mask array
        flowacc - flow accumulation array
        slope - slope array
        recharge - list of recharge arrays (drainage from spathy.bucketgrid)
        roff - pd dataframe of measured streamflow (evaluation data)
        spinup_end - 'yyyy-mm-dd'; end of spinup-period
        """
        self.ptop = ptop
        self.cellsize = cellsize
        self.cmask = cmask
        self.flowacc = flowacc
        self.slope = slope
        self.recharge = recharge
        self.forcingdates = forcingdates
        self.evaldates = roff.index.values
        self.evaldata = roff.values.tolist()
        
        """ --- initialize homogenous topmodel --- """
        self.top = Topmodel(self.ptop, self.cellsize**2, self.cmask,
                            self.flowacc, self.slope, outputs=False)
        #print('INIT M',self.top.M, 'init to', self.top.To)
        self.Nsteps = len(self.recharge)

        """
        create parameter distributions: MODIFY HERE ranges and parameters to be
        calibrated
        """
        # spotpy.parameter.distric('paramname',low,high,stepsize,optquess,minbound,maxbound)
        self.params = [spotpy.parameter.Uniform('m', 1e-4, 0.1, 0.001, 0.01),
                       spotpy.parameter.Uniform('ko', 0.001, 0.05, 0.01, 0.01)]

        self.simcount = 0

    def parameters(self):
        # samples parameters from distribution
        return spotpy.parameter.generate(self.params)

    def simulation(self, vector):
        """
        run simulation, vary parameters. vector is in order:
        vector[0]=m, [1]=ko
        """

        """ adjust parameter values and update spa components """
        self.ptop['m'] = vector[0]
        self.ptop['ko'] = vector[1]

        """ --- initialize homogenous topmodel --- """
        self.top = Topmodel(self.ptop, self.cellsize**2, self.cmask,
                            self.flowacc, self.slope, outputs=False)

        """ ---- run Topmodel and return modeled streamflow (mm/d)---"""
        self.simcount += 1
        res = {'Qmod': [None]*self.Nsteps}

        for k in range(0, self.Nsteps):
            # compute drainage to topmodel bucket; this is minimum of local 'potential'
            # drainage (from root zone) and local storage deficit (s) 
            # note: potential drainage input from full spathy-simulation is affected by topmodel
            # parameters (slightly); remains to be seen which is the impact. Maybe test in two rounds;
            # in second phase set m and ko in spathy.ini according to parameters from 1st round.
                    
            potdrain = self.recharge[k]
            s = self.top.local_s(self.top.S)
            s[s < 0] = 0.
            drain = np.minimum(potdrain, s)           
            infi_ex = np.nansum((potdrain - drain))*self.top.CellArea / self.top.CatchmentArea

            # catchment average ground water recharge [m per unit area]
            RR = np.nansum(drain)*self.top.CellArea / self.top.CatchmentArea        
            # RR = 0.0
            qb, qr, roff, satf, _ = self.top.run_timestep(RR)

            res['Qmod'][k] = 1e3*(qb + qr + roff + infi_ex)  # mm/timestep

        res = pd.DataFrame.from_dict(res)
        res.index = self.forcingdates
        simulations = res['Qmod'].loc[self.evaldates].values.tolist()

        return simulations  # returns Qmod

    def evaluation(self, evaldates=False):
        if evaldates:
            return self.evaldates
        else:
            # print('evallen',len(self.evaldata))
            return self.evaldata

    def objectivefunction(self, simulation, evaluation):
        # # normalize to zero mean, unit variance to remove bias
        # evaluation=evaluation - np.mean(evaluation)
        # evaluation=list(evaluation / np.std(evaluation))
        # simulation=simulation-np.mean(simulation)
        # simulation=list(simulation / np.std(simulation))

        # select objective function to be used
        # objectivefunction = _MNS(evaluation, simulation)
        # objectivefunction = -spotpy.objectivefunctions.rmse(evaluation,simulation)
        objectivefunction = spotpy.objectivefunctions.nashsutcliffe(evaluation, simulation)
        # objectivefunction = spotpy.objectivefunctions.lognashsutcliff(evaluation,simulation)

        def _MNS(evaluation, simulation):
            """ modified nash-suchliffe -criteria"""
            if len(evaluation) == len(simulation):
                s, e = np.array(simulation), np.array(evaluation)
                # s, e = simulation, evaluation
                mean_observed = np.mean(e)
                # compute numerator and denominator
                numerator = sum(abs(e - s))
                denominator = sum(abs(e - mean_observed))
                # compute coefficient
                return 1 - (numerator/denominator)

            else:
                print("Error: evaluation and simulation lists does not have the same length.")
                return np.nan

        return objectivefunction

""" ********** CALIBRATION SCRIPT FOR CanopyGrid *********** """


class CanopyGrid_cal(object):
    """
    initializes CanopyGrid for calibration against Hyytiala -data.
    Reads forcing data and evaluation data, returns class instance.
    """
    def __init__(self, paramfile, evalvar='ET'):
        # NOTE: if evalvars is different, check that evaldata and evaldates
        # order is correct here and in run_cal
        # evalvar='ET', 'Trfall', 'SWE'
        from iotools import read_HydeDaily
        # path to spathy folder
        spathy_path = os.path.join(os.path.expanduser('~'), 'projects', 'spathy')
        fname = os.path.join(spathy_path, 'work/cal', 'HydeDaily2000-2010.txt')
        # print fname
        """read parameter file into pp dict"""
        pgen, pcpy, pbu, ptop = read_setup(paramfile)
        self.params = pcpy
        self.dt=pgen['dt']

        """ read forcing data and evaluation data """
        dat, FORC = read_HydeDaily(fname)
        FORC['Prec'] = FORC['Prec'] / self.dt  # mms-1
        self.FORC = FORC
        del FORC
        
        self.Evalvar = evalvar
        if evalvar is 'Trfall':
            tmp = dat[dat['doy'].between(120, 273, inclusive=True)]['Trfall'].dropna()
            self.evaldata = tmp.values.tolist()
            self.evaldates = tmp.index
        elif evalvar is 'ET':
            tmp = dat[dat['Prec'].between(0,0.1, inclusive=True)]['ET'].dropna()
            self.evaldata = tmp.values.tolist()
            self.evaldates = tmp.index           
        elif evalvar is 'SWE':
            tmp = dat['SWE'].dropna()
            self.evaldata = tmp.values.tolist()
            self.evaldates = tmp.index

        self.Nsteps = len(self.FORC)
        self.simcount = 1
        
    
    def run_cal(self, paramset=None, gsref=None, f=None, kmelt=None, wmax=None, wmaxsnow=None, rw=None, rwmin=None, resu=False):
        """runs CanopyGrid for varying parameter values"""
        """
        initialize CanopyGrid: identical to pure Spathy except that now self.params are used
        """
        #adjust parameter dicts for new parameter guesses
        if paramset:
            if self.Evalvar is 'ET':
                self.params['gsref'] = paramset[0]
                self.params['f'] = paramset[1]
                self.params['rw'] = paramset[2]
                self.params['rwmin'] = paramset[3]
            if self.Evalvar is 'Trfall':
                self.params['Wmax'] = paramset[0]
            if self.Evalvar is 'SWE':
                self.params['Wmaxsnow'] = paramset[0]
                self.params['kmelt'] = paramset[1]
                # self.params['kmt'] = paramset[1]
                # self.params['kmr'] = paramset[2]
        else:
            if gsref: self.params['gsref'] = gsref
            if f: self.params['f'] = f
            if kmelt: self.params['kmelt'] = kmelt
            if wmax: self.params['wmax'] = wmax
            if wmaxsnow: self.params['wmaxsnow'] = wmaxsnow        
            if rw: self.params['rw'] = rw
            if rwmin: self.params['rwmin'] = rwmin

        # initialize canopygrid; make needs at least 2 cells;
        # later values extracted from [0]
        cpy = CanopyGrid(self.params, cmask=np.ones(1))
        # print cpy.X
        # print 'run : ' +str(self.simcount)
        self.simcount += 1

        res = {'ET': [None]*self.Nsteps, 'Evap': [None]*self.Nsteps,
               'Trfall': [None]*self.Nsteps, 'SWE':[None]*self.Nsteps,
               'fQ': [None]*self.Nsteps,'fD': [None]*self.Nsteps,
                'fRew': [None]*self.Nsteps, 'fPheno': [None]*self.Nsteps,
                'X': [None]*self.Nsteps, 'Mbe': [None]*self.Nsteps}        

        for k in range(0, self.Nsteps):
            #print 'k='+str(k)
        
            doy = self.FORC['doy'].iloc[k]; ta = self.FORC['Ta'].iloc[k]
            vpd = self.FORC['VPD'].iloc[k]; rg = self.FORC['Rg'].iloc[k];
            par = self.FORC['Par'].iloc[k]; prec = self.FORC['Prec'].iloc[k]
            rew = self.FORC['Rew'].iloc[k]
            
            # run CanopyGrid
            potinf, trfall, interc, evap,et, mbe = cpy.run_CanopyGrid(doy, self.dt, ta, prec, rg, par, vpd, Rew=rew, P=101300.0)

            res['ET'][k] = et[0]; res['Evap'][k] = evap[0]; 
            res['Trfall'][k] = trfall[0]; res['SWE'][k] = cpy.SWE[0]
            res['X'][k] = cpy.X; res['Mbe'][k] = mbe
            # res['fQ'][k]=fQ; res['fD'][k]=fD; res['fRew'][k]=fRew; res['fPheno'][k]=fPheno; 

        # return simulated, take same 'samples' as in evaldata
        res = pd.DataFrame.from_dict(res)
        res.index = self.FORC.index

#        simudata=[res['ET'].ix[self.evaldates].values.tolist(), res['TF'].ix[self.evaldates[1]].values.tolist(), 
#                  res['SWE'].ix[self.evaldates[2]].values.tolist() ]
        if resu:
            return res, cpy
        else:  # check what to return
            simudata = res[self.Evalvar].ix[self.evaldates].values.tolist()
            return simudata

class spotpy_setup_CanopyGrid(object):
    """ 
    spotpy parameter optimization module setup script for Spathy CanopyGrid
    """
    def __init__(self, paramfile, datapath, evalvar):
        """
        paramfile - path to Spathy paramfile 
        spinup_end - 'yyyy-mm-dd'
        """
        self.model = CanopyGrid_cal(paramfile, datapath, evalvar)
        # print len(self.model.evaldata)
        if self.model.Evalvar == 'ET':
            self.params = [spotpy.parameter.Uniform('gsref', 1.5e-3, 3e-3, 1e-4, 2e-3),  # ms-1
                           spotpy.parameter.Uniform('f', 0.0, 0.7, 1e-2, 0.6),  #-
                           spotpy.parameter.Uniform('rw', 0.1, 0.25, 0.001, 0.25),  # rew decay threshold
                           spotpy.parameter.Uniform('rwmin', 0.3, 0.6, 0.001, 0.35)]  # minimum transpiration ratio [-]
                                   #spotpy.parameter.distrib('paramname',low,high,stepsize,optquess,minbound,maxbound)
        if self.model.Evalvar == 'Trfall':
            self.params = [spotpy.parameter.Uniform('wmax', 1.0, 5.0, 0.01, 2.0)] #mm/unit of LAI]
        if self.model.Evalvar == 'SWE':
            self.params=[spotpy.parameter.Uniform('wmaxsnow', 2.0, 10.0, 0.1, 5.0),  # mm /unit of LAI
                         spotpy.parameter.Uniform('kmelt', 2e-5, 4e-5, 1e-6, 2.8934e-05) ]  # mm s-1            
#                         spotpy.parameter.Uniform('kmt', 1.0e-5, 3.5e-5, 1e-6, 2.3e-05),  # mm s-1degC-1                        
#                         spotpy.parameter.Uniform('kmr', 1e-7, 2e-7, 1e-8, 1.4e-07)]  # mmd-1

    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def simulation(self, vector):
        # run simulation, vary parameters
        simulations = self.model.run_cal(paramset=vector)

        return simulations

    def evaluation(self, evaldates=False):
        if evaldates:
            return self.model.evaldates
        else:
            # print('evallen',len(self.spathy.evaldata))
            return self.model.evaldata

    def objectivefunction(self, simulation, evaluation):
        obj = -spotpy.objectivefunctions.rmse(evaluation, simulation)
        # objectivefunction= -spotpy.objectivefunctions.nashsutcliff(evaluation,simulation)
        return obj


""" *********** joint calibration for carboage and hyytiälä -sites ******"""


class CanopyGrid_multisitecal(object):
    """
    initializes CanopyGrid for calibration against Hyytiala & carboage -dry-canopy ET data. 
    Reads forcing data and evaluation data, returns class instance.
    """
    
    def __init__(self, paramfile, datapath, lai=[4.5, 1.8, 0.6], hc=[15.0, 1.7, 0.4], cf=[0.6, 0.45, 0.25]):

        from iotools import read_HydeDaily, read_CageDaily
        """read parameter file into pp dict"""
        pgen, pcpy, pbu, ptop = read_setup(paramfile)
        self.params = pcpy        
        self.dt=pgen['dt']
        
        self.params = pcpy
        self.params['lai'] = lai
        self.params['hc'] = hc
        self.params['cf'] = cf
        
        """ read forcing data and evaluation data """        
        dat,FORC = read_HydeDaily(datapath + 'HydeDaily2000-2010.txt')
        FORC['Prec'] = FORC['Prec'] / self.dt  # mms-1
        self.FORC = FORC
        del FORC

        # evaluation data
        evaldata = []
        evaldates = []

        # evaldata[0] = hyde
        tmp = dat[dat['Prec'].between(0, 0.1, inclusive=True)]['ET'].dropna()
        evaldata.append(tmp.values.tolist())
        evaldates.append(tmp.index)
        del tmp

        # read cage data: evaldata[1]=cage12yr, evaldata[2]=cage4yr

        dat1, dat2, _, _ = read_CageDaily(datapath)
        tmp = dat2[dat2['Prec'].between(0, 0.1, inclusive=True)]['ET'].dropna()
        evaldata.append(tmp.values.tolist())
        evaldates.append(tmp.index)
        del tmp

        tmp = dat1[dat1['Prec'].between(0, 0.1, inclusive=True)]['ET'].dropna()
        evaldata.append(tmp.values.tolist())
        evaldates.append(tmp.index)
        del tmp

        self.evaldates = evaldates
        self.evaldata = evaldata
        del evaldates, evaldata

        self.Nsteps = len(self.FORC)
        self.simcount = 1

    def run_cal(self, gsref=None, f=None, kmelt=None, wmax=None, wmaxsnow=None, rw=None, rwmin=None):
        """runs CanopyGrid for varying parameter values"""
        
        """
        initialize CanopyGrid: identical to pure Spathy except that now self.params are used
        """
        #adjust parameter dicts for new parameter guesses

        if gsref: self.params['gsref'] = gsref
        if f: self.params['f'] = f
#        if kmelt: self.params['kmelt'] = kmelt
#        if wmax: self.params['wmax'] = wmax
#        if wmaxsnow: self.params['wmaxsnow'] = wmaxsnow        
        if rw: self.params['rw'] = rw
        if rwmin: self.params['rwmin'] = rwmin

        # initialize canopygri
        cpy = CanopyGrid(self.params, cmask=np.ones(3))
        # print cpy.X
        # print 'run : ' +str(self.simcount)
        self.simcount += 1
             
        res = {'ET0': [None]*self.Nsteps, 'ET1': [None]*self.Nsteps, 'ET2': [None]*self.Nsteps}
        for k in range(0, self.Nsteps):
            # print 'k='+str(k)
            # forcing
            doy=self.FORC['doy'].iloc[k]; ta=self.FORC['Ta'].iloc[k];vpd=self.FORC['VPD'].iloc[k]; rg=self.FORC['Rg'].iloc[k];
            par=self.FORC['Par'].iloc[k]; prec=self.FORC['Prec'].iloc[k]; rew=self.FORC['Rew'].iloc[k]
            rew2=np.ones(3); rew2[0]=rew
            # run CanopyGrid
            potinf, trfall, interc, evap,et, mbe = cpy.run_CanopyGrid(doy, self.dt, ta, prec, rg, par, vpd, Rew=rew2, P=101300.0)

            res['ET0'][k] = et[0]
            res['ET1'][k] = et[1]
            res['ET2'][k] = et[2]
            
        # return simulated, take same 'samples' as in evaldata
        res = pd.DataFrame.from_dict(res)
        res.index = self.FORC.index
        simudata = [res['ET0'].ix[self.evaldates[0]].values.tolist(),
                    res['ET1'].ix[self.evaldates[1]].values.tolist(),
                    res['ET2'].ix[self.evaldates[2]].values.tolist()]
        return simudata #,res

class spotpy_setup_multisite(object):
    """ 
    spotpy parameter optimization module setup script for Spathy CanopyGrid
    """
    def __init__(self, paramfile, datapath):
        """ 
        paramfile - path to Spathy paramfile 
        spinup_end - 'yyyy-mm-dd'
        """
        self.model = CanopyGrid_multisitecal(paramfile, datapath, 
                lai=[4.0, 1.8, 0.6], hc=[16.0, 1.7,0.4], cf=[0.6, 0.45, 0.25])  # create Spathy_cal object
        print len(self.model.evaldata)
        self.params = [spotpy.parameter.Uniform('gsref', 1e-3, 5e-3, 1e-4, 2.1e-3),  # ms-1
                       spotpy.parameter.Uniform('f',0.0,0.7,1e-2,0.6),  # -
                       spotpy.parameter.Uniform('rw',0.0,0.3,0.001,0.11),  # rew decay threshold
                       spotpy.parameter.Uniform('rwmin',0.0,0.5,0.001,0.34)  # minimum transpiration ratio [-]
                      ]
#                       spotpy.parameter.Uniform('kmelt',2e-5,4e-5,1e-6,2.8934e-05), #mmd-1
#                       spotpy.parameter.Uniform('wmax',2.0,5.0,0.1,2.0), #mm/unit of LAI
#                       spotpy.parameter.Uniform('wmaxsnow',2.0,10.0,0.1,5.0), #mm /unit of LAI
                       #spotpy.parameter.distric('paramname',low,high,stepsize,optquess,minbound,maxbound)

    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def simulation(self, vector):
        # run simulation
        simulations = self.model.run_cal(gsref=vector[0], f=vector[1], rw=vector[2], rwmin=vector[3])
        return simulations #returns Qmod

    def evaluation(self, evaldates=False):
        if evaldates:
            return self.model.evaldates
        else:
            # print('evallen',len(self.spathy.evaldata))
            return self.model.evaldata

    def objectivefunction(self, simulation, evaluation):
        obj1 = -spotpy.objectivefunctions.rmse(evaluation[0], simulation[0])
        obj2 = -spotpy.objectivefunctions.rmse(evaluation[1], simulation[1])
        obj3 = -spotpy.objectivefunctions.rmse(evaluation[2], simulation[2])
        # weighting: 0.5 Hyde, 0.25 each cage site
        obj = 0.5*obj1 + 0.25*obj2 + 0.25*obj3
        return obj
