import VBMicrolensing
import os
import glob

import tqdm
# import sys
# import numpy as np
# import math
# import pandas as pd
import tqdm
from joblib import Parallel, delayed
from .icgs_helpers import *
import json


# @delayed
# @wrap_non_picklable_objects
# def initialize_vbm(Tol, RelTol, turn_off_secondary_source,turn_off_secondary_lens,mass_luminosity_exponent,mass_radius_exponent,lens_mass_luminosity_exponent):
#     vbm = VBMicrolensing.VBMicrolensing()
#     vbm.Tol = Tol
#     vbm.RelTol = RelTol
#     vbm.SetMethod(VBMicrolensing.VBMicrolensing.Multipoly)
#     vbm.turn_off_secondary_source = turn_off_secondary_source
#     vbm.turn_off_secondary_lens = turn_off_secondary_lens
#     vbm.mass_luminosity_exponent = mass_luminosity_exponent
#     vbm.mass_radius_exponent = mass_radius_exponent
#     vbm.lens_mass_luminosity_exponent = lens_mass_luminosity_exponent
#     return vbm

def minimize_linear_pars(y, err, x):
    """
    Exact solution to minimize linear eq chi2
    """
    alpha1 = np.sum(x / err ** 2)
    alpha2 = np.sum(1 / err ** 2)
    alpha3 = np.sum(y / err ** 2)
    alpha4 = np.sum((x / err) ** 2)
    alpha5 = np.sum(x * y / err ** 2)
    a = (alpha5 * alpha2 - alpha1 * alpha3) / (alpha4 * alpha2 - alpha1 ** 2)
    b = (alpha4 * alpha3 - alpha5 * alpha1) / (alpha4 * alpha2 - alpha1 ** 2)
    return a, b


class ICGS:
    def __init__(self, eventname,
              satellitedir = '.', Tol = 0.01,RelTol=0.001, ncores = 1,grid_dictionary = {},overwrite=False):

        self.ncores = ncores
        self.nfil = 0
        self.lightcurves = []
        self.LCToFit = None
        self.satellites = [0]
        self.satellitedir = satellitedir
        self.eventname = eventname

        self.grid_dictionary = grid_dictionary
        self.model_type = self.grid_dictionary['modelcode']
        self.Tol = Tol
        self.RelTol = RelTol
        self.reference_chi2 = 0 # temp TODO add this!
        self.overwrite = overwrite


        # General information on models
        self.modelcodes= ['PS','PX','BS','BO','LS','LX','LO','LK','TS','TX','TO']
        self.modnumber =  self.modelcodes.index(self.model_type)
        self.npars=[4,6,7,12,7,9,12,14,10,12]
        self.logposs=[[0,1,3],
                 [1,3],
                 [0,1,6],
                 [0,1,6],
                 [0,1,4,5],
                 [0,1,4,5],
                 [0,1,4,5],
                 [0,1,4,5],
                 [0,1,4,5,7,8],
                 [0,1,4,5,7,8]]
        self.parnames = [['u0','tE','t0','rho'],
                    ['u0','tE','t0','rho','piN','piE'],
                    ['tE','FR','u01','u02','t0','t02','rho'],
                    ['tE','FR','u01','u02','t0','t02','rho','piN','piE','gamma1','gamma2','gammaz'],
                    ['s','q','u0','alpha','rho','tE','t0'],
                    ['s','q','u0','alpha','rho','tE','t0','piN','piE'],
                    ['s','q','u0','alpha','rho','tE','t0','piN','piE','gamma1','gamma2','gammaz'],
                    ['s','q','u0','alpha','rho','tE','t0','piN','piE','gamma1','gamma2','gammaz','sz_s','a_s3d'],
                    ['s','q','u0','alpha','rho','tE','t0','s2','q2','beta'],
                    ['s','q','u0','alpha','rho','tE','t0','s2','q2','beta','piN','piE']]
        self.astroparnames = ['muS_Dec','muS_RA','piS','thetaE']

        self.astrometric = False
        self.nlinpars = 2
        self.parstring = ''
        #initialize these to defaults.
        self.turn_off_secondary_source = False
        self.turn_off_secondary_lens = False

        self.mass_luminosity_exponent = 4.0
        self.mass_radius_exponent = 0.9
        self.lens_mass_luminosity_exponent = 4.0
        # run
    def run(self):
        self.readdata()
        self.readoptions()
        self.setup_grid()
        self.search_grid()
        self.create_initial_conditions()

    def readdata(self):
        if self.eventname== None:
            self.lightcurves = []
            self.telescopes = []
            self.limbdarkenings = [0]
            self.nfil = 0
            self.npoints = 0
        else:
            os.chdir(self.eventname)

            columns = ['filter', 'HJD', 'Flux', 'errFlux', 'satellite', 'Dec', 'errDec', 'RA', 'errRA']
            self.LCToFit = pd.read_csv('LCToFit.txt',skiprows=1, names = columns, sep='\s+')
            #print(self.LCToFit)
            self.npoints=self.LCToFit.shape[0]
            self.nfil = np.unique(self.LCToFit['filter']).shape[0]
            self.satellites = np.unique(self.LCToFit['filter'])
            #print(self.satellites) # test remove later TODO
            if self.LCToFit.iloc[0,-1] > 0:
                self.astrometric = True
                self.nlinpars = 4
                for i in range(len(self.npars)):
                    self.npars[i] += 4
                    self.parnames[i] += self.astroparnames

        self.lightcurves = []
        for i in range(self.nfil):
            lc = self.LCToFit[self.LCToFit['filter']==i].values
            lc_list = []
            for j in range(1,8):
                lc_list.append([lc[:,j]])
            self.lightcurves.append(lc_list)


        # self.lightcurves = [ [np.array([dl[0] for dl in d]),np.array([dl[1] for dl in d]),np.array([dl[2] for dl in d]) , np.array([dl[4] for dl in d]),np.array([dl[5] for dl in d]),np.array([dl[6] for dl in d]),np.array([dl[7] for dl in d]),d[0][3]] for d in data]
        with open('FilterToData.txt') as f:
            self.telescopes = f.readlines()
            for i in range(0,self.nfil):
                self.telescopes[i] = self.telescopes[i][0:self.telescopes[i].index('.')]
            self.telescopes = [tel.split('_')[0] for tel in self.telescopes]

        if(os.path.exists(self.eventname + '/Data/LimbDarkening.txt')):
            with open(self.eventname + '/Data/LimbDarkening.txt') as f:
                lines = f.readlines()
                self.limbdarkenings = [float(ld) for ld in lines]
        else:
            self.limbdarkenings = [0 for t in self.telescopes]

    def readoptions(self):
        if(self.eventname!=None and os.path.exists(self.eventname + '/ini/LevMar.ini')):
            with open(self.eventname + '/ini/LevMar.ini') as f:
                lines = f.readlines()
                for line in lines:
                    chunks = line.split()
                    if(chunks[0] == 'turn_off_secondary_source' and chunks[2] == 'True'):
                        self.turn_off_secondary_source = True
                    elif(chunks[0] == 'turn_off_secondary_lens' and chunks[2] == 'True'):
                        self.turn_off_secondary_lens = True
                    elif(chunks[0] == 'mass_luminosity_exponent'):
                        self.mass_luminosity_exponent = float(chunks[2])
                    elif(chunks[0] == 'mass_radius_exponent'):
                        self.mass_radius_exponent = float(chunks[2])
                    elif(chunks[0] == 'lens_mass_luminosity_exponent'):
                        self.lens_mass_luminosity_exponent = float(chunks[2])

    # @delayed
    # @wrap_non_picklable_objects
    def initialize_vbm(self):
        vbm = VBMicrolensing.VBMicrolensing()
        vbm.Tol = self.Tol
        vbm.RelTol = self.RelTol
        vbm.SetMethod(VBMicrolensing.VBMicrolensing.Multipoly)
        vbm.turn_off_secondary_source = self.turn_off_secondary_source
        vbm.turn_off_secondary_lens = self.turn_off_secondary_lens
        vbm.mass_luminosity_exponent = self.mass_luminosity_exponent
        vbm.mass_radius_exponent = self.mass_radius_exponent
        vbm.lens_mass_luminosity_exponent = self.lens_mass_luminosity_exponent
        if(self.modnumber == 1 or self.modnumber == 3 or self.modnumber > 4):
            vbm.SetObjectCoordinates(glob.glob('Data/*.coordinates')[0],self.satellitedir)
            vbm.parallaxsystem = 1

        return vbm
    def lightcurve(self, vbm, times, pars):
        results = []
        if(self.modnumber == 0):
            results = vbm.ESPLLightCurve(pars, times)
        elif(self.modnumber == 1):
            if(self.astrometric):
                results = vbm.ESPLAstroLightCurve(pars, times)
            else:
                self.results = vbm.ESPLLightCurveParallax(pars, times)
        elif(self.modnumber == 2):
            results = vbm.BinSourceExtLightCurve(pars, times)
        elif(self.modnumber == 3):
            if(self.astrometric):
                results = vbm.BinSourceAstroLightCurveXallarap(pars, times)
            else:
                results = vbm.BinSourceExtLightCurveXallarap(pars, times)
        elif(self.modnumber == 4):
            results = vbm.BinaryLightCurve(pars, times)
        elif(self.modnumber == 5):
            if(self.astrometric):
                results = vbm.BinaryAstroLightCurve(pars, times)
            else:
                results = vbm.BinaryLightCurveParallax(pars, times)
        elif(self.modnumber == 6):
            if(self.astrometric):
                results = vbm.BinaryAstroLightCurveOrbital(pars, times)
            else:
                results = vbm.BinaryLightCurveOrbital(pars, times)
        elif(self.modnumber == 7):
            if(self.astrometric):
                results = vbm.BinaryAstroLightCurveKepler(pars, times)
            else:
                results = vbm.BinaryLightCurveKepler(pars, times)
        elif(self.modnumber == 8):
            results = vbm.TripleLightCurve(pars, times)
        elif(self.modnumber == 9):
            if(self.astrometric):
                results = vbm.TripleAstroLightCurve(pars, times)
            else:
                results = vbm.TripleLightCurveParallax(pars, times)
        if not results:
            print(f'Results: {results}')
            raise ValueError("Lightcurve not calculated!")
        return results

    def get_chi2(self,flux, filter_num):
        residuals = (self.lightcurves[filter_num][1][0] - flux)/self.lightcurves[filter_num][2][0]
        chi2 = np.sum(residuals**2)
        return chi2

    def calculate(self,vbm,pars):
        lightcurve_results = []
        lightcurve_chi2 = []
        #print(pars)
        for i in range(len(self.satellites)):
            if (self.modnumber == 1 or self.modnumber == 3 or self.modnumber > 4):
                vbm.satellite = self.satellites[i]
                #print(f'sat: {vbm.satellite}')
            vbm.a1 = self.limbdarkenings[i]
            #print(f'a1: {vbm.a1}')
            times = self.lightcurves[i][0][0]
            results = self.lightcurve(vbm,times,pars)
            mul, lin = minimize_linear_pars(self.lightcurves[i][1][0],self.lightcurves[i][2][0],results[0])
            #print(mul,lin)
            flux = mul * np.array(results[0]) + lin
            lightcurve_results.append(flux)
            chi2 = self.get_chi2(flux,filter_num=i)
            lightcurve_chi2.append(chi2)
        #print('lc chi2',lightcurve_chi2)
        total_chi2 = np.sum(lightcurve_chi2)
        #print('total chi2',total_chi2)
        return total_chi2
    ### TODO (possibly) return chi2 in each filter.
    # @delayed
    # @wrap_non_picklable_objects
    def calculate_lightcurve_wrapper(self,grid_pars):
        vbm = self.initialize_vbm()
        pars = self.pars_template
        result_list = []

        for grid_index in tqdm.trange(grid_pars.shape[0],colour='red'):
            # assign grid parameters to pars
            for i in range(len(self.scan_list)):
                #check if they should be log. TODO put this outside of loop?
                if i in self.logposs[self.modnumber]:
                    pars[self.scan_list[i]] = np.log(grid_pars[grid_index, i])
                else:
                    pars[self.scan_list[i]] = grid_pars[grid_index, i]
            total_chi2  = self.calculate(vbm,pars)
            results = pars.copy()
            #print(results)
            results.append(total_chi2)
            #print(total_chi2)
            result_list.append(results)
        return result_list

    def setup_grid(self):
        if self.model_type not in self.modelcodes:
            raise ValueError("Specified model does not exist")


        if (self.modnumber == 1 or self.modnumber == 3 or self.modnumber > 4):
            #if parallax model, extract the grid dictionaries for +/- u0
            u0_plus_grid_dictionary = self.grid_dictionary['u0_plus']
            u0_minus_grid_dictionary = self.grid_dictionary['u0_minus']


        self.scan_list = [] # list of indices for parameters that will be looped over
        iterable_parameters_list = [] # list of parameter grids that get looped over
        ### if any of the grid parameter arrays are empty, fill them from predecessor model (PS -> LS ...)


        for key in self.grid_dictionary:
            if key != 'modelcode' and key != 'chi2' and len(self.grid_dictionary[key]) == 0:
                predecessor_model = FindFixedParameters(self.eventname, self.modnumber, self.nfil, self.grid_dictionary)
                self.grid_dictionary = predecessor_model.grid_dictionary
                break
        self.reference_chi2 = self.grid_dictionary['chi2']

        # now set up
        parameter_index = 0
        self.n_gridpoints=1
        for key in self.grid_dictionary:
            if key!='modelcode' and key!='chi2':
                self.n_gridpoints*=len(self.grid_dictionary[key]) #to get number of grid models.
                if len(self.grid_dictionary[key])>1:
                    self.scan_list.append(parameter_index)
                    iterable_parameters_list.append(self.grid_dictionary[key])
                parameter_index +=1
        grid_tuple = np.meshgrid(*iterable_parameters_list)
        self.grid_array = np.zeros((self.n_gridpoints,len(self.scan_list)))
        for i in range(len(self.scan_list)):
            self.grid_array[:,i] = grid_tuple[i].flatten()

    def search_grid(self):
        self.pars_template = [self.grid_dictionary[key][0] for key in self.parnames[self.modnumber]] # initialze a parameter array.
        for i in self.logposs[self.modnumber]:
            self.pars_template[i] = np.log(self.pars_template[i])

        parallel = Parallel(n_jobs=self.ncores,verbose=0)
        #split into self.ncores subarrays for parallel calculation

        # split the grid points uniformly between each process by assigning them to different arrays.
        # prevents one process from being much slower by getting all of the difficult parameters.
        indices = []
        for i in range(self.ncores):
            indices.append(np.arange(i, self.grid_array.shape[0], step=self.ncores))

        #self.grid_results = self.calculate_lightcurve_wrapper(self.grid_array)

        self.grid_results = parallel(delayed(self.calculate_lightcurve_wrapper)(self.grid_array[indices[i]]) for i in range(self.ncores))
        self.grid_results = np.concatenate(self.grid_results)
        #del self.grid_array # just to save memory.
        column_names = self.parnames[self.modnumber] +['chi2']
        self.grid_results = pd.DataFrame(np.array(self.grid_results),columns=column_names)

        # And that's the grid :)

    def create_initial_conditions(self):
        #['PS', 'PX', 'BS', 'BO', 'LS', 'LX', 'LO', 'LK', 'TS', 'TX']
        init_conds_selector = SelectInitConds(self.modnumber, self.grid_results,self.eventname,self.reference_chi2)
        self.icgs_initial_conditions = init_conds_selector.best_grid_models

        # take out of log space
        for i in self.logposs[self.modnumber]:
            self.icgs_initial_conditions.iloc[:,i] = np.exp(self.icgs_initial_conditions.iloc[:,i])


        os.mkdir(f"{self.eventname}/GridOutputs")
        self.grid_results.to_csv(f"{self.eventname}/GridOutputs/{self.model_type}ICGS.csv",index=None,)
        self.icgs_initial_conditions.to_csv(f"{self.eventname}/GridOutputs/{self.model_type}ICGS_best_models.csv",index=None)

        temppath = f'{self.eventname}/InitCond/InitCond{self.model_type}.txt.tmp'
        initpath = f'{self.eventname}/InitCond/InitCond{self.model_type}.txt'

        if self.overwrite:
            with open(temppath, 'w') as tempinit:
                with open(initpath, 'r') as oldinit:
                    npeaks_nconds = oldinit.readline()
                    npeaks = int(npeaks_nconds.split(' ')[0])
                    #old_ninit = int(npeaks_nconds.split(' ')[1])
                    npeaks_nconds = f'{npeaks} {self.icgs_initial_conditions.shape[0]}\n'
                    tempinit.write(npeaks_nconds)
                    # write peaks from old file
                    for n in range(npeaks):
                        peak = oldinit.readline()
                        tempinit.write(peak)
            # now write initial conditons to file
            self.icgs_initial_conditions.iloc[:, 0:self.npars[self.modnumber]].to_csv(temppath, sep=' ', mode='a', header=False, index=False)
        else:
            with open(temppath, 'w') as tempinit:
                with open(initpath, 'r') as oldinit:
                    npeaks_nconds = oldinit.readline()
                    npeaks = int(npeaks_nconds.split(' ')[0])
                    old_ninit = int(npeaks_nconds.split(' ')[1])
                    npeaks_nconds = f'{npeaks} {old_ninit + self.icgs_initial_conditions.shape[0]}\n'
                    tempinit.write(npeaks_nconds)
                    # write peaks from old file
                    for n in range(npeaks):
                        peak = oldinit.readline()
                        tempinit.write(peak)
                    for n in range(old_ninit):
                        line = oldinit.readline()
                        tempinit.write(line)
                        # now write initial conditons to file
            self.icgs_initial_conditions.iloc[:, 0:self.npars[self.modnumber]].to_csv(temppath, sep=' ', mode='a', header=False,
                                                                    index=False)
        os.remove(initpath)
        os.rename(temppath, initpath)



















