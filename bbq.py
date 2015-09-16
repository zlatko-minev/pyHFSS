# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 11:46:56 2015

@author: Zaki
"""

from hfss import *
import numpy as np
import h5py
import time
import os
from scipy.constants import *

root_dir = 'Y:/Data/PumpingCats/HFSS/Analyzed'
fluxQ = hbar / (2*e)  

def fact(n):
    if n <= 1:
        return 1
    return n * fact(n-1)

def nck(n, k):
    return fact(n)/(fact(k)*fact(n-k))

class Bbq(object):
    """ 
    This class defines a BBQ object which calculates and saves
    Hamiltonian parameters from an HFSS simulation
    """
    
    def __init__(self, project, design):
        self.project = project
        self.design = design
        self.setup = design.get_setup()
        self.nmodes = int(self.setup.n_modes)
        self.listvariations = design._solutions.ListVariations(str(self.setup.solution_name))
        self.nominalvariation = design.get_nominal_variation()        
        self.nvariations = np.size(self.listvariations)
        self.solutions = self.setup.get_solutions()
        
        self.setup_data()
        
    def setup_data(self):
        data_dir = root_dir + '/' + self.project.name 
        if not os.path.isdir(data_dir):
            os.makedirs(data_dir)
        self.data_dir = data_dir
        self.data_filename = self.data_dir + '/' + self.design.name + '_' + time.strftime('%Y%m%d_%H%M%S', time.localtime()) + '.hdf5'
        self.h5file = h5py.File(self.data_filename, 'w')   

    def get_p_j(self, modes=None, variation=None):
        lv = self.get_lv(variation)            
        if modes is None:
            modes = range(self.nmodes)

        pjs = {}
        for m in modes:
            self.solutions.set_mode(m+1, 0)
            fields = self.setup.get_fields()
            P_J = fields.P_J
            pjs['pj_'+str(m)]=P_J.evaluate(lv=lv)
        self.pjs = pjs
        return pjs
    
    def get_freqs_bare(self):
        freqs_bare = {}
        for m in range(self.nmodes):
            freqs_bare['freq_bare_'+str(m)] = 1e9*self.solutions.eigenmodes()[m]
        self.freqs_bare = freqs_bare
        return freqs_bare
        
    def get_lv(self, variation):
        if variation is None:
            lv = self.nominalvariation
            lv = self.parse_listvariations(lv)
        else:
            lv = self.listvariations[variation]
            lv = self.parse_listvariations(lv)
        return lv
    
    def parse_listvariations(self,lv):
        lv = str(lv)
        lv = lv.replace("=",":=,")
        lv = lv.replace(' ',',')
        lv = lv.replace("'","")
        lv = lv.split(",")
        return lv
        
    def get_variables(self,variation=None):
        lv = self.get_lv(variation)
        variables={}
        for ii in range(np.size(lv)/2):
            variables[lv[2*ii][:-2]]=lv[2*ii+1]
        self.variables = variables
        return variables
    
    def save_data(self, data, variation):
        group = self.h5file.create_group(str(variation))
        for name, val in data.items():
            group[name] = val
            
    def get_variable_variations(variablename):
        variable_variations = []
        for variation in self.variations:
            variable_variations = data_list[variation]
    
    def get_Hparams(self, freqs, pjs, lj):
        Hparams = {}
        fzpfs = []
        
        # calculate Kerr and fzpf
        for m in range(self.nmodes):
            omega = 2*pi*freqs[m]
            lj *= 1e-9
            ej = fluxQ**2/lj
            pj = pjs[m]
            fzpf = np.sqrt(pj*hbar*omega/ej)
            fzpfs.append(fzpf)
            Hparams['fzpf_'+str(m)] = fzpf
            alpha = 2*ej/fact(4)*nck(4,2)*(fzpf**4)/hbar
            Hparams['alpha_'+str(m)] = alpha
            Hparams['freq_'+str(m)]=(omega-alpha)/2/pi

        # calculate chi
        for m in range(self.nmodes):
            for n in range(m):
                chi_mn = ej/hbar*(fzpfs[m]*fzpfs[n])**2
                Hparams['chi_'+str(m)+'_'+str(n)] = chi_mn

        return Hparams
        
    def do_bbq(self, LJvariablename, variations=None):
        # list of data dictionaries. One dict per optimetric sweep.        
        data_list = []
        data = {}
        
        # A variation is a combination of project/design 
        # variables in an optimetric sweep
        if variations is None:
            variations = range(self.nvariations)
        self.variations = variations
            
        for variation in variations: 
            # get variable values (e.g $pad_length etc.)
            data.update(self.get_variables(variation=variation))
            
            # get LJ value
            lj = eval(data[LJvariablename][:-2])
            
            # calculate participation ratio for each mode for this variation
            #pjs = self.get_p_j(variation=variation)  
            pjs = {'pj_1': 0.465358987252528,
                   'pj_2': 0.0667792007677874,
                   'pj_3': 0.00380017499629142,
                   'pj_4': 0.00117714836568518}
            self.pjs = pjs


            data.update(pjs)
            data['nmodes'] = self.nmodes

            # get bare freqs from HFSS
            data.update(self.get_freqs_bare())            

            # get Kerrs and chis
            data.update(self.get_Hparams(self.freqs_bare.values(), self.pjs.values(), lj))
            
            self.data = data
            data_list.append(data)
            self.data_list = data_list
            
            # save hdf5 containing data = variables, chis, kerrs, freqs, 
            self.save_data(data, variation)
        
        self.h5file.close()
        plot_Hparams(data_filename=self.data_filename, variations=self.variations)
        return 

class BbqAnalysis(data_filename, variations = None):
    self.data_filename = data_filename
    self.h5data = h5py.File(data_filename, 'r')

    if variations is None:
        variations = self.h5data.keys()
    self.variations = variations
    
    self.nmodes = h5data[self.variations[0]]['nmodes']        
        
    def get_swept_variables(self):
        swept_variables_names = []
        swept_variables_values = []
        for name in self.h5data[self.variations[0]].keys():
            variables = []
            for variation in self.variations:
                variables.append(self.h5data[variation][name])
            if len(set(variables))>1:
                swept_variables_names.append(name)
                swept_variables_values.append(list(set(variables)))
        return swept_variables_names, swept_variables_values
        
    def get_variable_variations(self, variablename):
        variables = []
        for variation in self.variations:
            variables.append(h5data[variation][variablename].value())
        return variables
        
    def plot_Hparams(filename, variations):
        
        fig1 = plt.subplots()
        ax1 = fig1.add_subplot(221)
        ax2 = fig1.add_subplot(222)
        ax3 = fig1.add_subplot(223)
        ax4 = fig1.add_subplot(224)
        
        for m in range(self.nmodes):        
            freq_m = 'freq_'+str(m)
            Kerr_m = 'alpha_'+str(m)
            ax1.plot(self.variations, self.get_variable_variations(freq_m), label=freq_m)
            ax4.plot(self.variations, self.get_variable_variations(Kerr_m)/2/pi, label = Kerr_m+'/2pi')
            
            for n in range(m):
                chi_m_n = 'chi_'+str(m)+'_'+str(n)
                ax3.plot(self.variations, self.get_variable_variations(chi_m_n)/2/pi, label = chi_m_n+'/2pi')
        
        swept_variables_names, swept_variables_values = self.get_swept_variables()     
        for variable in swept_variables_names:
            fig1 = plt.subplots()
            ax1 = fig1.add_subplot(221)
            ax.scatter()
        return
        
def keys(f):
    return [key for key in f.keys()]