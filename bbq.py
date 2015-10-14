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
import matplotlib.pyplot as plt
from stat import S_ISREG, ST_CTIME, ST_MODE
import sys
import shutil
import time

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
    
    def __init__(self, project, design, verbose=True, append_analysis=False, calculate_H=True):
        self.project = project
        self.design = design
        self.setup = design.get_setup()
        self.nmodes = int(self.setup.n_modes)
        self.listvariations = design._solutions.ListVariations(str(self.setup.solution_name))
        self.nominalvariation = design.get_nominal_variation()
        self.nvariations = np.size(self.listvariations)
        self.solutions = self.setup.get_solutions()
        self.verbose = verbose
        self.calculate_H = calculate_H
        self.append_analysis = append_analysis
        
        self.setup_data()
        
        print 'Number of modes : ' + str(self.nmodes)
        print 'Number of variations : ' + str(self.nvariations)
        
        self.get_latest_h5()
        if self.latest_h5_path is not None and self.append_analysis:
            latest_bbq_analysis = BbqAnalysis(self.latest_h5_path)
            print 'Varied variables and values : '
            print latest_bbq_analysis.get_swept_variables()
            print 'Variations : '
            print latest_bbq_analysis.variations

    def get_latest_h5(self):
        dirpath = self.data_dir
        
        # get all entries in the directory w/ stats
        entries1 = (os.path.join(dirpath, fn) for fn in os.listdir(dirpath))
        entries2 = ((os.stat(path), path) for path in entries1)
        
        # leave only regular files, insert creation date
        entries3 = ((stat[ST_CTIME], path)
                   for stat, path in entries2 if S_ISREG(stat[ST_MODE]) and path[-4:]=='hdf5')
        #NOTE: on Windows `ST_CTIME` is a creation date 
        #  but on Unix it could be something else
        #NOTE: use `ST_MTIME` to sort by a modification date
        
        paths_sorted = []
        for cdate, path in sorted(entries3):
            paths_sorted.append(path)
            #print time.ctime(cdate), os.path.basename(path)
        if len(paths_sorted) > 0:
            self.latest_h5_path = paths_sorted[-1]
            print 'This simulations has been analyzed, latest data in ' + self.latest_h5_path
        else:
            self.latest_h5_path = None
            print 'This simulation has never been analyzed'
        
    def setup_data(self):
        data_dir = root_dir + '/' + self.project.name + '/' + self.design.name
        if not os.path.isdir(data_dir):
            os.makedirs(data_dir)
        self.data_dir = data_dir
        self.data_filename = self.data_dir + '/' + self.design.name + '_' + time.strftime('%Y%m%d_%H%M%S', time.localtime()) + '.hdf5'

    def get_p_j(self, modes=None, variation=None):
        lv = self.get_lv(variation)
        if modes is None:
            modes = range(self.nmodes)

        pjs = {}
        for ii, m in enumerate(modes):
            print 'Calculating p_j for mode ' + str(m) + ' (' + str(ii) + '/' + str(np.size(modes)) + ')'
            self.solutions.set_mode(m+1, 0)
            fields = self.setup.get_fields()
            P_J = fields.P_J
            pjs['pj_'+str(m)] = P_J.evaluate(lv=lv)
        self.pjs = pjs
        if self.verbose: print pjs
        return pjs
    
    def get_freqs_bare(self, variation):
        freqs_bare_vals = []
        freqs_bare_dict = {}
        freqs, bws = self.solutions.eigenmodes(str(self.listvariations[eval(variation)]))
        for m in range(self.nmodes):
            freqs_bare_dict['freq_bare_'+str(m)] = 1e9*freqs[m]
            freqs_bare_vals.append(1e9*freqs[m])
            if bws is not None:
                freqs_bare_dict['Q_'+str(m)] = freqs[m]/bws[m]
        self.freqs_bare = freqs_bare_dict
        if self.verbose: print freqs_bare_dict
        return freqs_bare_dict, freqs_bare_vals
        
        
    def get_lv(self, variation):
        if variation is None:
            lv = self.nominalvariation
            lv = self.parse_listvariations(lv)
        else:
            lv = self.listvariations[eval(variation)]
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
            variables['_'+lv[2*ii][:-2]]=lv[2*ii+1]
        self.variables = variables
        return variables
    
    def save_data(self, data, variation):
        group = self.h5file.create_group(variation)
        for name, val in data.items():
            group[name] = val
    
    def get_Hparams(self, freqs, pjs, lj):
        Hparams = {}
        fzpfs = []
        
        # calculate Kerr and fzpf
        for m in range(self.nmodes):
            omega = 2*pi*freqs[m]
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
        
    def do_bbq(self, LJvariablename, variations=None, plot_fig=True):
        
        if self.latest_h5_path is not None and self.append_analysis:
            shutil.copyfile(self.latest_h5_path, self.data_filename)
        
        self.h5file = h5py.File(self.data_filename)
        
        # list of data dictionaries. One dict per optimetric sweep.        
        data_list = []
        data = {}
        
        # A variation is a combination of project/design 
        # variables in an optimetric sweep
        if variations is None:
            variations = [str(i) for i in range(self.nvariations)]
        self.variations = variations

        for ii, variation in enumerate(variations):
            print 'variation : ' + variation + ' (' + str(ii) + ' / ' + str(len(self.variations)) + ')'
            
            if variation in self.h5file.keys() and self.append_analysis:
                print 'variation previously analyzed ...'
                continue
            print 'variation ' + variation + ' NOT analyzed'
            time.sleep(1)

            #return
            # get variable values (e.g $pad_length etc.)
            data.update(self.get_variables(variation=variation))
            
            # get LJ value
            if self.calculate_H:
                lj = eval(data['_'+LJvariablename][:-2])
                lj *= 1e-9
                
                #calculate participation ratio for each mode for this variation
                pjs = self.get_p_j(variation=variation)
    
                data.update(pjs)
            data['nmodes'] = self.nmodes

            # get bare freqs from HFSS
            freqs_bare_dict, freqs_bare_vals = self.get_freqs_bare(variation)
            data.update(freqs_bare_dict)

            # get Kerrs and chis
            if self.calculate_H:
                data.update(self.get_Hparams(freqs_bare_vals, self.pjs.values()[::-1], lj))
            
            self.data = data
            data_list.append(data)
            self.data_list = data_list
            
            # save hdf5 containing data = variables, chis, kerrs, freqs, 
            self.save_data(data, variation)
        
        self.h5file.close()
        self.bbq_analysis = BbqAnalysis(self.data_filename, variations=self.variations)
        if plot_fig:
            self.bbq_analysis.plot_Hparams()
        return

class BbqAnalysis(object):
    ''' defines an analysis object which loads and plots data from a h5 file
    This data is obtained using e.g bbq.do_bbq
    '''    
    def __init__(self, data_filename, variations = None):
        self.data_filename = data_filename
        self.h5data = h5py.File(data_filename, 'r')
    
        if variations is None:
            variations = self.h5data.keys()
        self.variations = variations
        
        self.nmodes = self.h5data[self.variations[0]]['nmodes'].value      
        
    def get_swept_variables(self):
        swept_variables_names = []
        swept_variables_values = []
        for name in self.h5data[self.variations[0]].keys():
            if '_'==name[0]: # design variables all start with _
                variables = []
                for variation in self.variations:
                    variables.append(self.h5data[variation][name].value)
                if len(set(variables))>1:
                    swept_variables_names.append(name)
                    swept_variables_values.append(list(set(variables)))
            else:
                pass
        return swept_variables_names, swept_variables_values
        
    def get_variable_variations(self, variablename):
        variables = []
        for variation in self.variations:
            variables.append(self.h5data[variation][variablename].value)
        return np.asarray(variables)
    
    def get_float_units(self, variable_name, variation='0'):
        variable_value = self.h5data[variation][variable_name].value
        n = 1
        try:
            float(variable_value)
            return float(variable_value), ''
        except ValueError:
            while True:
                try:
                    float(variable_value[:-n])
                    return float(variable_value[:-n]), variable_value[len(variable_value)-n:]
                except:
                    n+=1
                    
    def plot_Hparams(self, variable_name = None):
        fig, ax = plt.subplots(2,2, figsize=(24,10))

        if variable_name == None:
            xaxis = self.variations
        else:
            xaxis = []
            for variation in self.variations:
                xaxis.append(self.get_float_units(variable_name, variation)[0])
    
        for m in range(self.nmodes):        
            freq_m = 'freq_'+str(m)
            Kerr_m = 'alpha_'+str(m)
            Q_m = 'Q_'+str(m)
            if freq_m not in self.h5data[self.variations[0]].keys():
                freq_m = 'freq_bare_'+str(m)
            else:
                pass
            if Kerr_m in self.h5data[self.variations[0]].keys():
                ax[0][1].plot(xaxis, self.get_variable_variations(Kerr_m)/2/pi/1e6, 'o', label = str(m))
            else:
                pass
    
            ax[0][0].plot(xaxis, self.get_variable_variations(freq_m)/1e9, 'o', label=str(m))            
            if Q_m in self.h5data[self.variations[0]].keys():             
                ax[1][1].plot(xaxis, self.get_variable_variations(Q_m), 'o', label = str(m))
            else:
                pass
            
            for n in range(m):
                chi_m_n = 'chi_'+str(m)+'_'+str(n)
                if chi_m_n in self.h5data[self.variations[0]].keys():
                    ax[1][0].plot(xaxis, self.get_variable_variations(chi_m_n)/2/pi/1e6, 'o', label=str(m)+','+str(n))
        
        ax[0][0].legend()
        ax[0][0].set_ylabel('freq (GHz)')
        
        ax[0][1].legend()
        ax[0][1].set_ylabel('Kerr/2pi (MHz)')
        ax[0][1].set_yscale('log')
        
        ax[1][0].legend()
        ax[1][0].set_ylabel('Chi/2pi (MHz)')
        ax[1][0].set_yscale('log')
        
        ax[1][1].legend()
        ax[1][1].set_ylabel('Q')
        ax[1][1].set_yscale('log')
        
        if variable_name == None:
            swept_variables_names, swept_variables_values = self.get_swept_variables()
            xticks = []
            for variation in xaxis:
                xtick = ''
                for name in swept_variables_names:
                    xtick += name[1:] + ' = ' + self.h5data[variation][name].value + '\n'
                xticks.append(str(xtick))
            ax[1][0].set_xticks([int(v) for v in xaxis])
            ax[1][0].set_xticklabels(xticks, rotation='vertical')
            ax[1][1].set_xticks([int(v) for v in xaxis])
            ax[1][1].set_xticklabels(xticks, rotation='vertical')
            
            ax[0][0].set_xticklabels([])
            ax[0][1].set_xticklabels([])
        else:
            xlabel = variable_name + ' (' + self.get_float_units(variable_name, self.variations[0])[1] + ')'
            ax[1][0].set_xlabel(xlabel)
            ax[1][1].set_xlabel(xlabel)

        fig.subplots_adjust(bottom=0.3)
        fig.suptitle(self.data_filename)
        fig.savefig(self.data_filename[:-5]+'.jpg')

        return fig, ax
            
            
            
#        for variable in swept_variables_names:
#            fig1 = plt.subplots()
#            ax1 = fig1.add_subplot(221)
#            ax.scatter()
#        return