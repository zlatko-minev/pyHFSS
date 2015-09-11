# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 11:46:56 2015

@author: Steven
"""

from hfss import *
import numpy as np

class Bbq(object):
    def __init__(self, design):
        self.design = design
        self.setup = design.get_setup
        self.nmodes = int(self.setup.n_modes)
        self.listvariations = design._solutions.ListVariations(str(self.setup.solution_name))
        self.nvariations = np.size(self.listvariations)
        self.solutions = self.setup.get_solutions()
        
    def get_p_j(self, modes=None, lv=None):
        if lv is None:
            lv=self.listvariations
        if modes is None:
            modes = [m+1 for m in range(self.nmodes)]

        P_J_vec = []
        for mode in modes:
            self.solutions.set_mode(mode, 0)
            fields = self.setup.get_fields()
            P_J = fields.P_J
            P_J.evaluate(lv=lv)
            P_J_vec.append(P_J.evaluate())
        return P_J_vec
    
    def get_freqs(self):
        return self.solutions.eigenmodes()
    #P_J_vec = get_P_J()
    
    #oModule = design._solutions
    oModule = solutions
    listvariations = oModule.ListVariations(str(self.setup.solution_name))
    #print listvariations
    
    def parselistvariations(self,lv):
        lv = str(lv)
        lv = lv.replace("=",":=,")
        lv = lv.replace(' ',',')
        lv = lv.replace("'","")
        lv = lv.split(",")
        return lv
        
    def get_vardict(self,lv):
        vardict={}
        for ii in range(np.size(lv)/2):
            vardict[lv[2*ii][:-2]]=lv[2*ii+1]
        return vardict
        
    def append_variableVec(self, vardict, variableVec, variablename):
        for ii in np.size(variableVec):
            vardict[variablename+'_'+str(ii)] = variableVec[ii]
        return vardict
        
    def saveHparams():
        
    def plotHparams():
        
    def getHparams():
        listHparams = []
        for variation in range(self.nvariations):
            lv = self.parselistvariations(self.listvariations[variation])
            vardict = self.get_vardict(lv)
            
            P_J_vec = self.get_p_j(lv)            
            vardict = self.append_variableVec(vardict, P_J_vec, 'P_J')

            freqs = self.get_freqs()
            vardict = self.append_variableVec(vardict, freqs, 'freq')
            
            
            
            listHparams.append(vardict)
        self.saveHparams()       
        
        return listHparams
        