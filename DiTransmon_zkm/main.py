# -*- coding: utf-8 -*-
"""
see http://arrc.ou.edu/~cody/hfsslib/hfss_examples/#
"""
#%%
import sys;  IMP_PATH = r'C:\\Users\\rslqulab\\Desktop\\zkm\\pyHFSS\\';
if ~(IMP_PATH in sys.path): sys.path.append(IMP_PATH);
    
import hfss;
from hfss import CalcObject
from hfss import load_HFSS_project
import bbq, matplotlib.pyplot as plt, numpy as np; 
#from scipy.constants import *
#plt.close('all')
    
if 1:    
    proj_name    = r'2016_03_28_Di_Transmon from Antonio' 
    project_path = 'C:\\Users\\rslqulab\\Desktop\\zkm\\2016_qm_jumps\\DiTransmon_Asymetric\\'
    app, desktop, project = load_HFSS_project(proj_name, project_path)
    design       = project.get_active_design() #get_design("CoaxCav")
    
    bbp = bbq.Bbq(project, design, append_analysis=False, calculate_H=True)

if 1:
    junc_rects    = ['juncV','juncH']
    junc_LJ_names = ['LJ1', 'LJ2'];
    junc_lens     = [0.0001]*2 
    bbp.do_bbq('LJ1', calculate_H = False,  plot_fig = False,
               Pj_from_current= True, junc_rect=junc_rects, junc_len = junc_lens, junc_LJ_var_name = junc_LJ_names )
    sol = bbp.PJ_multi_sol
              
if 0: # for debug 
    variation = '0'; self = bbp; junc_rect = 'juncV';
    print_color(' Setup: ' + self.setup.name)
    lv = bbp.get_lv(variation)
    calc = CalcObject([],self.setup)
    calc = calc.getQty("Jsurf").mag().integrate_surf(name = junc_rect)
 
    #bbp.calc_avg_current_J_surf_mag('0','juncV',1)