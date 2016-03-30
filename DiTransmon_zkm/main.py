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
import bbq, matplotlib.pyplot as plt, numpy as np;  from bbq import print_color
from IPython.display import display 
#from scipy.constants import *
#plt.close('all')
    
if 1:    
    proj_name    = r'2016_03_28_Di_Transmon from Antonio' 
    project_path = 'C:\\Users\\rslqulab\\Desktop\\zkm\\2016_qm_jumps\\DiTransmon_Asymetric\\'
    app, desktop, project = load_HFSS_project(proj_name, project_path)
    design       = project.get_active_design() #get_design("CoaxCav")
    
    bbp = bbq.Bbq(project, design, append_analysis=False, calculate_H=True)

if 0:
    junc_rects    = ['juncV','juncH'] # ['qubit1','qubit2']
    junc_LJ_names = ['LJ1', 'LJ2'];
    junc_lens     = [0.0001]*2 
    bbp.do_bbq(calculate_H = False,  plot_fig = False,
               Pj_from_current= True, junc_rect=junc_rects, junc_len = junc_lens, junc_LJ_var_name = junc_LJ_names )
    sol = bbp.PJ_multi_sol


if 0:
    def comp_CHI_matrix(s):
        ''' s      = sol['0'];  CHI, PJ, Om, EJ, diff = comp_CHI_matrix(s)'''
        import  scipy;    Planck  = scipy.constants.Planck
        LJs        = s.loc[0,s.keys().str.contains('LJs')] # LJ in nH
        EJs        = (bbq.fluxQ**2/LJs/Planck*10**-9).astype(np.float)        # EJs in GHz
        PJ_Jsu     = s.loc[:,s.keys().str.contains('pJ')]  # EPR from Jsurf avg
        PJ_Jsu_sum = PJ_Jsu.apply(sum, axis = 1)           # sum of participations as calculated by avg surf current 
        PJ_glb_sum = (s['U_E'] - s['U_H'])/(2*s['U_E'])    # sum of participations as calculated by global UH and UE  
        diff       = (PJ_Jsu_sum-PJ_glb_sum)/PJ_glb_sum*100# debug
        #TODO: figure out the systematic
        # print '% diff b/w Jsurf_avg & global Pj:'; display(diff)
        if 1:# Renormalize: to sum to PJ_glb_sum; so that PJs.apply(sum, axis = 1) - PJ_glb_sum =0
            PJs = PJ_Jsu.divide(PJ_Jsu_sum, axis=0).mul(PJ_glb_sum,axis=0)
        else: PJs = PJ_Jsu
        PJ    = np.mat(PJs.values)
        Om    = np.mat(np.diagflat(s['freq'].values)) 
        EJ    = np.mat(np.diagflat(EJs.values))
        CHI   = Om * PJ * EJ.I * PJ.T * Om * 1000 # MHz
        return CHI, PJ, Om, EJ, diff
    #TODO: phi_zpf -> full BBQ 
    
#==============================================================================
#     Plot results for sweep
#==============================================================================
    swpvar='join_h'    
    RES = []; SWP = [];
    for key, s in sol.iteritems():
#        bbp.
        varz  = bbp.get_variables(variation=key)
        SWP  += [ eval(varz['_'+swpvar][:-2]) ] 
        RES  += [ comp_CHI_matrix(s) ]
    import matplotlib.gridspec as gridspec;
    
    fig = plt.figure(num = 1, figsize=(15,5)) 
    gs1 = gridspec.GridSpec(1, 3, width_ratios=[1,1,1]); gs1.update(left=0.05, right=0.95)  # wspace=0.05
    ax1 = plt.subplot(gs1[0]); ax2 = plt.subplot(gs1[1]); ax3 = plt.subplot(gs1[2])
    
    ax = ax1
    args = {'lw':0,'marker':'o','ms':5}
    ax.plot(SWP, [r[0][0,1]for r in RES], label = '$\\chi_{DB}$', **args)
    ax.plot(SWP, [r[0][0,2]for r in RES], label = '$\\chi_{DC}$', **args)
    ax.plot(SWP, [r[0][1,2]for r in RES], label = '$\\chi_{BC}$', **args)
    ax.set_ylim([0.01,10**2]); ax.set_xlabel(swpvar); ax.set_ylabel('$\\chi$ (MHz)'); ax.legend(loc = 0)
    ax.set_yscale('log')
    ax = ax2
    ax.plot(SWP, [r[0][0,0]/2 for r in RES], label = '$\\alpha_{D}$', **args)
    ax.plot(SWP, [r[0][1,1]/2 for r in RES], label = '$\\alpha_{B}$', **args)
    ax.plot(SWP, [r[0][2,2]/2 for r in RES], label = '$\\alpha_{C}$', **args)
    ax.set_ylim([100 +0.01,3.5*10**2]); ax.set_xlabel(swpvar); ax.set_ylabel('$\\alpha$ (MHz)'); ax.legend(loc = 0)
    ax.set_yscale('linear')
    ax = ax3
    ax.plot(SWP, [np.diag(r[2]) for r in RES],  **args)
    ax.set_xlabel(swpvar); ax.set_ylabel('Freq. (GHz)');  ax.legend(['D','B','C'], loc= 0)
    
        
if 0: 
    variation = '0';  pJ_method = 'J_surf_mag';
    #pJ_mj_series = bbp.calc_Pjs_from_I_for_mode(variation, bbp.U_H,bbp.U_E, bbp.LJs, junc_rects, junc_lens, method = pJ_method) # to be implemented          
    res = bbp.calc_avg_current_J_surf_mag(variation,junc_rects[0], junc_lens[0])
    
if 0: # for debug 
    variation = '0'; junc_rect = 'juncV';
    print_color(' Setup: ' + bbp.setup.name)
    lv = bbp.get_lv(variation)
    calc = CalcObject([],bbp.setup)
    calc = calc.getQty("Jsurf").mag().integrate_surf(name = junc_rect)
 
    #bbp.calc_avg_current_J_surf_mag('0','juncV',1)