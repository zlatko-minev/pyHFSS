# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 14:58:53 2016

@author: Zlatko
"""

import bbq, matplotlib.pyplot as plt
from bbq import BbqAnalysis
import pandas as pd, matplotlib.pyplot as plt, numpy as np;
import hfss, bbq, bbqNumericalDiagonalization
from hfss import ureg
from bbq  import eBBQ_Pmj_to_H_params, print_color, print_matrix
import matplotlib.gridspec as gridspec;

#file_name = u'C:\\Users\\rslqulab\\Desktop\\Lysander\\\\/pin_position_sweep(perfect conductor)_4-20-16/10. SHANTANU FAB 1 [April13 2016]/10. SHANTANU FAB 1 [April13 2016]_20160427_143612.hdf5'
folder   = u'\\\\?\\C:\\Users\\rslqulab\\Desktop\\zkm\\2016_qm_jumps\\DiTransmon_Asymetric\\'
folder  += u'qubit_separation_sweep(0-150microns)_6-1-16\\11. SHANTANU FAB SM22 Col 1 Row 9 Measured [May 6 2016]\\' ; import os; os.chdir(folder);
filename = 'qubit_separation_sweep(0-150microns)_6-1-16';         #print os.listdir(folder )
plot_title    = 'Qubit Separation Sweep (0-150 microns)'

load_data   = True
plot_Qs     = True
analyze_BBQ = True
plot_chis   = True
plot_Fs     = True
plot_Pmj    = False
cos_trunc = 8;   fock_trunc  = 9;
plt.close('all')

if load_data:
    file_name = folder + filename + '.hdf5'
    swp_var   = 'qubit_distance'
    bba = BbqAnalysis(file_name)
    hfss_variables = bba.hfss_variables
    sols           = bba.sols
    meta_datas     = bba.meta_data

#%%==============================================================================
# Plot Qs
#==============================================================================
if plot_Qs:
    Qs  = bba.get_Qs(swp_var=swp_var)
    args = {'lw':0,'marker':'o','ms':5}
    Qs.plot(**args); plt.legend(['D', 'B','G' ], loc = 0)
    ax = plt.gca();
    ax.set_xlabel(swp_var);
    ax.set_title('Qs- ' + plot_title); 
    ax.set_ylabel('Q'); 
    ax.set_yscale('log'); 
    plt.gcf().savefig(folder + filename +' Q.png')
    
if plot_Fs:
    Fs  = bba.get_Fs(swp_var=swp_var)
    args = {'lw':0,'marker':'o','ms':5}
    Fs.plot(**args); plt.legend(['D', 'B','G' ], loc = 0)
    ax = plt.gca();
    ax.set_xlabel(swp_var);
    ax.set_title('Frequencies- ' + plot_title); 
    ax.set_ylabel('F (GHz)'); 
    plt.gcf().savefig(folder + filename + ' F.png')
    plt.figure(); 
    plt.subplot(311); Fs[0].plot(**args)
    plt.subplot(312); Fs[1].plot(**args)
    plt.subplot(313); Fs[2].plot(**args)
        
        

#%% 
if analyze_BBQ:  
    RES = []; SWP = [];
    for key, s in sols.iteritems():
        print '\r Analyzing ', key,
        varz  = hfss_variables[key]
        SWP  += [ ureg.Quantity(varz['_'+swp_var]).magnitude ]  
        RES  += [ eBBQ_Pmj_to_H_params(s, meta_datas[key], cos_trunc = cos_trunc, fock_trunc = fock_trunc) ]
        
#%%==============================================================================
# Plot Chis
#==============================================================================
if plot_chis:
    fig, (ax1,ax2,ax3) = plt.subplots(3, 1, figsize = (10,6), sharex = True)    
    ax3.set_xlabel('Sweep '  + swp_var)
    args = {'lw':0,'marker':'o','ms':5}
    ID = 1; 
    ax1.plot(SWP, [r[ID][0,1]for r in RES], **args); ax1.set_ylabel('$\\chi_{DB}$')
    ax2.plot(SWP, [r[ID][0,2]for r in RES], **args); ax2.set_ylabel('$\\chi_{DC}$') 
    ax3.plot(SWP, [r[ID][1,2]for r in RES], **args); ax3.set_ylabel('$\\chi_{BC}$')
    ax1.set_title('Chi- ' + plot_title);     
    fig.tight_layout()
    fig.savefig(folder + filename + ' chis.png')
    
    
#%%
if plot_Pmj:
   ### Participation Plot ###
    fig = plt.figure(num = 2, figsize=(15,5)) 
    gs1 = gridspec.GridSpec(2, 3, width_ratios=[1,1,1], height_ratios=[5,1]); gs1.update(left=0.05, right=0.95)  # wspace=0.05
    ax4 = plt.subplot(gs1[0,0]); ax5 = plt.subplot(gs1[0,1]); ax6 = plt.subplot(gs1[0,2])
    ax = ax4; ID = 2
    ax.plot(SWP, [r[ID][0,0] for r in RES], label = '$P_{DH}$', **args)
    ax.plot(SWP, [r[ID][0,1] for r in RES], label = '$P_{BV}$', **args)
    ax.set_ylabel('Participation'); ax.legend(loc = 0)
    ax = ax5
    ax.plot(SWP, [r[ID][0,1] for r in RES], label = '$P_{DV}$', **args)
    ax.plot(SWP, [r[ID][1,0] for r in RES], label = '$P_{BH}$', **args)
    ax.set_xlabel(swp_var); ax.set_ylabel('Participation'); ax.legend(loc = 0)
    ax = ax6
    ax.plot(SWP, [r[ID][2,0] for r in RES], label = '$P_{CH}$', **args)
    ax.plot(SWP, [r[ID][2,1] for r in RES], label = '$P_{CV}$', **args)
    ax.set_yscale('log'); ax.set_xlabel(swp_var); ax.set_ylabel('Participation'); ax.legend(loc = 0)
    plt.figure()
    plt.subplot(211)
    plt.plot(SWP, [r[ID][2,0] for r in RES], label = '$P_{CH}$', **args); plt.legend(loc = 0)
    plt.subplot(212)
    plt.plot(SWP, [r[ID][2,0] for r in RES], label = '$P_{CH}$', **args); plt.legend(loc = 0)
  
#    ax   = plt.subplot(gs1[1,0]); ID =1;
#    chiDC = np.array([r[ID][0,2] for r in RES])
#    chiDB = np.array([r[ID][0,1] for r in RES])
#    #print_color("chiDB/chiDC ratios:"); print  chiDB/chiDC
   # ax.plot(SWP, chiDB/chiDC, **args); ax.locator_params(nbins=4); ax.grid(); ax.set_ylabel('$\\chi_{DB}/\\chi_{DC}$')
    #ax.set_xlabel(swp_var); 