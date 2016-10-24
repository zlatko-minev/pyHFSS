# -*- coding: utf-8 -*-
"""
Created on Fri Jul 01 12:16:56 2016

@author: rslqulab
"""

if 1:
    import bbq, matplotlib.pyplot as plt
    from bbq import BbqAnalysis
    import pandas as pd, matplotlib.pyplot as plt, numpy as np;
    import hfss, bbq, bbqNumericalDiagonalization
    from hfss import ureg
    from bbq  import eBBQ_Pmj_to_H_params, print_color, print_matrix
    import matplotlib.gridspec as gridspec;
    import glob, os
    import matplotlib as mpl
    import matplotlib.gridspec as gridspec
    
    mpl.rc('axes', grid = True)

    swp_var   = 'scale_factor'
    swp_slice = [0.9,1.01]
    pass_per_file = 2
    folder   = u'\\\\?\\C:\\Users\\rslqulab\\Desktop\\zkm\\2016_qm_jumps\\DiTransmon_Asymetric\\'
    folder  += u'mesh_analysis_project\\20pass_7k_seed- 11. SHANTANU FAB SM22 Col 1 Row 9 Measured [May 6 2016]\\' ; import os; os.chdir(folder);
    
    ### Get files in directory sorted from oldest to newest
    files = glob.glob("*.hdf5")
    files.sort(key=os.path.getmtime);  #files.reverse() 
    #print("\n".join(files))
    
    files   = files[-11:-1]
    
    ### Load all files
    A = []
    for i, filename in enumerate(files):
        print filename
        file_name = folder + filename 
        A += [ BbqAnalysis(file_name) ] 
        
if 1: ### Calculates percent error in freq 
    plt.close('all')    
    plt.figure(2)
    gs = gridspec.GridSpec(len(A), 1)
    modei, rng = 0, [4.8,5]
    
    frq_mean_vs_pass = {}
    frq_std_vs_pass  = {}
    mesh_tot_vs_pass = {}
    for i, a in enumerate(A):
        pass_num = pass_per_file * i +1 
        fs = a.get_Fs(swp_var)[swp_slice[0]:swp_slice[1]]
        print fs.shape
        frq_mean_vs_pass[pass_num]  = [fs[ii].mean() for ii in range(fs.shape[1])]
        frq_std_vs_pass[pass_num]   = [fs[ii].std()  for ii in range(fs.shape[1])]
        mesh_tot_vs_pass[pass_num]  = np.mean(a.get_mesh_tot().values())
        ax = plt.subplot(gs[i]); 
        ax.hist(fs[modei], bins = 30, histtype='bar', stacked=False, range = rng, label = str(pass_num),\
                    alpha = 0.6); 
        ax.legend(prop={'size': 10})
    frq_mean_vs_pass = pd.DataFrame(frq_mean_vs_pass).T
    frq_std_vs_pass  = pd.DataFrame(frq_std_vs_pass).T
    frq_err_vs_pass  = 100 *  frq_std_vs_pass / frq_mean_vs_pass
    mesh_tot_vs_pass = pd.Series(mesh_tot_vs_pass)

    plt.figure(1)
    gs = gridspec.GridSpec(2, 2)
    #ax1 = plt.subplot(gs[0, 0]); 
    gs1 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs[0, 0])
    ax1 = [plt.subplot(gs1[i]) for i in range(3)]
    ax2 = plt.subplot(gs[1, 0]); ax3 = plt.subplot(gs[0, 1]); ax4 = plt.subplot(gs[1, 1])
    for i in range(3): ax1[i].plot(frq_mean_vs_pass[i]);  
    ax1[0].set_title('Mean Freq.');
    ax2.plot(frq_err_vs_pass);   ax2.set_title('Freq. err %'); ax2.set_xlabel('Pass Number');
    ax3.plot(mesh_tot_vs_pass);  ax3.set_title('Total tets');
    ax4.set_xlabel('Pass Number');
    
    
    
    