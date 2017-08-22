'''
    @author: Zlatko Minev
    Updated 2017-02, created 2016
'''
from __future__ import print_function    # Python 2.7 and 3 compatibility
import bbq
import numpy  as np
import matplotlib.pyplot as plt
import pandas as pd

from bbq     import BbqAnalysis, print_color, CalcObject, ureg
from toolbox import DataFrame_col_diff, isint

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## HFSS design: select
project_path  = r"C:\\Users\\rslqulab\Desktop\\Lysander\participation_ratio_project\\Shyam's autonomous stabilization simulations\\"
project_name  = r'2017-08 Zlakto - (2013-11-21-2013-04-24_cooldown_2qubit)' # r'2013-11-21-2013-04-24_cooldown_2qubit'
design_name   = r'01b - 2qubit device EM'

## HFSS desgin: describe junction parameters
junc_rects    = ['qubitAlice','qubitBob'] #['juncH','juncV']           # Name of junction rectangles in HFSS
junc_lines    = ['alice_line','bob_line']#['juncH_line','juncV_line']  # Name of lines in HFSS used to define the current orientation for each junction
junc_LJ_names = ['LJAlice','LJBob'] #['LJ2','LJ1']               # Name of junction inductance variables in HFSS. DO NOT SUSE Global names that start with $.
junc_lens     = [0.0001]*2                   # This is in meters

## Analaysis:
swp_var       = 'indshiftB' #'inductance_shift'  # the name of the swept variable set to None if none
cos_trunc     = 10
fock_trunc    = 7

if 1:
    ## Run:         Connect to HFSS & analyze
    app, desktop, project = bbq.load_HFSS_project(project_name, project_path)
    design = project.get_design(design_name) if design_name != None else project.get_active_design()

    bbp = bbq.Bbq(project, design, append_analysis=False)
    bbp.do_eBBQ(junc_rect=junc_rects, junc_lines = junc_lines,  junc_len = junc_lens, junc_LJ_var_name = junc_LJ_names)


if 1:
    ## Load results
    bba            = BbqAnalysis(bbp.data_filename) # load the data (alternativly, could use  bbp.bbq_analysis)
    sol            = bba.sols
    meta_datas     = bba.meta_data
    hfss_variables = bba.hfss_variables

if 0: # Analyze a single variation
    variation = None
    CHI_O1, CHI_ND, PJ, Om, EJ, diff, LJs, SIGN, f0s, f1s, fzpfs, Qs, varz = \
        bba.analyze_variation('0', cos_trunc = cos_trunc, fock_trunc  = fock_trunc)

#%% ANALYZE ALL VARIATIONS  ==============================================================================
if 1:
    if swp_var is not None:
        RES = []
        SWP = []
        for key, s in sol.items():
            print( '\r Analyzing ', key,)
            try:
                SWP  += [ ureg.Quantity(hfss_variables[key]['_'+swp_var]).magnitude ]
                RES  += [ bbq.eBBQ_Pmj_to_H_params(s,  meta_datas[key],
                                                   cos_trunc  = cos_trunc,
                                                   fock_trunc = fock_trunc,
                                                   use_1st_order = False) ]
            except Exception as e:
                print_color(" ! ERROR %s" % e)


#%% PLOT ALL VARIATIONS  ==============================================================================
if 1:
    import matplotlib.gridspec as gridspec;
    fig = plt.figure(num = 1, figsize=(19,5))
    gs1 = gridspec.GridSpec(1, 4, width_ratios=[2,2,2,1]); gs1.update(left=0.05, right=0.98)  # wspace=0.05
    ax1 = plt.subplot(gs1[0]); ax2 = plt.subplot(gs1[1]); ax3 = plt.subplot(gs1[2]); ax3b = plt.subplot(gs1[3])

    ax = ax1; ID = 1;
    args = {'lw':0,'marker':'o','ms':5}
    ax.plot(SWP, [r[ID][0,1]for r in RES], label = '$\\chi_{DB}$', **args)
    ax.plot(SWP, [r[ID][0,2]for r in RES], label = '$\\chi_{DC}$', **args)
    ax.plot(SWP, [r[ID][1,2]for r in RES], label = '$\\chi_{BC}$', **args)
    ax.set_ylim([0.01,10**2]); ax.set_xlabel(swp_var); ax.set_title('cross-Kerr'); ax.set_ylabel('$\\chi$ (MHz)'); ax.legend(loc = 0)
    ax.set_yscale('log');   ax.set_ylim(0.1,100)
    ax1.axhspan(5.5,6.5, alpha =0.4, color= 'b')
    ax1.axhline(0.5, alpha =0.4, color= 'b')
    ax1.axhspan(85,150, alpha =0.4, color= 'b')
    ax = ax2;
    ax.plot(SWP, [r[ID][0,0] for r in RES], label = '$\\alpha_{D}$', **args)
    ax.plot(SWP, [r[ID][1,1] for r in RES], label = '$\\alpha_{B}$', **args)
    ax.plot(SWP, [r[ID][2,2] for r in RES], label = '$\\alpha_{C}$', **args)
    ax.set_ylim([100 +0.01,3.5*10**2]); ax.set_xlabel(swp_var); ax.set_title('Anharmonicity');ax.set_ylabel('$\\alpha$ (MHz)'); ax.legend(loc = 0)
    ax.set_yscale('linear')
    ax.axhline(5.238)
    ax = ax3;
    ax.plot(SWP, [r[9][:2]*10**-9 for r in RES],  **args)
    ax.plot(SWP, [r[8][:2]        for r in RES],  **{'lw':0,'marker':'x','ms':3})
    ax.set_xlabel(swp_var); ax.set_ylabel('Freq1 (GHz)'); ax.set_title('Freq.'); ax.legend(['D','B','C'], loc= 0); ax.grid(axis='y',color='gray', linestyle='-', linewidth=0.8, alpha =0.4)
    ax = ax3b
    ax.plot(SWP, [r[11] for r in RES],  **args)
    ax.set_xlabel(swp_var); ax.set_ylabel('Q'); ax.legend(['D','B','C'], loc = 0)
    try:
        ax.set_yscale('log')
    except Exception as e:
        pass
    ax.set_title('Quality')

if 0:
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

    ax   = plt.subplot(gs1[1,0]); ID =1;
    chiDC = np.array([r[ID][0,2] for r in RES])
    chiDB = np.array([r[ID][0,1] for r in RES])
    #print_color("chiDB/chiDC ratios:"); print  chiDB/chiDC
    ax.plot(SWP, chiDB/chiDC, **args); ax.locator_params(nbins=4); ax.grid(); ax.set_ylabel('$\\chi_{DB}/\\chi_{DC}$')
    ax.set_xlabel(swp_var);

    # plot the chis again
    plt.close(3);   ID = 1;
#    fig, (ax7,ax8,ax9) = plt.subplots(3,1,sharex = True, num = 3, figsize=(6,7)) ;
#
#    ax7.plot(SWP, [r[ID][0,1]for r in RES], label = '$\\chi_{DB}$', c = 'b', **args); ax7.set_ylabel('$\\chi_{DB}$ (MHz)');
#    ax9.plot(SWP, [r[ID][0,2]for r in RES], label = '$\\chi_{DC}$', c = 'g', **args); ax9.set_ylabel('$\\chi_{DC}$ (MHz)');
#    ax8.plot(SWP, [r[ID][1,2]for r in RES], label = '$\\chi_{BC}$', c = 'r', **args); ax8.set_ylabel('$\\chi_{BC}$ (MHz)');
#    ax9.set_xlabel(swp_var); ax7.set_title('cross-Kerr');
#    #ax7.axhspan(85,150,  alpha =0.4, color= 'b')
#    ax8.axhspan(5.5,6.5, alpha =0.4, color= 'b')
#    ax9.axhline(0.5,     alpha =0.4, color= 'b')
#    fig.tight_layout()
#

if 1:
    ID = 1;
    fig, (ax7,ax8,ax9) = plt.subplots(3,1,sharex = True, num = 3, figsize=(6,7)) ;

    #RES = RES0
    args = {'lw':0,'marker':'o','ms':4}
    ax7.plot(SWP, [r[ID][0,1]for r in RES], label = '$\\chi_{DB}$', c = 'b', **args); ax7.set_ylabel('$\\chi_{DB}$ (MHz)');
    ax9.plot(SWP, [r[ID][0,2]for r in RES], label = '$\\chi_{DC}$', c = 'g', **args); ax9.set_ylabel('$\\chi_{DC}$ (MHz)');
    ax8.plot(SWP, [r[ID][1,2]for r in RES], label = '$\\chi_{BC}$', c = 'r', **args); ax8.set_ylabel('$\\chi_{BC}$ (MHz)');
    #RES = RES1
    args = {'lw':0,'marker':'x','ms':10}
    ax7.plot(SWP, [r[ID][0,1]for r in RES], label = '$\\chi_{DB}$', c = 'b', **args); ax7.set_ylabel('$\\chi_{DB}$ (MHz)');
    ax9.plot(SWP, [r[ID][0,2]for r in RES], label = '$\\chi_{DC}$', c = 'g', **args); ax9.set_ylabel('$\\chi_{DC}$ (MHz)');
    ax8.plot(SWP, [r[ID][1,2]for r in RES], label = '$\\chi_{BC}$', c = 'r', **args); ax8.set_ylabel('$\\chi_{BC}$ (MHz)');
    ax9.set_xlabel(swp_var); ax7.set_title('cross Kerr');
    #ax7.axhspan(85,150,  alpha =0.4, color= 'b')
    ax8.axhspan(5.5,6.5, alpha =0.4, color= 'b')
    ax9.axhline(0.5,     alpha =0.4, color= 'b')


#%%
if 1: # plot mesh
    fig = plt.figure(8);  fig.clf()
    tets = bba.get_convergences_max_tets()
    varsz  = bba.get_variable_vs(swp_var)
    Y = {}
    for key in tets.keys():
        Y[varsz[key]] = tets[key]
    y =  pd.Series(Y) #.values(), index = varsz.values())
    y.plot(marker = '*', ms = 20)
    ax7t = ax7.twinx()
    ax7t.plot(y, marker = '*', ms = 10, c = 'g')
#%%
if 1:
    fig = plt.figure(21); fig.clf()
    tts = bba.get_convergences_Tets_vs_pass()
    for key, x in tts.items():
        #np.log10(x).plot(label = varsz[key])
        x.plot(label = varsz[key])
    plt.legend(loc = 0)

#%%
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
