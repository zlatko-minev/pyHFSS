#%% @author: Zlatko
import sys;  IMP_PATH = r'C:\\Users\\rslqulab\\Desktop\\zkm\\pyHFSS\\';
if ~(IMP_PATH in sys.path): sys.path.insert(0,IMP_PATH);

import pandas as pd, matplotlib.pyplot as plt, numpy as np;
import hfss, bbq, bbqNumericalDiagonalization
from hfss import CalcObject, ureg, load_HFSS_project
from bbq  import eBBQ_Pmj_to_H_params, print_color, print_matrix


if 1:    
    proj_name    = r'DTW1_shantanu_12_17_2015' 
    project_path = 'C:\\Users\\rslqulab\\Desktop\\Lysander\\participation_ratio_project\\2q_1w\\'
    #proj_name    = r'pin_shift_sweep(circuit=-1500microns_perfect conductor)_5-2-16' 
    #project_path = 'C:\\Users\\rslqulab\\Desktop\Lysander\\'
    app, desktop, project = load_HFSS_project(proj_name, project_path)
    design       = project.get_design("5- SMA Adapter SMA Pin 3Pad design Moved sapphire1") #project.get_active_design(10. SHANTANU FAB 1 [April13 2016])  #
    
    bbp = bbq.Bbq(project, design, append_analysis=False)

if 1:
    junc_rects    = ['qubit1','qubit2']#['juncV','juncH'] 
    junc_lines    = ['line_qubit1','line_qubit2']#['juncV_line','juncH_line'] #
    junc_LJ_names = ['LJ1', 'LJ2'];
    junc_lens     = [0.0001]*2     #this is in meters!                                                  # this can soon be replaced by intgrating over junc_lines 
    bbp.do_eBBQ(junc_rect=junc_rects, junc_lines = junc_lines,  junc_len = junc_lens, junc_LJ_var_name = junc_LJ_names)
    bba = bbp.bbq_analysis 
    
    sol           = bba.sols
    meta_datas    = bba.meta_data
    hfss_variables= bba.hfss_variables

#%%
if 1:
    reload(bbq)   # so that we can make changes to the class and reload quickly 
    from bbq  import BbqAnalysis
    bba = BbqAnalysis(bbp.data_filename)
    
    cos_trunc = 10;   fock_trunc  = 7;
    CHI_O1, CHI_ND, PJ, Om, EJ, diff, LJs, SIGN, f0s, f1s, fzpfs, Qs, varz = \
        bba.analyze_variation(variation = '0', cos_trunc = cos_trunc,   fock_trunc  = fock_trunc)
#    s         = sol[variation];   
#    meta_data = meta_datas[variation]
#    varz      = hfss_variables[variation]    
#    CHI_O1, CHI_ND, PJ, Om, EJ, diff, LJs, SIGN, f0s, f1s, fzpfs, Qs = \
#        eBBQ_Pmj_to_H_params(s, meta_data, cos_trunc = cos_trunc, fock_trunc = fock_trunc)
#    
#    print '\nPJ=\t(renorm.)';        print_matrix(PJ*SIGN, frmt = "{:7.4f}")
#    #print '\nCHI_O1=\t PT. [alpha diag]'; print_matrix(CHI_O1,append_row ="MHz" )
#    print '\nf0={:6.2f} {:7.2f} {:7.2f} GHz'.format(*f0s)
#    print '\nCHI_ND=\t PJ O(%d) [alpha diag]'%(cos_trunc); print_matrix(CHI_ND, append_row ="MHz")
#    print '\nf1={:6.2f} {:7.2f} {:7.2f} GHz'.format(*(f1s*1E-9))   
#    print 'Q={:8.1e} {:7.1e} {:6.0f}'.format(*(Qs))
    #print pd.Series({ key:varz[key] for key in ['_join_w','_join_h','_padV_width', '_padV_height','_padH_width', '_padH_height','_scaleV','_scaleH', '_LJ1','_LJ2'] })

#%%==============================================================================
#     Plot results for sweep
#==============================================================================
if 1:
    swpvar='LJ'    
    use_1st_order = True  # use 1st O PT  to identify correct eigenvectors in ND
    RES = []; SWP = [];
    for key, s in sol.iteritems():     
        print '\r Analyzing ', key,
        try:
            varz  = hfss_variables[key]
            SWP  += [ ureg.Quantity(varz['_'+swpvar]).magnitude ]  
            RES  += [ eBBQ_Pmj_to_H_params(s, meta_datas[key], cos_trunc = cos_trunc, fock_trunc = fock_trunc, use_1st_order = False) ]
        except Exception as e:
            print_color(" !ERROR %s" % e)
    import matplotlib.gridspec as gridspec;
    #%%
    fig = plt.figure(num = 1, figsize=(19,5)) 
    gs1 = gridspec.GridSpec(1, 4, width_ratios=[2,2,2,1]); gs1.update(left=0.05, right=0.98)  # wspace=0.05
    ax1 = plt.subplot(gs1[0]); ax2 = plt.subplot(gs1[1]); ax3 = plt.subplot(gs1[2]); ax3b = plt.subplot(gs1[3])
    
    ax = ax1; ID = 1; 
    args = {'lw':0,'marker':'o','ms':5}
    ax.plot(SWP, [r[ID][0,1]for r in RES], label = '$\\chi_{DB}$', **args)
    ax.plot(SWP, [r[ID][0,2]for r in RES], label = '$\\chi_{DC}$', **args)
    ax.plot(SWP, [r[ID][1,2]for r in RES], label = '$\\chi_{BC}$', **args)
    ax.set_ylim([0.01,10**2]); ax.set_xlabel(swpvar); ax.set_title('cross-Kerr'); ax.set_ylabel('$\\chi$ (MHz)'); ax.legend(loc = 0)
    ax.set_yscale('log');   ax.set_ylim(0.1,100)    
    ax1.axhspan(5.5,6.5, alpha =0.4, color= 'b')
    ax1.axhline(0.5, alpha =0.4, color= 'b')
    ax1.axhspan(85,150, alpha =0.4, color= 'b')
    ax = ax2; 
    ax.plot(SWP, [r[ID][0,0] for r in RES], label = '$\\alpha_{D}$', **args)
    ax.plot(SWP, [r[ID][1,1] for r in RES], label = '$\\alpha_{B}$', **args)
    ax.plot(SWP, [r[ID][2,2] for r in RES], label = '$\\alpha_{C}$', **args)
    ax.set_ylim([100 +0.01,3.5*10**2]); ax.set_xlabel(swpvar); ax.set_title('Anharmonicity');ax.set_ylabel('$\\alpha$ (MHz)'); ax.legend(loc = 0)
    ax.set_yscale('linear')
    ax = ax3;  
    ax.plot(SWP, [r[9][:2]*10**-9 for r in RES],  **args)
    ax.plot(SWP, [r[8][:2]        for r in RES],  **{'lw':0,'marker':'x','ms':3})
    ax.set_xlabel(swpvar); ax.set_ylabel('Freq1 (GHz)'); ax.set_title('Freq.'); ax.legend(['D','B','C'], loc= 0); ax.grid(axis='y',color='gray', linestyle='-', linewidth=0.8, alpha =0.4)
    ax = ax3b;
    ax.plot(SWP, [r[11] for r in RES],  **args)
    ax.set_xlabel(swpvar); ax.set_ylabel('Q'); ax.legend(['D','B','C'], loc = 0)
    ax.set_yscale('log'); ax.set_title('Quality')
    
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
    ax.set_xlabel(swpvar); ax.set_ylabel('Participation'); ax.legend(loc = 0)
    ax = ax6
    ax.plot(SWP, [r[ID][2,0] for r in RES], label = '$P_{CH}$', **args)
    ax.plot(SWP, [r[ID][2,1] for r in RES], label = '$P_{CV}$', **args)
    ax.set_yscale('log'); ax.set_xlabel(swpvar); ax.set_ylabel('Participation'); ax.legend(loc = 0)
  
    ax   = plt.subplot(gs1[1,0]); ID =1;
    chiDC = np.array([r[ID][0,2] for r in RES])
    chiDB = np.array([r[ID][0,1] for r in RES])
    #print_color("chiDB/chiDC ratios:"); print  chiDB/chiDC
    ax.plot(SWP, chiDB/chiDC, **args); ax.locator_params(nbins=4); ax.grid(); ax.set_ylabel('$\\chi_{DB}/\\chi_{DC}$')
    ax.set_xlabel(swpvar);
    
    # plot the chis again     
    plt.close(3);   ID = 1;
#    fig, (ax7,ax8,ax9) = plt.subplots(3,1,sharex = True, num = 3, figsize=(6,7)) ; 
#    
#    ax7.plot(SWP, [r[ID][0,1]for r in RES], label = '$\\chi_{DB}$', c = 'b', **args); ax7.set_ylabel('$\\chi_{DB}$ (MHz)');
#    ax9.plot(SWP, [r[ID][0,2]for r in RES], label = '$\\chi_{DC}$', c = 'g', **args); ax9.set_ylabel('$\\chi_{DC}$ (MHz)');
#    ax8.plot(SWP, [r[ID][1,2]for r in RES], label = '$\\chi_{BC}$', c = 'r', **args); ax8.set_ylabel('$\\chi_{BC}$ (MHz)');
#    ax9.set_xlabel(swpvar); ax7.set_title('cross-Kerr');   
#    #ax7.axhspan(85,150,  alpha =0.4, color= 'b')
#    ax8.axhspan(5.5,6.5, alpha =0.4, color= 'b')
#    ax9.axhline(0.5,     alpha =0.4, color= 'b')
#    fig.tight_layout()
#    

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
    ax9.set_xlabel(swpvar); ax7.set_title('cross-Kerr');   
    #ax7.axhspan(85,150,  alpha =0.4, color= 'b')
    ax8.axhspan(5.5,6.5, alpha =0.4, color= 'b')
    ax9.axhline(0.5,     alpha =0.4, color= 'b')
    
    
#%%
if 1: # plot mesh
    fig = plt.figure(8);  fig.clf()
    tets = bba.get_convergences_max_tets()
    varsz  = bba.get_variable_vs(swpvar)
    Y = {}
    for key in tets.keys():
        Y[varsz[key]] = tets[key]
    y =  pd.Series(Y ) #.values(), index = varsz.values())
    y.plot(marker = '*', ms = 20)
    ax7t = ax7.twinx()
    ax7t.plot(y, marker = '*', ms = 10, c = 'g')
#%%
if 1:
    fig = plt.figure(21); fig.clf()
    tts = bba.get_convergences_Tets_vs_pass()
    for key, x in tts.iteritems():
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
