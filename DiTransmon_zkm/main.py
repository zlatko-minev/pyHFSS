#%%
import sys;  IMP_PATH = r'C:\\Users\\rslqulab\\Desktop\\zkm\\pyHFSS\\';
if ~(IMP_PATH in sys.path): sys.path.append(IMP_PATH);
    
import hfss;   import bbqNumericalDiagonalization; import pandas as pd
from hfss import CalcObject, ureg
from hfss import load_HFSS_project
import bbq, matplotlib.pyplot as plt, numpy as np;  from bbq import print_color, divide_diagonal_by_2, print_matrix
from IPython.display import display 

    
if 1:    
    proj_name    = r'2016_03_28_Di_Transmon from Antonio' 
    project_path = 'C:\\Users\\rslqulab\\Desktop\\zkm\\2016_qm_jumps\\DiTransmon_Asymetric\\'
    app, desktop, project = load_HFSS_project(proj_name, project_path)
    design       = project.get_active_design() #get_design("CoaxCav")
    
    bbp = bbq.Bbq(project, design, append_analysis=False, calculate_H=True)

if 1:
    junc_rects    = ['juncV',     'juncH'] 
    junc_lines    = ['juncV_line','juncH_line'] 
    junc_LJ_names = ['LJ1', 'LJ2'];
    junc_lens     = [0.0001]*2                                                       # this can soon be replaced by intgrating over junc_lines 
    bbp.do_eBBQ(calculate_H = False,  plot_fig = False,
               Pj_from_current= True, junc_rect=junc_rects, junc_lines = junc_lines,  junc_len = junc_lens, junc_LJ_var_name = junc_LJ_names )
    sol = bbp.PJ_multi_sol

#%%
def eBBQ_ND(freqs, PJ, Om, EJ, LJs, SIGN, cos_trunc = 6, fock_trunc  = 7):
    ''' numerical diagonalizaiton for energy BBQ
        fzpfs: reduced zpf  ( in units of \phi_0
    '''    
    from bbqNumericalDiagonalization import bbq_hmt, make_dispersive, fqr
    
    fzpfs = np.zeros(PJ.T.shape)
    for junc in xrange(fzpfs.shape[0]):
        for mode in xrange(fzpfs.shape[1]):
            fzpfs[junc, mode] = np.sqrt(PJ[mode,junc] * Om[mode,mode] /  EJ[junc,junc] ) #*0.001
    fzpfs = fzpfs * SIGN.T
    
    H     = bbq_hmt(freqs*10**9, LJs.values.astype(np.float), fqr*fzpfs, cos_trunc, fock_trunc)
    f1s, CHI_ND, fzpfs, f0s  = make_dispersive(H, fock_trunc, fzpfs, freqs)  # f0s = freqs
    CHI_ND= -1*CHI_ND *1E-6;
    return f1s, CHI_ND, fzpfs, f0s; 
    
def eBBQ_participation2_H_params(s, cos_trunc = None, fock_trunc = None):
    '''   
    returns the CHIs as MHz with anharmonicity alpha as the diagonal  (with - sign)
        f1: qubit dressed freq
        f0: qubit linear freq (eigenmode) 
        and an overcomplete set of matrcieis
    '''
    import  scipy;    Planck  = scipy.constants.Planck
    f0s        = s['freq'].values    
    Qs         = s.loc[:,'modeQ']
    LJs        = s.loc[0,s.keys().str.contains('LJs')] # LJ in nH
    EJs        = (bbq.fluxQ**2/LJs/Planck*10**-9).astype(np.float)        # EJs in GHz
    PJ_Jsu     = s.loc[:,s.keys().str.contains('pJ')]  # EPR from Jsurf avg
    PJ_Jsu_sum = PJ_Jsu.apply(sum, axis = 1)           # sum of participations as calculated by avg surf current 
    PJ_glb_sum = (s['U_E'] - s['U_H'])/(2*s['U_E'])    # sum of participations as calculated by global UH and UE  
    diff       = (PJ_Jsu_sum-PJ_glb_sum)/PJ_glb_sum*100# debug
    if 1:  # Renormalize: to sum to PJ_glb_sum; so that PJs.apply(sum, axis = 1) - PJ_glb_sum =0
           #TODO: figure out the systematic   # print '% diff b/w Jsurf_avg & global Pj:'; display(diff)
        PJs = PJ_Jsu.divide(PJ_Jsu_sum, axis=0).mul(PJ_glb_sum,axis=0)
    else: PJs = PJ_Jsu
    SIGN  = s.loc[:,s.keys().str.contains('sign_')]
    PJ    = np.mat(PJs.values)
    Om    = np.mat(np.diagflat(f0s)) 
    EJ    = np.mat(np.diagflat(EJs.values))
    CHI_O1= Om * PJ * EJ.I * PJ.T * Om * 1000       # MHz
    CHI_O1= divide_diagonal_by_2(CHI_O1)            # Make the diagonals alpha 
    f1s   = f0s - np.diag(CHI_O1)                   # 1st order PT expect freq to be dressed down by alpha 
    if cos_trunc is not None:
        f1s, CHI_ND, fzpfs, f0s = eBBQ_ND(f0s, PJ, Om, EJ, LJs, SIGN, cos_trunc = cos_trunc, fock_trunc = fock_trunc)                
    else: CHI_ND, fzpfs = None, None
    return CHI_O1, CHI_ND, PJ, Om, EJ, diff, LJs, SIGN, f0s, f1s, fzpfs, Qs
#%%
if 1:
    variation = '0'; s           = sol[variation];  
    cos_trunc = 6;   fock_trunc  = 7;
    CHI_O1, CHI_ND, PJ, Om, EJ, diff, LJs, SIGN, f0s, f1s, fzpfs, Qs = \
        eBBQ_participation2_H_params(s, cos_trunc = cos_trunc, fock_trunc = fock_trunc)
    print '\nPJ=\t(renorm.)';        print_matrix(PJ*SIGN, frmt = "{:7.4f}")
    #print '\nCHI_O1=\t PT. [alpha diag]'; print_matrix(CHI_O1,append_row ="MHz" )
    print '\nf0={:6.2f} {:7.2f} {:7.2f} GHz'.format(*f0s)
    print '\nCHI_ND=\t PJ O(%d) [alpha diag]'%(cos_trunc); print_matrix(CHI_ND, append_row ="MHz")
    print '\nf1={:6.2f} {:7.2f} {:7.2f} GHz'.format(*(f1s*1E-9))   
    varz = bbp.get_variables(variation=variation)     
    print pd.Series({ key:varz[key] for key in ['_join_w','_join_h','_padV_width', '_padV_height','_padH_width', '_padH_height','_scaleV','_scaleH'] })
#%%==============================================================================
#     Plot results for sweep
#==============================================================================
if 1:
    swpvar='join_h'    
    RES = []; SWP = [];
    for key, s in sol.iteritems():
        varz  = bbp.get_variables(variation=key)
        SWP  += [ ureg.Quantity(varz['_'+swpvar]).magnitude ]  
        RES  += [ eBBQ_participation2_H_params(s, cos_trunc = cos_trunc, fock_trunc = fock_trunc) ]
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