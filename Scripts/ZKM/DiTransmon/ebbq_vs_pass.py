#%% @author: Zlatko
import sys;  IMP_PATH = r'C:\\Users\\rslqulab\\Desktop\\zkm\\pyHFSS\\';
if ~(IMP_PATH in sys.path): sys.path.insert(0,IMP_PATH);

import pandas as pd, matplotlib.pyplot as plt, numpy as np;
import hfss, bbq, bbqNumericalDiagonalization
from hfss import CalcObject, ureg, load_HFSS_project
from bbq  import eBBQ_Pmj_to_H_params, print_color, print_matrix


if 1:    
    proj_name    = r'mesh_analysis_project' 
    project_path = 'C:\\Users\rslqulab\Desktop\Lysander\\'
    app, desktop, project = load_HFSS_project(proj_name, project_path)
    design       = project.get_design("20pass_7k_seed- 11. SHANTANU FAB SM22 Col 1 Row 9 Measured [May 6 2016]") 
    
    junc_rects    = ['juncV',     'juncH'] 
    junc_lines    = ['juncV_line','juncH_line'] 
    junc_LJ_names = ['LJ1', 'LJ2'];
    junc_lens     = [0.0001]*2                                   
    
    setup_names   = design.get_setup_names()
    setup         = design.get_setup(setup_names[0])
    setup.delta_f = '0.00001'                                       # set max delta f convergence 
    ### Set min # passes = 1 & min conv passes = 1
    # wrap in function  / class eventually 
    #delte all data?
    optimtrc_nms  = design._optimetrics.GetSetupNames()
    optimtrc_nms  = [optimtrc_nms[i] for i in  range(optimtrc_nms.Count)]
    #assert len(optimtrc_nms) == 1, "THis was made for optimetric sweeps which have only 1 optimetric swep "
    optimtrc_name = optimtrc_nms[0]  #'scale_factor_tiny'          # optimetrics sweep  -- doesnt have to be 
    
    for max_passes in range(1, 20, 2):
        print_color("Running optimetric sweep '%s' up to %d passes" % (optimtrc_name, max_passes))
        setup.passes  = max_passes
        design._optimetrics.SolveSetup(optimtrc_name)        # wait for the optimetrics to finish  --> would be nice to add a timer here and call asynch
        
        bbp = bbq.Bbq(project, design, append_analysis=False)
        bbp.do_eBBQ(junc_rect=junc_rects, junc_lines = junc_lines,  junc_len = junc_lens, junc_LJ_var_name = junc_LJ_names, save_mesh_stats=True)
        bba = bbp.bbq_analysis 