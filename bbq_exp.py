from __future__ import print_function    # Python 2.7 and 3 compatibility  
from hfss import get_active_project
import bbq
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import *
plt.close('all')

project = get_active_project()
design = project.get_design("Design")

bbq_exp = bbq.Bbq(project, design, append_analysis=False, calculate_H=True)

bbq_exp.do_bbq("LJ", surface=True, seams=["seam1"])

#bbq_exp.bbq_analysis.plot_Hparams()


#bbq_exp.get_Qseam_sweep('seam1', 0, '0', 'seamz',np.linspace(0,35,36), 'mm')

#for k in range(39):
#    bbq_exp.get_Qseam('seam1', 0, str(k))

#ba=bbq_exp.bbq_analysis
#bbq_exp.bbq_analysis.print_Hparams(variation='0',modes=[0,1,2])

#data_filename = r'Y:\Data\PumpingCats\HFSS\Analyzed\2015-10-02_closerPads\Dump1\Dump1_20151014_204956.hdf5'
#bbq_analysis_exp = bbq.BbqAnalysis(data_filename)
#fig, ax = bbq_analysis_exp.plot_Hparams(variable_name='_$r_pad')
#bbq_analysis_exp.plot_Hparams(variable_name='_$LJ_Dump')