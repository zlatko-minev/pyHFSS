from hfss import get_active_project
import bbq
import matplotlib.pyplot as plt
import numpy as np
plt.close('all')

project = get_active_project()
design = project.get_design("Storage_Readout")
#
bbq_exp = bbq.Bbq(project, design, append_analysis=True, calculate_H=True)

#bbq_exp.do_bbq('$LJ_Qubit')
#bbq_exp.bbq_analysis.plot_Hparams(variable_name='$LJ_Qubit')

data_filename = r'Y:\Data\PumpingCats\HFSS\Analyzed\2015-10-08_3_Coax_cav\Storage_Readout\Storage_Readout_20151013_181350.hdf5'
bbq_analysis_exp = bbq.BbqAnalysis(data_filename)
fig, ax = bbq_analysis_exp.plot_Hparams(variable_name='_$r_pad')
bbq_analysis_exp.plot_Hparams(variable_name='_$LJ_Qubit')