# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 14:58:53 2016

@author: Zlatko
"""

import bbq, matplotlib.pyplot as plt
from bbq import BbqAnalysis
file_name = u'C:\\Users\\rslqulab\\Desktop\\zkm\\2016_qm_jumps\\DiTransmon_Asymetric\\\\/pin_position_sweep(perfect conductor)_4-20-16/10. SHANTANU FAB 1 [April13 2016]/10. SHANTANU FAB 1 [April13 2016]_20160427_143612.hdf5'
swp_var   = 'pin_shift'
bba = BbqAnalysis(file_name)
Qs  = bba.get_Qs(swp_var=swp_var)

args = {'lw':0,'marker':'o','ms':5}
Qs.plot(**args); plt.legend(['D', 'B','G' ], loc = 0)
ax = plt.gca();
ax.set_xlabel(swp_var);
ax.set_title('Qs'); 
ax.set_ylabel('Q'); 
ax.set_yscale('log'); 