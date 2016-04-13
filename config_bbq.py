# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 14:21:45 2015

@author: Zaki
"""

root_dir = r'C:\Users\rslqulab\Desktop\zkm\2016_qm_jumps\DiTransmon_Asymetric\\'

# seams:
# ref: http://arxiv.org/pdf/1509.01119.pdf
gseam = 1.0e3 # per Ohm meter: seam conductance

# surfaces:
# ref: http://arxiv.org/pdf/1509.01854.pdf
th = 3e-9 # dirt thickness on dielectric
eps_r = 10 # dielectric constant of dirt
tan_delta_surf = 1e-3 # tan(delta) for surfaces 

# bulk:
# ref: http://arxiv.org/pdf/1509.01854.pdf
tan_delta_sapp = 1e-6 # tan(delta) for bulk surface
epsi = 10 # dielectric