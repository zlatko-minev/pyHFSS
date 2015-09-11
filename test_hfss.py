from hfss import get_active_project
import numpy as np

project = get_active_project()
design = project.get_design("Qubit1")


freqs = solutions.eigenmodes()
nmodes = int(setup.n_modes)
nmodes = 1


    