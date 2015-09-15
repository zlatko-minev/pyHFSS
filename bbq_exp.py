from hfss import get_active_project
from bbq import Bbq

project = get_active_project()
design = project.get_design("Qubit")

bbq = Bbq(project, design)
bbq.do_bbq('$LJ_Qubit')