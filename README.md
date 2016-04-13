Energy participation BBQ
======
Zlatko's eBBQ setup tips for Eigenmode simulation 
1. Geometry & B.C. definitions:
1A.  Define junction rectangles & boundary condition. E.g. rectangle name: juncV. Define the lenght of the junction and junction Lj by a local variable; e.g., junc_len and LJ_V.   It is good to keep the junc rectangle on the order of 50 x 100 um. 
1B. Define junction line to define current flow direction (polyline). Make it a dmodel objecta and make sure it spans the length of the junction, so use the same variables to define it. Tip: Define an Object coordinate axis before you make the polyline, which is based on the junc rect (e.g., juncV); then define the polyline in the Object CS ( make sure to give the CS a proper name juncV_cs). Now, give the polyline a good name: E.g.,juncV_line.
4. Meshing tips:
4A. Lighly mesh the pads and junc rectangles with roughly 4 initial tets across the smallest dimension. 
5. Simulation setup 
5A. Advisable, not necessary to used mixed order solutions. 

pyHFSS
======

HFSS scripting interface in python

Create a Design
---------------

```python
from hfss import get_active_project
proj = get_active_project()
design = proj.new_dm_design("Test")
```
    
Or Get an Existing Design
-------------------------

```python
from hfss import get_active_design
design = get_active_design()
```

Creating Variables
------------------

```python
bx = design.set_variable("Box_X", "3mm")
by = design.set_variable("Box_Y", "6mm")
bz = design.set_variable("Box_Z", "1mm")
```
    

3D Modeler
----------

```python
modeler = design.modeler
modeler.draw_box_center([0,0,0], [bx, by, bz], material="silicon")
```

Setup Analysis
--------------

```python
setup = design.create_dm_setup(freq_ghz=5)
sweep = setup.insert_sweep(4, 10, count=1000)
setup.analyze()
freqs, (S12, Y11) = sweep.get_network_data("S12,Y11")
```

Fields Calculator
-----------------

```python
fields = setup.get_fields()
Mag_E_Sq = fields.Mag_E ** 2
Surface_E = Mag_E_Sq.integrate_surf("Object Name")
print Surface_E.evaluate()
```


Keyword Arguments for Drawing Commands
--------------------------------------

  - name: str
  - nonmodel: bool
  - color: (int, int, int) each in [0...255]
  - transparency: float in [0, 1]
  - material: str (matching existing material name)

HFSS refuses to close
---------------------

If your script terminates improperly, this can happen. pyHFSS tries to
catch termination events and handle them. Your safety should be
guaranteed however, if you call `hfss.release()` when you have finished

Requires
---------------------

pint     - pip install pint
pandas   - conda install pandas

do_BBQ:
---------------------
    calculate_H:  
        True: 1 junction method of Pj calculation based on U_H-U_E global. 
        
    Pj_from_current:
        Multi-junction calculation of energy participation ratio matrix based on <I_J>. Current is integrated average of J_surf by default: (zkm 3/29/16)
        Will calculate the Pj matrix for the selected modes for the given junctions junc_rect array & length of juuncs
        
        junc_rect = ['junc_rect1', 'junc_rect2'] name of junc rectangles to integrate H over
        junc_len = [0.0001]   specify in SI units; i.e., meters
        junc_LJ_var_name = ['LJ1', 'LJ2']
        pJ_method = 'J_surf_mag'   - currently only 1 implemented - takes the avg. Jsurf over the rect. Make sure you have seeded lots of tets here. i recommend starting with 4 across smallest dimension.

        Assumptions:
            Low dissipation (high-Q). 
            Right now, we assume that there are no lumped capcitors to simply calculations. Not required. 
            We assume that there are only lumped inductors, so that U_tot = U_E+U_H+U_L    and U_C =0, so that U_tot = 2*U_E;
        Results in:
            self.PJ_multi_sol - a Pandas DataFrame of all the information
    
    Other parameters:
        seams = ['seam1', 'seam2']  (seams needs to be a list of strings)
        variations = ['0', '1']
