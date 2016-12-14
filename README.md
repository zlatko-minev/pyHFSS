# pyEPR: Energy Participation Ratio Approach to Quantum Circuits   
Authors (alphabetical): Zaki Leghtas & Zlatko Minev.   
Short about pyEPR here. Key features here.   


## Installation
-------------
If you have a python 2.7 enviornment setup, with qutip installed, fork the repository and keep up to date with SourceTee. If you are starting from scratch, follow the guide below.

##### Installation on windows from scratch.
Recommended procedure.   
1. Install [Anaconda CE](https://www.continuum.io/downloads) (Python XY would work too, but we will follow the CE install). This is your python environment.  64 bit or 32 bit should work---32 bit tends to lead to less installation conflicts. If there is a previous Python install, either delete it, or in your System Path variable set the Anaconda search path at the top. Install in C:\Anaconda2 <br /> 
2. Install required packages, from a command prompt (terminal):   
```sh
        pip install pint 
        conda install pandas
```   
3. Install [Qutip](http://qutip.org/), which is used in dealing with quantum objects. Follow the instruction on their website, or the rest of this bullet. Qutip is only required to do more complicated quantum analysis on the data from HFSS.   
First, you need to install a C compiler, since Qutip uses Cython. If you dont have VS9, gcc, or mingw installed, the following works:   
```sh
        pip install -i https://pypi.anaconda.org/carlkl/simple mingwpy
```   
This is the compiler anaconda will use. Let anaconda know to use it by creating the file C:\Anaconda2\Lib\distutils\distutils.cfg with the following content   
```
    [build]
    compiler = mingw32    
    [build_ext]
    compiler = mingw32
```   
Next, install qutip, you can choose to use conda intall or pip install, or pull from the git directly   
```sh
        conda install git
        pip install git+https://github.com/qutip/qutip.git
 ```
 4. Use [SourceTree](https://www.sourcetreeapp.com/) or git to fork this repository.   
 5. Add the repository folder to your python search path.    
 6. Edit the config.py  to set the data saving directory correctly. :+1:     
  
 
##### Installation on Mac/Linux from scratch.   
Follow the windows instructions, but dont install mingw and make distutils.cfg, your distribution should come with gcc as the default compiler.    

## HFSS Setup for EPR 
-------------
Eigenmode setup tips   
1. Geometry & boundary condition (BC) definitions.   
  1.1 Define junction rectangles & boundary condition. E.g. create a rectangle for the Josephson junction, name it "junc1" (for instance). Define its length, which youll put into the script, and give it a lumped RLC BC with a local variable "Lj1" (for instance) for the inductance. This name will be input into the script.   Recomended junc rectangle size is 50 x 100 um, and much smaller will work as well.    
  1.2 Over this rectangle draw a model polyline to define current flow direction. It should spans the length of the junc. rectangle, and be defined ny the same variables to used to define junc1---so that they move together if the geometry is altered.    
  1.2 Tip: Define an object coordinate axis before you make the polyline, which is based on the junc rect, junc1. Next, define the polyline in the Object CS ( make sure to give the CS a proper name, such as junc1_cs). Next, name the polyline sinsibly; e.g., "junc1_line."   
4. Meshing.   
 4.1 Lighly mesh the pads and junc rectangles with roughly 4 initial tets across the smallest dimension.    
5. Simulation setup    
 5.1 Advisable, not necessary to used mixed order solutions.    

## Features
---------------------
TBA   

## Troubleshooting
---------------------

###### COM Error on opening HFSS 
Check the project and design file names carefully. Make sure that the filepath doesnt have apostrophies or other bad charecters, such as in C:\\Minev's PC\\my:Project.  Check that HFSS hasn't popped up an error dialogue, such as "File locked." Manually open HFSS and the file.    

###### COM error on calculation of expression
Either HFSS popped an error dialogue, froze up, or you mistyped the name of something.    

###### HFSS refuses to close
If your script terminates improperly, this can happen. pyHFSS tries to catch termination events and handle them. Your safety should be guaranteed however, if you call `hfss.release()` when you have finished. Use the Taskmanager (Activity Monitor on MAC) to kill HFSS if you want.   
