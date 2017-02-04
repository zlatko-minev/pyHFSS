Welcome to pyEPR!
===================
###Energy-Participation-Ratio (EPR) Design of Quantum Circuits. 
#####By Zlatko Minev & Zaki Leghtas.    

pyEPR is ... and is amazing because ...  (short blurb about pyEPR here). 

# Features
---------------------
TBA   


## Installation
-------------
If you are starting from scratch, follow the installation guide below to setup a Python 2.7 environment and fork this repository. To keep up to date, you may use SourceTree, a git-gui manager. 

**Recommended procedure.**   <br /> 

 1. Install a Python 2.7 environment.  
   * We recommend [Anaconda CE](https://www.continuum.io/downloads), and have tested this installation procedure with Anaconda v4.2 32-bit edition. Other environments, such as Python XY, and 64 bit work as well, but are not recommend. Install Anaconda in "C:\Anaconda2".
   * Set the Windows System PATH variable. In Control Panel, search for Environment Variables (in System), and open it. In the section System Variables, find the PATH environment variable and select it. Click Edit.  Place`C:\Anaconda2;C:\Anaconda2\Scripts;C:\Anaconda2\Library\bin;` at the beginning. If you have a previous Python installation this step is *very* important, especially to compile the qutip module. You may verity your path using the following command in the Command Prompt (terminal):
      ``` sh
      $ echo %PATH%
      ```  
    
 2. Install required packages, form a terminal
 ```sh 
 pip install pint
 conda install pandas
 ```
 3. Fork this pyEPR repository on GitHub with your GitHub account. You may clone the fork to your PC and manage it using the [SourceTree](https://www.sourcetreeapp.com/) git-gui manager.
 4. Add the pyEPR repository folder to your python search path.
 5. Edit pyEPR module `config.py`  to set your data-saving directory and other parameters of interest.   
 6. **ENJOY! **  :+1:  

#### Note for Mac/Linux.   
Follow the same instructions above. You shouldn't have to install mingw or modify distutils.cfg, since your distribution should come with gcc as the default compiler.    

####Optional package installation
You may also choose to install the optional qutip package for some advanced analysis. 

Optionally, you may install [Qutip](http://qutip.org/), used in handling quantum objects. Follow the instruction on their website, or the rest of this bullet. Here, we summarize trick points in this. First, you need to install a C compiler, since Qutip uses Cython. If you dont have VS9, gcc, or mingw installed, the following works:   
```sh
	pip install -i https://pypi.anaconda.org/carlkl/simple mingwpy
```
Let anaconda know to use this compiler by creating the file `C:\Anaconda2\Lib\distutils\distutils.cfg` with the following content   
```
    [build]
    compiler = mingw32    
    [build_ext]
    compiler = mingw32
```   
Next, let's install qutip. You can choose to use conda intall or pip install, or pull from the git directly  as done here: 
```sh
        conda install git
        pip install git+https://github.com/qutip/qutip.git
```

## pyEPR Project Setup in HFSS
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


## Troubleshooting
---------------------

###### COM Error on opening HFSS 
Check the project and design file names carefully. Make sure that the filepath doesnt have apostrophies or other bad charecters, such as in C:\\Minev's PC\\my:Project.  Check that HFSS hasn't popped up an error dialogue, such as "File locked." Manually open HFSS and the file.    

###### COM error on calculation of expression
Either HFSS popped an error dialogue, froze up, or you mistyped the name of something.    

###### HFSS refuses to close
If your script terminates improperly, this can happen. pyHFSS tries to catch termination events and handle them. Your safety should be guaranteed however, if you call `hfss.release()` when you have finished. Use the Taskmanager (Activity Monitor on MAC) to kill HFSS if you want.   
