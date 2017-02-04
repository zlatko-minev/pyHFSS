Welcome to pyEPR!
===================
###Energy-Participation-Ratio (EPR) Design of Quantum Circuits. 
#####By Zlatko Minev & Zaki Leghtas.    

About pyEPR. TBA.

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
      ```sh
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
#### Eigenmode Design --- How to set up junctions
You may find an advised workflow and some setup tips here.

 1. Define circuit geometry & electromagnetic boundary condition (BC).   
   1. Junction rectangles and BC: Create a rectangle for each Josephson junction and give it a good name; e.g., `junc1`. We recommend 50 x 100 um rectangle for a simple simulation, although orders of magnitude smaller rectangles work as well. Note the length of this junction, you will supply it to pyEPR. Assign a `Lumped RLC` BC on this rectangle surface, with an inductance value given by a local variable, `Lj1` for instance. The name of this variable will also be supplied to the pyEPR. 
   2. Over each junction rectangle draw a model `polyline` to define give a sense of the junction current-flow direction. This line should spans the length of the full junction rectangle. Define it using an object coordinate system on the junction rectangle (so that they move together when the geometry is altered). The name of this line will be supplied to the pyEPR module.
 2. Meshing.   
   1. Lightly mesh thin-film metal BC. Lightly mesh the junction rectangles, so that each has roughly 4 tets across its smallest dimension in an initial mesh seed.    
 3. Simulation setup    
   1. We recommend `mixed order` solutions.    


## Troubleshooting
---------------------

###### COM Error on opening HFSS 
Check the project and design file names carefully. Make sure that the file-path doesn't have apostrophes or other bad characters, such as in C:\\Minev's PC\\my:Project.  Check that HFSS hasn't popped up an error dialogue, such as "File locked." Manually open HFSS and the file.    

###### COM error on calculation of expression
Either HFSS popped an error dialogue, froze up, or you miss-typed the name of something.    

###### HFSS refuses to close
If your script terminates improperly, this can happen. pyHFSS tries to catch termination events and handle them. Your safety should be guaranteed however, if you call `hfss.release()` when you have finished. Use the Task-manager (Activity Monitor on MAC) to kill HFSS if you want.   
