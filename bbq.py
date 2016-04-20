from hfss import *
from hfss import CalcObject
import time, os, shutil, matplotlib.pyplot as plt, numpy as np, pandas as pd
from stat import S_ISREG, ST_CTIME, ST_MODE
from pandas import HDFStore, Series, DataFrame
from scipy.constants import *; from scipy.constants import hbar, e as e_el, epsilon_0, pi;  # not sure what else ened sto be imported, idellay we should get rid of all *
from config_bbq      import root_dir, gseam, th, eps_r, tan_delta_surf, tan_delta_sapp
from pint import UnitRegistry; 

ureg = UnitRegistry(system='mks')

#import h5py

#==============================================================================
# Utility functions and difinitions
#==============================================================================
fluxQ = hbar / (2*e_el)

def fact(n):
    if n <= 1:
        return 1
    return n * fact(n-1)

def nck(n, k):
    return fact(n)/(fact(k)*fact(n-k))
    
import warnings
def deprecated(func):
    """This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emmitted
    when the function is used."""
    def newFunc(*args, **kwargs):
        warnings.simplefilter('always', DeprecationWarning) #turn off filter 
        warnings.warn("Call to deprecated function {}.".format(func.__name__), category=DeprecationWarning, stacklevel=2)
        warnings.simplefilter('default', DeprecationWarning) #reset filter
        return func(*args, **kwargs)
    newFunc.__name__ = func.__name__
    newFunc.__doc__ = func.__doc__
    newFunc.__dict__.update(func.__dict__)
    return newFunc


def print_matrix(M, frmt = "{:7.2f}", append_row = ""):
    M = np.mat(M)
    for row in np.array(M.tolist()):
        print ' ',
        for chi in row:
            print frmt.format(chi),
        print append_row+"\n",
        
def divide_diagonal_by_2(CHI0):
    CHI = CHI0.copy();
    CHI[np.diag_indices_from(CHI)] /= 2
    return CHI;
    
def print_NoNewLine(text):
    print(text),

def print_color(text, style = 0, fg=24, bg = 43, newline = True):
    '''style 0..8;   fg  30..38;  bg  40..48'''
    format = ';'.join([str(style), str(fg), str(bg)])
    s = '\x1b[%sm %s \x1b[0m' % (format, text)
    if newline: print s 
    else: print s,


#==============================================================================
# Main compuation class & interface with HFSS
#==============================================================================
class Bbq(object):
    """ 
    This class defines a BBQ object which calculates and saves
    Hamiltonian parameters from an HFSS simulation
    """
    
    def __init__(self, project, design, verbose=True, append_analysis=False, setup_name = None):
        '''  calculate_H is the single-jucntion method using UH-Ue '''
        self.project = project
        self.design  = design
        self.setup   = design.get_setup(name=setup_name)
        self.fields  = self.setup.get_fields()
        self.nmodes  = int(self.setup.n_modes)
        self.listvariations   = design._solutions.ListVariations(str(self.setup.solution_name))
        self.nominalvariation = design.get_nominal_variation()
        self.nvariations      = np.size(self.listvariations)
        self.solutions        = self.setup.get_solutions()
        self.verbose          = verbose
        self.append_analysis  = append_analysis
        self.hfss_variables   = {}                             # container for eBBQ list of varibles  
        self.sols             = {}                             # container for eBBQ solutions; could make a Panel
        self.meta_data        = {}                             # container for eBBQ metadata
        
        self.setup_data()
        if self.verbose: print '       # Modes: ' + str(self.nmodes), '\n  # Variations: ' + str(self.nvariations)
        
        self.get_latest_h5()
        if self.latest_h5_path is not None and self.append_analysis:
            latest_bbq_analysis = BbqAnalysis(self.latest_h5_path)
            if self.verbose: print 'Varied variables and values : ', latest_bbq_analysis.get_swept_variables(), \
                                   'Variations : ', latest_bbq_analysis.variations

    def get_latest_h5(self):
        dirpath = self.data_dir
        
        entries1 = (os.path.join(dirpath, fn) for fn in os.listdir(dirpath))     # get all entries in the directory w/ stats
        entries2 = ((os.stat(path), path) for path in entries1)
        entries3 = ((stat[ST_CTIME], path)                                       # leave only regular files, insert creation date
                   for stat, path in entries2 if S_ISREG(stat[ST_MODE]) and path[-4:]=='hdf5')
        #NOTE: on Windows `ST_CTIME` is a creation date but on Unix it could be something else
        #NOTE: use `ST_MTIME` to sort by a modification date
        
        paths_sorted = []
        for cdate, path in sorted(entries3):
            paths_sorted.append(path)
            #print time.ctime(cdate), os.path.basename(path)
        if len(paths_sorted) > 0:
            self.latest_h5_path = paths_sorted[-1]
            if self.verbose: print 'This simulations has been analyzed, latest data in ' + self.latest_h5_path
        else:
            self.latest_h5_path = None
            if self.verbose: print 'This simulation has never been analyzed'
        
    def setup_data(self):
        data_dir = root_dir + '/' + self.project.name + '/' + self.design.name
        if self.verbose: print data_dir
        if not os.path.isdir(data_dir):
            os.makedirs(data_dir)
        self.data_dir = data_dir
        self.data_filename = self.data_dir + '/' + self.design.name + '_' + time.strftime('%Y%m%d_%H%M%S', time.localtime()) + '.hdf5'
        if self.verbose: print "Data will be saved in " + str(data_dir)
        
    @deprecated
    def calc_p_j(self, modes=None, variation=None):
        '''
        Calculates the p_j for all the modes. 
        Requires a calculator expression called P_J.
        '''
        lv = self.get_lv(variation)
        if modes is None:
            modes = range(self.nmodes)

        pjs = {}
        for ii, m in enumerate(modes):
            print 'Calculating p_j for mode ' + str(m) + ' (' + str(ii) + '/' + str(np.size(modes)-1) + ')'
            self.solutions.set_mode(m+1, 0)
            self.fields = self.setup.get_fields()
            P_J = self.fields.P_J
            pjs['pj_'+str(m)] = P_J.evaluate(lv=lv)
        self.pjs = pjs
        if self.verbose: print pjs
        return pjs

    def get_p_j(self, mode):
        pj = {}
        pj_val = (self.U_E-self.U_H)/(2*self.U_E)
        pj['pj_'+str(mode)] = np.abs(pj_val)
        print '    p_j_' + str(mode) + ' = ' + str(pj_val)
        return pj
    
    def get_freqs_bare(self, variation):
        #str(self.get_lv(variation))
        freqs_bare_vals = []
        freqs_bare_dict = {}
        freqs, kappa_over_2pis = self.solutions.eigenmodes(self.get_lv_EM(variation))
        for m in range(self.nmodes):
            freqs_bare_dict['freq_bare_'+str(m)] = 1e9*freqs[m]
            freqs_bare_vals.append(1e9*freqs[m])
            if kappa_over_2pis is not None:
                freqs_bare_dict['Q_'+str(m)] = freqs[m]/kappa_over_2pis[m]
        self.freqs_bare = freqs_bare_dict
        self.freqs_bare_vals = freqs_bare_vals
        return freqs_bare_dict, freqs_bare_vals
        
        
    def get_lv(self, variation):
        ''' variation is a string #; e.g., '0'
            returns array of var names and var values '''
        if variation is None:
            lv = self.nominalvariation
            lv = self.parse_listvariations(lv)
        else:
            lv = self.listvariations[ ureg(variation) ]
            lv = self.parse_listvariations(lv)
        return lv
    
    def get_lv_EM(self, variation):
        if variation is None:
            lv = self.nominalvariation
            #lv = self.parse_listvariations_EM(lv)
        else:
            lv = self.listvariations[ ureg(variation) ]
            #lv = self.parse_listvariations_EM(lv)
        return str(lv)
    
    def parse_listvariations_EM(self,lv):
        lv = str(lv)
        lv = lv.replace("=",":=,")
        lv = lv.replace(' ',',')
        lv = lv.replace("'","")
        lv = lv.split(",")
        return lv
        
    def parse_listvariations(self,lv):
        lv = str(lv)
        lv = lv.replace("=",":=,")
        lv = lv.replace(' ',',')
        lv = lv.replace("'","")
        lv = lv.split(",")
        return lv
        
    def get_variables(self,variation=None):
        lv = self.get_lv(variation)
        variables={}
        for ii in range(len(lv)/2):
            variables['_'+lv[2*ii][:-2]]=lv[2*ii+1]
        self.variables = variables
        return variables
    
#    @deprecated # TODO: delete this 
#    def save_data(self, data, variation):
#        group = self.h5file.create_group(variation)
#        for name, val in data.items():
#            group[name] = val
   
    def get_Qseam(self, seam, mode, variation):
        '''
        caculate the contribution to Q of a seam, by integrating the current in
        the seam with finite conductance: set in the config file
        ref: http://arxiv.org/pdf/1509.01119.pdf
        '''
        lv = self.get_lv(variation)
        Qseam = {}
        print 'Calculating Qseam_'+ seam +' for mode ' + str(mode) + ' (' + str(mode) + '/' + str(self.nmodes-1) + ')'
        j_2_norm = self.fields.Vector_Jsurf.norm_2() # overestimating the loss by taking norm2 of j, rather than jperp**2
        int_j_2 = j_2_norm.integrate_line(seam)
        int_j_2_val = int_j_2.evaluate(lv=lv, phase=90)
        yseam = int_j_2_val/self.U_H/self.omega
        Qseam['Qseam_'+seam+'_'+str(mode)] = gseam/yseam
        print 'Qseam_' + seam + '_' + str(mode) + str(' = ') + str(gseam/yseam)
        return Series(Qseam)

    def get_Qseam_sweep(self, seam, mode, variation, variable, values, unit, pltresult=True):
        # values = ['5mm','6mm','7mm']
        # ref: http://arxiv.org/pdf/1509.01119.pdf
        
        self.solutions.set_mode(mode+1, 0)
        self.fields = self.setup.get_fields()
        freqs_bare_dict, freqs_bare_vals = self.get_freqs_bare(variation)
        self.omega = 2*np.pi*freqs_bare_vals[mode]
        print variation
        print type(variation)
        print ureg(variation)
        self.U_H = self.calc_U_H(variation)
        lv = self.get_lv(variation)
        Qseamsweep = []
        print 'Calculating Qseam_'+ seam +' for mode ' + str(mode) + ' (' + str(mode) + '/' + str(self.nmodes-1) + ')'
        for value in values:
            self.design.set_variable(variable, str(value)+unit)
            
            j_2_norm = self.fields.Vector_Jsurf.norm_2() # overestimating the loss by taking norm2 of j, rather than jperp**2
            int_j_2 = j_2_norm.integrate_line(seam)
            int_j_2_val = int_j_2.evaluate(lv=lv, phase=90)
            yseam = int_j_2_val/self.U_H/self.omega
            Qseamsweep.append(gseam/yseam)
#        Qseamsweep['Qseam_sweep_'+seam+'_'+str(mode)] = gseam/yseam
            #Cprint 'Qseam_' + seam + '_' + str(mode) + str(' = ') + str(gseam/yseam)
        if pltresult:
            fig, ax = plt.subplots()
            ax.plot(values,Qseamsweep)
            ax.set_yscale('log')
            ax.set_xlabel(variable+' ('+unit+')')
            ax.set_ylabel('Q'+'_'+seam)
        return Qseamsweep

    def get_Qdielectric(self, dielectric, mode, variation):
        Qdielectric = {}
        print 'Calculating Qdielectric_'+ dielectric +' for mode ' + str(mode) + ' (' + str(mode) + '/' + str(self.nmodes-1) + ')'
        
        U_dielectric = self.calc_U_E(variation, volume=dielectric)
        p_dielectric = U_dielectric/self.U_E
        Qdielectric['Qdielectric_'+dielectric+'_'+str(mode)] = 1/(p_dielectric*tan_delta_sapp)
        print 'p_dielectric'+'_'+dielectric+'_'+str(mode)+' = ' + str(p_dielectric)
        return Series(Qdielectric)

    def get_Qsurface(self, mode, variation):
        '''
        caculate the contribution to Q of a dieletric layer of dirt on all surfaces
        set the dirt thickness and loss tangent in the config file
        ref: http://arxiv.org/pdf/1509.01854.pdf
        '''
        lv = self.get_lv(variation)
        Qsurf = {}
        print 'Calculating Qsurface for mode ' + str(mode) + ' (' + str(mode) + '/' + str(self.nmodes-1) + ')'
#        A = self.fields.Mag_E**2
#        A = A.integrate_vol(name='AllObjects')
#        U_surf = A.evaluate(lv=lv)
        calcobject=CalcObject([],self.setup)
        vecE=calcobject.getQty("E")
        A=vecE
        B=vecE.conj()
        A=A.dot(B)
        A=A.real()
        A=A.integrate_surf(name='AllObjects')
        U_surf = A.evaluate(lv=lv)
        U_surf *= th*epsilon_0*eps_r
        p_surf = U_surf/self.U_E
        Qsurf['Qsurf_'+str(mode)] = 1/(p_surf*tan_delta_surf)
        print 'p_surf'+'_'+str(mode)+' = ' + str(p_surf)
        return Series(Qsurf)
    
    def get_Hparams(self, freqs, pjs, lj):
        Hparams = {}
        fzpfs = []
        
        # calculate Kerr and fzpf
        for m in self.modes:
            omega = 2*pi*freqs[m]
            ej = fluxQ**2/lj
            pj = pjs['pj_'+str(m)]
            fzpf = np.sqrt(pj*hbar*omega/ej)
            fzpfs.append(fzpf)
            Hparams['fzpf_'+str(m)] = fzpf
            alpha = 2*ej/fact(4)*nck(4,2)*(fzpf**4)/hbar
            Hparams['alpha_'+str(m)] = alpha
            Hparams['freq_'+str(m)]=(omega-alpha)/2/pi

        # calculate chi
        for m in self.modes:
            for n in self.modes:
                if n<m:
                    chi_mn = ej/hbar*(fzpfs[m]*fzpfs[n])**2
                    Hparams['chi_'+str(m)+'_'+str(n)] = chi_mn

        return Hparams
       
    def calc_U_E(self, variation, volume=None):
        ''' This is 2 * the peak electric energy.(since we do not divide by 2, and use the peak phasors) '''
        lv = self.get_lv(variation)
        if volume is None:
            volume = 'AllObjects'
        else:
            pass
        calcobject=CalcObject([],self.setup)
        vecE=calcobject.getQty("E")
        A=vecE.times_eps()
        B=vecE.conj()
        A=A.dot(B)
        A=A.real()
        A=A.integrate_vol(name=volume)
        return A.evaluate(lv=lv)
        
    def calc_U_H(self, variation, volume=None):
        lv = self.get_lv(variation)
        if volume is None:
            volume = 'AllObjects'
        else:
            pass
        calcobject=CalcObject([],self.setup)
        vecH=calcobject.getQty("H")
        A=vecH.times_mu()
        B=vecH.conj()
        A=A.dot(B)
        A=A.real()
        A=A.integrate_vol(name=volume)
        return A.evaluate(lv=lv)
        
    def calc_current(self, fields, line ):
        '''Function to calculate Current based on line. Not in use 
            line = integration line between plates - name 
        '''
        self.design.Clear_Field_Clac_Stack()
        comp = fields.Vector_H
        exp  = comp.integrate_line_tangent(line)
        I    = exp.evaluate(phase = 90)
        self.design.Clear_Field_Clac_Stack()
        return I
        
    def calc_avg_current_J_surf_mag(self, variation, junc_rect, junc_len):
        ''' Peak current I_max for mdoe J in junction J  
            The avg. is over the surface of the junction. I.e., spatial. '''
        lv   = self.get_lv(variation)
        calc = CalcObject([],self.setup)
        calc = calc.getQty("Jsurf").mag().integrate_surf(name = junc_rect)
        I    = calc.evaluate(lv=lv) / junc_len #phase = 90
        #self.design.Clear_Field_Clac_Stack()
        return  I
    
    def calc_line_current(self, variation, junc_line_name):
        lv   = self.get_lv(variation)
        calc = CalcObject([],self.setup)
        calc = calc.getQty("H").imag().integrate_line_tangent(name = junc_line_name)
        #self.design.Clear_Field_Clac_Stack()
        return calc.evaluate(lv=lv)
        
    def calc_Pjs_from_I_for_mode(self,variation, U_H,U_E, LJs, junc_rects,junc_lens, method = 'J_surf_mag' , 
                                 freq = None, calc_sign = None):
        ''' Expected that you have specified the mode before calling this 
            Expected to precalc U_H and U_E for mode, will retunr pandas series object 
                junc_rect = ['junc_rect1', 'junc_rect2'] name of junc rectangles to integrate H over
                junc_len  = [0.0001]   specify in SI units; i.e., meters
                LJs       = [8e-09, 8e-09] SI units
                calc_sign = ['junc_line1', 'junc_line2']    used to define sign of ZPF
            Potential errors:  If you dont have a line or rect by the right name you will prob get an erorr o the type:
                com_error: (-2147352567, 'Exception occurred.', (0, None, None, None, 0, -2147024365), None)
        '''
        dat = {} 
        for i, junc_rect in enumerate(junc_rects):
            print_NoNewLine('     ' + junc_rect)
            if method is 'J_surf_mag':
                I_peak = self.calc_avg_current_J_surf_mag(variation, junc_rect, junc_lens[i])                         
            else: 
                print 'Not yet implemented.'
            if LJs is None: print_color(' -----> ERROR: Why is LJs passed as None!?')
            #dat['I_'  +junc_rect] = I_peak # stores the phase information as well
            dat['pJ_' +junc_rect] = LJs[i] * I_peak**2 / (2*U_E) 
            if calc_sign is not None:
                Idum = self.calc_line_current(variation, calc_sign[i])
                dat['sign_'+junc_rect] = +1 if Idum > 0 else -1
                print   '  %+.5f' %(dat['pJ_' +junc_rect] * dat['sign_'+junc_rect] )
            else: print '  %0.5f' %(dat['pJ_' +junc_rect])
        return pd.Series(dat) 
        
    def do_eBBQ(self, variations= None, plot_fig  = False, modes      = None,
               Pj_from_current  = True, junc_rect = [],    junc_lines = None,  junc_len = [],  junc_LJ_var_name = [],    
               dielectrics      = None, seams     = None,  surface    = False, 
               calc_Hamiltonian = False,pJ_method =  'J_surf_mag'):
        """               
            Pj_from_current:
                Multi-junction calculation of energy participation ratio matrix based on <I_J>. Current is integrated average of J_surf by default: (zkm 3/29/16)
                Will calculate the Pj matrix for the selected modes for the given junctions junc_rect array & length of juuncs
                
                junc_rect  = ['junc_rect1', 'junc_rect2'] name of junc rectangles to integrate H over
                junc_lines = ['junc_line1', 'junc_line2'] used to define the current flow direction, arbitrary, doesnt really matter that much, just need a line there
                junc_len   = [0.0001]                     lenght of junc = lenght of junc_line #TODO: could now get rid of this and use the line     [specify in SI units; i.e., meters]
                junc_LJ_var_name = ['LJ1', 'LJ2']
                pJ_method  = 'J_surf_mag'   - takes the avg. Jsurf over the rect. Make sure you have seeded lots of tets here. i recommend starting with 4 across smallest dimension.

                Assumptions:
                    Low dissipation (high-Q). 
                    Right now, we assume that there are no lumped capcitors to simply calculations. Not required. 
                    We assume that there are only lumped inductors, so that U_tot = U_E+U_H+U_L    and U_C =0, so that U_tot = 2*U_E;
            
            Other parameters:
                seams = ['seam1', 'seam2']  (seams needs to be a list of strings)
                variations = ['0', '1']
            
            A variation is a combination of project/design variables in an optimetric sweep
        """

        self.Pj_from_current = Pj_from_current;  meta_data = {};  assert(type(junc_LJ_var_name) == list), "Please pass junc_LJ_var_name as a list "
        if Pj_from_current        :  print_color(' Setup: ' + self.setup.name); self.PJ_multi_sol = {} # this is where the result will go             
        if seams       is not None:  self.seams       = seams;       meta_data['seams']       = seams;    
        if dielectrics is not None:  self.dielectrics = dielectrics; meta_data['dielectrics'] = dielectrics;
        if variations      is None:  variations = (['-1'] if self.listvariations == (u'',)  else [str(i) for i in range(self.nvariations)] )
        if modes           is None:  modes = range(self.nmodes)
        if self.latest_h5_path is not None and self.append_analysis:shutil.copyfile(self.latest_h5_path, self.data_filename);
        self.h5file     = hdf = pd.HDFStore(self.data_filename); 
        self.variations = variations;  self.modes = modes; self.njunc = len(junc_rect)
        meta_data['junc_rect'] = junc_rect; meta_data['junc_lines'] = junc_lines; meta_data['junc_len'] = junc_len; meta_data['junc_LJ_var_name'] = junc_LJ_var_name; meta_data['pJ_method'] = pJ_method;

        for ii, variation in enumerate(variations):
            print_color( 'variation : ' + variation + ' / ' + str(self.nvariations-1), bg = 44, newline = False )
            self.lv = self.get_lv(variation)
            if (variation+'/hfss_variables') in hdf.keys() and self.append_analysis: print_NoNewLine('  previously analyzed ...\n');  \
                continue;    

            print_NoNewLine( ' NOT analyzed\n' );  time.sleep(0.5)
            hdf[variation+'/hfss_variables'] = self.hfss_variables[variation] = varz \
                                             = pd.Series(self.get_variables(variation=variation))
            freqs_bare_dict, freqs_bare_vals = self.get_freqs_bare(variation)   # get bare freqs from HFSS

            self.pjs={}; var_sol_accum = [] 
            for mode in modes:
                sol = Series({'freq' : freqs_bare_vals[mode]*10**-9, 'modeQ' : freqs_bare_dict['Q_'+str(mode)] })
                self.omega  = 2*np.pi*freqs_bare_vals[mode] # this should really be passed as argument  to the functions rather than a property of the calss I would say 
                print ' Mode  \x1b[0;30;46m ' +  str(mode) + ' \x1b[0m / ' + str(self.nmodes-1)+'  calculating:'
                self.solutions.set_mode(mode+1, 0)
                self.fields = self.setup.get_fields()

                print_NoNewLine('   U_H ...');     sol['U_H'] = self.U_H = self.calc_U_H(variation)
                print_NoNewLine('   U_E');         sol['U_E'] = self.U_E = self.calc_U_E(variation)
                print(  "   =>   U_L = %.3f%%" %( (self.U_E - self.U_H )/(2*self.U_E)) )
                
                if self.Pj_from_current:
                    self.LJs    = [ ureg.Quantity(varz['_'+LJvar_nm]).to_base_units().magnitude  for LJvar_nm in junc_LJ_var_name]
                    meta_data['LJs'] = dict(zip(junc_LJ_var_name, self.LJs))
                    print '   I -> p_{mJ} ...'
                    sol_PJ = self.calc_Pjs_from_I_for_mode(variation, self.U_H, self.U_E, self.LJs, junc_rect, junc_len, 
                                                                 method = pJ_method, freq = freqs_bare_vals[mode]*10**-9,
                                                                 calc_sign = junc_lines)
                    sol = sol.append(sol_PJ)
                
                if self.njunc == 1:             # Single-junction method using global U_H and U_E; 
                    assert(type(junc_LJ_var_name) == list and len(junc_LJ_var_name) == 1), "Please pass junc_LJ_var_name as array of 1 element for a single junction; e.g., junc_LJ_var_name = ['junc1']" 
                    #lj  = 1E-3*ureg.Quantity(varz['_'+junc_LJ_var_name]).to_base_units().magnitude                        
                    sol['pj1'] = self.get_p_j(mode)
                    self.pjs.update(sol['pj1'])        # convinience function for single junction case
                    
                if seams is not None:           # get seam Q
                    for seam in seams: sol = sol.append(self.get_Qseam(seam,mode,variation))

                if dielectrics is not None:     # get Q dielectric      
                    for dielectric in dielectrics: sol = sol.append(self.get_Qdielectric(dielectric, mode, variation))
                               
                if surface is True:             # get Q surface                              
                    sol = sol.append( self.get_Qsurface(mode, variation) )
                    
                var_sol_accum +=[sol]
            
            #TODO: add metadata to the Dataframe & save it
            #      such as what are the junc_rect names and Lj values etc.  (e.g., http://stackoverflow.com/questions/29129095/save-additional-attributes-in-pandas-dataframe/29130146#29130146)
            hdf[variation+'/eBBQ_solution']  = self.sols[variation]  \
                                             = pd.DataFrame(var_sol_accum, index = modes)            
            hdf[variation+'/meta_data']      = self.meta_data[variation]  \
                                             = Series(meta_data)
            if calc_Hamiltonian:  raise('Not implemented'); #for 1 junct: self.get_Hparams(freqs_bare_vals, self.pjs, lj))
            
        self.h5file.close()
        self.bbq_analysis = BbqAnalysis(self.data_filename, variations=self.variations)
#TODO: to be implemented below
#        if plot_fig:
#            self.bbq_analysis.plot_Hparams(modes=self.modes)
#            self.bbq_analysis.print_Hparams(modes=self.modes)
        return
    

def eBBQ_ND(freqs, PJ, Om, EJ, LJs, SIGN, cos_trunc = 6, fock_trunc  = 7):
    ''' numerical diagonalizaiton for energy BBQ
        fzpfs: reduced zpf  ( in units of \phi_0
    '''    
    from bbqNumericalDiagonalization import bbq_hmt, make_dispersive, fqr
    
    fzpfs = np.zeros(PJ.T.shape)
    for junc in xrange(fzpfs.shape[0]):
        for mode in xrange(fzpfs.shape[1]):
            fzpfs[junc, mode] = np.sqrt(PJ[mode,junc] * Om[mode,mode] /  EJ[junc,junc] ) #*0.001
    fzpfs = fzpfs * SIGN.T
    
    H     = bbq_hmt(freqs*10**9, LJs.astype(np.float), fqr*fzpfs, cos_trunc, fock_trunc)
    f1s, CHI_ND, fzpfs, f0s  = make_dispersive(H, fock_trunc, fzpfs, freqs)  # f0s = freqs
    CHI_ND= -1*CHI_ND *1E-6;
    return f1s, CHI_ND, fzpfs, f0s;
    
def eBBQ_Pjm_to_H_params(s, meta_data, cos_trunc = None, fock_trunc = None):
    '''   
    returns the CHIs as MHz with anharmonicity alpha as the diagonal  (with - sign)
        f1: qubit dressed freq
        f0: qubit linear freq (eigenmode) 
        and an overcomplete set of matrcieis
        ask zkm for info.
    '''
    import  scipy;    Planck  = scipy.constants.Planck
    f0s        = np.array( s['freq'] )
    Qs         = s['modeQ']
    LJs        = np.array(meta_data['LJs'].values())                     # LJ in H
    EJs        = (fluxQ**2/LJs/Planck*10**-9).astype(np.float)        # EJs in GHz
    PJ_Jsu     = s.loc[:,s.keys().str.contains('pJ')]  # EPR from Jsurf avg
    PJ_Jsu_sum = PJ_Jsu.apply(sum, axis = 1)           # sum of participations as calculated by avg surf current 
    PJ_glb_sum = (s['U_E'] - s['U_H'])/(2*s['U_E'])    # sum of participations as calculated by global UH and UE  
    diff       = (PJ_Jsu_sum-PJ_glb_sum)/PJ_glb_sum*100# debug
    if 1:  # Renormalize
        PJs = PJ_Jsu.divide(PJ_Jsu_sum, axis=0).mul(PJ_glb_sum,axis=0)
    else: PJs = PJ_Jsu
    SIGN  = s.loc[:,s.keys().str.contains('sign_')]
    PJ    = np.mat(PJs.values)
    Om    = np.mat(np.diagflat(f0s)) 
    EJ    = np.mat(np.diagflat(EJs))
    CHI_O1= Om * PJ * EJ.I * PJ.T * Om * 1000       # MHz
    CHI_O1= divide_diagonal_by_2(CHI_O1)            # Make the diagonals alpha 
    f1s   = f0s - np.diag(CHI_O1)                   # 1st order PT expect freq to be dressed down by alpha 
    if cos_trunc is not None:
        import bbqNumericalDiagonalization
        from bbqNumericalDiagonalization import eBBQ_ND;
        f1s, CHI_ND, fzpfs, f0s = eBBQ_ND(f0s, PJ, Om, EJ, LJs, SIGN, cos_trunc = cos_trunc, fock_trunc = fock_trunc)                
    else: CHI_ND, fzpfs = None, None
    return CHI_O1, CHI_ND, PJ, Om, EJ, diff, LJs, SIGN, f0s, f1s, fzpfs, Qs
    # the return could be made clener, or dictionary 


#%%    
class BbqAnalysis(object):
    ''' defines an analysis object which loads and plots data from a h5 file
    This data is obtained using e.g bbq.do_bbq

    ''' 
    def __init__(self, data_filename, variations=None):
        #raise('not implemented')
        self.data_filename = data_filename
        with HDFStore(data_filename, mode = 'r') as hdf:  # = h5py.File(data_filename, 'r')        
            # i think we should open & close the file here, i dont see why we need to keep it open & keep accessing it. It is small in memeory, just load it into the RAM.
            # all the data will be stored in 3 objects.        
            if variations is None:
                import re
                variations = []
                for key in hdf.keys():
                    if 'hfss_variables' in key:
                        variations += re.findall(r'\b\d+\b', key)
            self.variations     = variations
            self.hfss_variables = {}
            self.sols           = {}
            self.meta_data      = {}
            for variation in variations:
                self.hfss_variables[variation] = hdf[variation+'/hfss_variables']
                self.sols[variation]           = hdf[variation+'/eBBQ_solution']  
                self.meta_data[variation]      = hdf[variation+'/meta_data']
            self.nmodes         = self.sols[variations[0]].shape[0] 
            
    @deprecated  
    def get_swept_variables(self):
        #TODO: needs to be updated to new standard; currently borken
        swept_variables_names = []
        swept_variables_values = []
        for name in self.h5data[self.variations[0]].keys():
            if '_'==name[0]: # design variables all start with _
                variables = []
                for variation in self.variations:
                    variables.append(self.h5data[variation][name].value)
                if len(set(variables))>1:
                    swept_variables_names.append(name)
                    swept_variables_values.append(list(set(variables)))
            else:
                pass
        return swept_variables_names, swept_variables_values
    
    @deprecated
    def get_variable_variations(self, variablename):
        variables = []
        for variation in self.variations:
            variables.append(self.h5data[variation][variablename].value)
        return np.asarray(variables)
    
    @deprecated
    def get_float_units(self, variable_name, variation='0'):
        variable_value = self.h5data[variation][variable_name].value
        n = 1
        try:
            float(variable_value)
            return float(variable_value), ''
        except ValueError:
            while True:
                try:
                    float(variable_value[:-n])
                    return float(variable_value[:-n]), variable_value[len(variable_value)-n:]
                except:
                    n+=1
    @deprecated
    def print_Hparams(self, variation=None, modes=None):
        #TODO: needs to be updated to new standard; currently borken
        if modes==None:
            modes = range(self.nmodes)
        else:
            pass
        if variation == None:
            variation = self.variations[-1]
        else:
            pass
        swept_variables_names, swept_variables_values = self.get_swept_variables()

        for vname in swept_variables_names:
            print vname + ' = ' + self.h5data[variation][vname].value
        for ii, m in enumerate(modes):        
            freq_m = 'freq_'+str(m)
            Kerr_m = 'alpha_'+str(m)
            Q_m = 'Q_'+str(m)
            if freq_m not in self.h5data[variation].keys():
                freq_m = 'freq_bare_'+str(m)
            else:
                pass
            if Kerr_m in self.h5data[variation].keys():
                print Kerr_m + ' = ' +str(self.h5data[variation][Kerr_m].value/2/pi/1e6) + ' MHz'
            else:
                pass
 
            print freq_m +' = ' + str(self.h5data[variation][freq_m].value/1e9) + ' GHz'   
            if Q_m in self.h5data[variation].keys():             
                print Q_m  + ' = ' + str(self.h5data[variation][Q_m].value)
            else:
                pass
            
            for n in modes[0:ii]:
                chi_m_n = 'chi_'+str(m)+'_'+str(n)
                if chi_m_n in self.h5data[variation].keys():
                    print chi_m_n + ' = ' + str(self.h5data[variation][chi_m_n].value/2/pi/1e6) + ' MHz'
       
    @deprecated
    def plot_Hparams(self, variable_name=None, modes=None):
        #TODO: needs to be updated to new standard; currently borken
        fig, ax = plt.subplots(2,2, figsize=(24,10))

        if variable_name == None:
            xaxis = self.variations
        else:
            xaxis = []
            for variation in self.variations:
                xaxis.append(self.get_float_units(variable_name, variation)[0])

        if modes==None:
            modes = range(self.nmodes)
        else:
            pass
    
        for ii, m in enumerate(modes):        
            freq_m = 'freq_'+str(m)
            Kerr_m = 'alpha_'+str(m)
            Q_m = 'Q_'+str(m)
            Qsurf_m = 'Qsurf_'+str(m)
            
            if freq_m not in self.h5data[self.variations[0]].keys():
                freq_m = 'freq_bare_'+str(m)
            else:
                pass
            if Kerr_m in self.h5data[self.variations[0]].keys():
                ax[0][1].plot(xaxis, self.get_variable_variations(Kerr_m)/2/pi/1e6, 'o', label = str(m))
            else:
                pass
    
            ax[0][0].plot(xaxis, self.get_variable_variations(freq_m)/1e9, 'o', label=str(m))     
            
            if Q_m in self.h5data[self.variations[0]].keys():             
                ax[1][1].plot(xaxis, self.get_variable_variations(Q_m), 'o', label = Q_m)
            else:
                pass
            
            if Qsurf_m in self.h5data[self.variations[0]].keys():             
                ax[1][1].plot(xaxis, self.get_variable_variations(Qsurf_m), 'o', label = Qsurf_m)
            else:
                pass            
            
            if 'seams' in self.h5data[self.variations[0]].keys():
                for seam in self.h5data[self.variations[0]]['seams'].value:
                    Qseam_m = 'Qseam_'+seam+'_'+str(m)
                    if Qseam_m in self.h5data[self.variations[0]].keys():             
                        ax[1][1].plot(xaxis, self.get_variable_variations(Qseam_m), 'o', label = Qseam_m)
                    else:
                        pass
                    
            if 'dielectrics' in self.h5data[self.variations[0]].keys():
                for dielectric in self.h5data[self.variations[0]]['dielectrics'].value:
                    Qdielectric_m = 'Qdielectric_'+dielectric+'_'+str(m)
                    if Qdielectric_m in self.h5data[self.variations[0]].keys():             
                        ax[1][1].plot(xaxis, self.get_variable_variations(Qdielectric_m), 'o', label = Qdielectric_m)
                    else:
                        pass            
            
            for n in modes[0:ii]:
                chi_m_n = 'chi_'+str(m)+'_'+str(n)
                if chi_m_n in self.h5data[self.variations[0]].keys():
                    ax[1][0].plot(xaxis, self.get_variable_variations(chi_m_n)/2/pi/1e6, 'o', label=str(m)+','+str(n))
        
        ax[0][0].legend()
        ax[0][0].set_ylabel('freq (GHz)')
        
        ax[0][1].legend()
        ax[0][1].set_ylabel('Kerr/2pi (MHz)')
        ax[0][1].set_yscale('log')
        
        ax[1][0].legend()
        ax[1][0].set_ylabel('Chi/2pi (MHz)')
        ax[1][0].set_yscale('log')
        
        ax[1][1].legend()
        ax[1][1].set_ylabel('Q')
        ax[1][1].set_yscale('log')
        
        if variable_name == None:
            swept_variables_names, swept_variables_values = self.get_swept_variables()
            xticks = []
            for variation in xaxis:
                xtick = ''
                for name in swept_variables_names:
                    xtick += name[1:] + ' = ' + self.h5data[variation][name].value + '\n'
                xticks.append(str(xtick))
            ax[1][0].set_xticks([int(v) for v in xaxis])
            ax[1][0].set_xticklabels(xticks, rotation='vertical')
            ax[1][1].set_xticks([int(v) for v in xaxis])
            ax[1][1].set_xticklabels(xticks, rotation='vertical')
            
            ax[0][0].set_xticklabels([])
            ax[0][1].set_xticklabels([])
        else:
            xlabel = variable_name + ' (' + self.get_float_units(variable_name, self.variations[0])[1] + ')'
            ax[1][0].set_xlabel(xlabel)
            ax[1][1].set_xlabel(xlabel)

        fig.subplots_adjust(bottom=0.3)
        fig.suptitle(self.data_filename)
        fig.savefig(self.data_filename[:-5]+'.jpg')

        return fig, ax
    
            
            
#        for variable in swept_variables_names:
#            fig1 = plt.subplots()
#            ax1 = fig1.add_subplot(221)
#            ax.scatter()
#        return