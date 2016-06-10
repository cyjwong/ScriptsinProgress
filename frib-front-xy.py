"""
Python input script for a Warp xy slice simulation of 
the FRIB front end.  For help and/or additional information contact:

     Steve Lund     lund@frib.msu.edu    (517) 908-7291
     Jonathan Wong  wong@nscl.msu.edu    (517) 908-7465 

For documentation see the Warp web-site:

     https://warp.lbl.gov

All code inputs are mks with the exception of particle kinetic energy (eV).   

To run this Warp script:
  Interactive (return to interpreter after executing):

    % python -i frib-front-xy.py 

  Non-interactive (exit interpreter after executing):

    % python frib-front-xy.py

This setup scripts for these simulations are large and have been broken 
apart into functional sections to aid maintainability.
 ** To be implemented **   

    frib-front-xy.py             :  Main setup/control for xy slice simulations 
    frib-front-xy-parms-name.py  :  Parameters for runs, name = ion type, U etc

    frib-front-lat.py       :  Lattice setup/description (3D, xy etc) 
    frib-front-lat-diag.py  :  Diagnostics for lattice     

    frib-front-xy-load.py   :  Initial beam specification 
    frib-front-xy-diag.py   :  Diagnostics 

    frib-front-env.py       :  Axisymmetric Envelope solver
    frib-front-env-diag.py  :  Diagnostics for envelope solver  

"""

# Load Warp and various script packages 
from warp import *                                # Warp code 
from warp.utils.errorcheck  import checksymmetry  # Check for input errors
#from warp.utils.runcounters import *             # Counter for parametric runs 

# Set informational labels included on all output cgm plots.   
top.pline2   = "xy Slice Simulation: FRIB Front End" 
top.pline1   = " "   # Add more info, if desired.  

# Parameters 
mr = 0.001  # milli-rad unit 

# Invoke setup routine for graphics and output files (THIS IS MANDATORY)
setup()

# Set runmaker - included in informational labels on output plots
top.runmaker = "Wong and Lund"

# Read in parameters describing beam to be simulated
#   - Use appropriate name input file for type of beam to be simulated. 
execfile("frib-front-xy-params-U.py") 


# --- Setup for variable particle weights in slice code. 
#     Notes: 
#      * All particles will carry individual weight set with v-velocity for full slice 
#          consistency using pid arrays to store the initial axial velocity. This is 
#          setup by hand in slice code for multi-species due to code structure.  
#      * Way done will be consistent with axial velocity spread.   
#      * Will be put in weight adjustment each step using beforeloadrho() 
#      * pid array elements hold particle properties and are consistently 
#          mirrored when particles scraped  
#      * nextpid() gets next pid array index 
#  
# Via Dave Grote: need to setup scaling by hand for multi-species 
#  top.wpid = nextpid()     # setup/allocate pid array on generate 
#
#  species.getw()  = s.w   gets weights of species  
#                          [equivalent to species.w] (top.wpid must be set to work)  
#  species.getweights() returns product of species.sw*species.w 
# 
#  beam.w * beam.sw = total weight = beam.getweights() 
# 

top.wpid = nextpid()       # pid index for variable weights: initialize on generate  
uzp0pid  = nextpid() - 1   # pid index for initial uz to scale weights: initialize on generate  

# Reference species
#  * Used to set gap accel and dipole bends 
#  * Set consistent with options in parameter specification file. 
#  * Corresponding reference particle energy to set bend strength set in lattice 
#    specification script 
#  * Logic sets:
#      m_ref = mass   ref particle [kg]    A_ref = m_ref/amu 
#      q_ref = charge ref particle [C]     Q_ref = q_ref/echarge 
#    These need not correspond to physical particle species (charge states 
#    can be non-integer).  
#

# --- reference defined by average target species 
if ref_mode == 0:
	A_ref = 0. 
	Q_ref = 0. 
	for ii in sp_target:
		s = sp[ii]
		A_ref += s.mass/amu 
		Q_ref += s.charge/echarge 
	A_ref = A_ref/len(sp_target) 
	Q_ref = Q_ref/len(sp_target) 
# --- reference defined by current weighted average target species
elif ref_mode == 2:
	A_ref = 0. 
	Q_ref = 0. 
	target_current = 0.
	for ii in sp_target:
		s = sp[ii]
		A_ref += s.mass*ibeam[ii]/amu 
		Q_ref += s.charge*ibeam[ii]/echarge 
		target_current += ibeam[ii]
	A_ref = A_ref/target_current 
	Q_ref = Q_ref/target_current 
# --- reference defined by named species 
elif ref_mode in sp.keys():
  A_ref = sp[ref_mode].mass/amu  
  Q_ref = sp[ref_mode].charge/echarge 
# --- reference A_ref and Q_ref input 
elif ref_mode == 3:
  pass  
# --- error trap 
else:
  raise Exception("Error: ref_mode not set properly")   

m_ref = A_ref*amu
q_ref = Q_ref*echarge


# Calculate vbeam and other species quantities for defined species 
derivqty()

# ---  Calculate and printout Q/M by species and store in a dictionary 
sp_qovm = {}
print("Species Charge to Mass Ratios:")
for ii in sort(sp.keys()):
  s = sp[ii]
  qovm = (s.charge/echarge)/(s.mass/amu)
  sp_qovm.update({ii:qovm})
  print("   Species: "+ii+" Q/A = %s"%qovm)

# ---  Calculate and printout injected rigidity [B rho] by species and store in a dictionary 
sp_brho = {}
print("Species Injected Rigidity:")
for ii in sort(sp.keys()):
  s = sp[ii]
  gamma = 1./sqrt(1.-(s.vbeam/clight)**2)
  brho  = gamma*s.mass*s.vbeam/s.charge   # Define rigidity with mean axial velocity 
  sp_brho.update({ii:brho})
  print("   Species: "+ii+" [B rho] = %s T-m"%brho)


# --- Set properties of initial species load 
#     * See notes for logic of formulas.  Note: Warp uses rms edge measures in loads. 
#         Care must be taken to consitently set various factors in coefficients.  

for s in sp.keys():
  ekin_i  = ekin[s]
  betab_i = sqrt(2.*jperev*ekin_i/sp[s].sm)/clight   # NR beta associated the KE 
  #
  rb_i    = sqrt(rx[s]*ry[s])          # take mean measure in case inhomogeneous
  #
  emitn_i = sqrt(emitnx[s]*emitny[s])  # take mean measure in case inhomogeneous  
  temp_i  = temp_p[s]  
  if init_ps_spec == "emitn":
    emitn_edge_i = 4.*emitn_i 
  elif init_ps_spec == "temp": 
    emitn_edge_i = betab_i*((2.*rb_i)/sqrt(2.))*sqrt(temp_i/ekin_i)
  else:
    raise Exception("Error: init_emit_spec not set properly") 
  #
  vthz_i = sqrt(2.*jperev*temp_z[s]/sp[s].sm)   
  # --- Energy and Current 
  sp[s].ekin   = ekin[s]             # kinetic energy of beam particle [eV]
  sp[s].vbeam  = 0.                  # beam axial velocity [m/sec] (set from ekin if 0) 
  sp[s].ibeam  = ibeam[s]            # beam current [amps] 
  # --- Perp PS areas and Par temp 
  sp[s].emitnx = emitn_edge_i        # beam x-emittance, rms edge [m-rad] 
  sp[s].emitny = emitn_edge_i        # beam y-emittance, rms edge [m-rad]
  sp[s].vthz   = vthz_i              # axial velocity spread [m/sec] 
  # --- centroid 
  sp[s].x0  = xc[s]  #   x0:   initial x-centroid xc = <x> [m]
  sp[s].y0  = yc[s]  #   y0:   initial y-centroid yc = <y> [m]
  sp[s].xp0 = xcp[s] #   xp0:  initial x-centroid angle xc' = <x'> = d<x>/ds [rad]
  sp[s].yp0 = ycp[s] #   yp0:  initial y-centroid angle yc' = <y'> = d<y>/ds [rad]
  # --- envelope 
  sp[s].a0   = rx[s]   #   a0:   initial x-envelope edge a = 2*sqrt(<(x-xc)^2>) [m]
  sp[s].b0   = ry[s]   #   b0:   initial y-envelope edge b = 2*sqrt(<(y-yc)^2>) [m]
  sp[s].ap0  = rxp[s]  #   ap0:  initial x-envelope angle ap = a' = d a/ds [rad]
  sp[s].bp0  = ryp[s]  #   bp0:  initial y-envelope angle bp = b' = d b/ds [rad]


# --- Max beam size measures to set distribution load properties 
r_x = maxnd(rx.values())
r_y = maxnd(ry.values())


#
# Setup Lattice  
#
#  Read in lattice description from auxillary script file 

execfile("frib-front-lat.py") 


# Define transverse simulation grid

sym_x = 1.
sym_y = 1.
if w3d.l4symtry:
  sym_x = 0.5
  sym_y = 0.5
elif w3d.l2symtry:
  sym_x = 0.5 

w3d.nx = int(sym_x*n_grid) 
w3d.ny = int(sym_y*n_grid)

# ---- Grid bounds 
#      Some bounds will be reset to zero by code on generation
#      if symmetry options are set.
l_grid = 2.*r_grid            # length edge of simulation grid [m]      
w3d.xmmax =  l_grid/2.        # x-grid max limit [m] 
w3d.xmmin = -l_grid/2.        # x-grid min limit [m] 
w3d.ymmax =  l_grid/2.        # y-grid max limit [m] 
w3d.ymmin = -l_grid/2.        # y-grid min limit [m] 

# --- grid increments to use before code generation in setup
dx = l_grid/float(n_grid)


# Particle loading
#
# Set simulation macro-particle number (top.npmax) by specifying the
#  number of macro-particles to load per xy grid cell (nppg) using an
# rms equivalent uniform density beam measure.  This number is (re)set
# consistently with the symmetry options.  
#   sp['name'].np = particle #
#   top.nplive = total number live particles = len(sp)*top.npmax at time of load 
# 

top.npmax = int(nppg*pi*(r_x*r_y)/dx**2*sym_x*sym_y) # max initial particles loaded (each species) 


top.zbeam = z_launch    # present z of simulation, to launch, reset consistently 

# --- random number options to use in loading 
w3d.xrandom  = "pseudo" # "digitrev"    # load x,y,z  with digitreverse random numbers 
w3d.vtrandom = "pseudo" # "digitrev"    # load vx, vy with digitreverse random numbers
w3d.vzrandom = "pseudo" # "digitrev"    # load vz     with digitreverse random numbers 
w3d.cylinder = true          # load a cylinder

# Particle moving
#

top.lrelativ   =  false    # turn off relativistic kinematics
top.relativity = 0         # turn off relativistic self-field correction
                           #   to account for approx diamagnetic B-field of beam

wxy.ds = ds                # ds for part adv [m] 
wxy.lvzchang = true        # Use iterative stepping, which is needed if
                           # the vz of the particles changes.
                           #  ... must change even in linear lattice 
                           #          for high-order energy conservation 
top.ibpush   = 2           # magnetic field particle push, 
                           #   0 - off, 1 - fast, 2 - accurate 


# Setup field solver using 2d multigrid field solver. 

w3d.boundxy = 0              # Neuman boundary conditions on edge of grid.
w3d.solvergeom = w3d.XYgeom  # fieldsolve type to 2d multigrid 

# --- Uncomment to turn off space-charge deposition for simulation of particles 
#     moving in applied field 
 
if single_particle: top.depos = "none"


# Turn on x-window plots, if desired; use winkill() to close interactively.  
#winon()


################################
# Particle simulation
################################

# Generate the xy PIC code.  
#   In the generate, particles are allocated and
#   loaded consistent with initial conditions and load parameters
#   set previously.  Particles are advanced with the step() command later
#   after various diagnostics are setup.

package("wxy"); generate()

# Make diagnostics for applied lattice fields 
execfile("frib-front-lat-diag.py") 

# Install conducting aperture on mesh for fieldsolves 
#for i in aperture:
  #installconductors(i,dfill=largepos) # will trigger "list has no attribute" error if aperture is used in installconductors

# Install conductors
for i in scraperlist:
  installconductors(i) # will trigger "list has no attribute" error if scraperlist is used as argument directrly

# Check that inputs are consistent with symmetries (errorcheck package function)
checksymmetry()

# Setup variable weight species needs for neutralization and acceleration 

# --- set initial weight factors consistent with neutralization factor 
#       w0 = initial weight (same for all particles in species) 
#       species.w = array of variable weight factors 
for sp_name in sp.keys():
  s = sp[sp_name]
  #       
  s.w0  = 1.-rho_neut_f(top.zbeam,sp_name)
  s.sw0    = s.sw       # save initial sw    (in case later changed)  
  s.vbeam0 = s.vbeam    # save initial vbeam (in case later changed)  

# --- save initial uzp for all species at once 
top.pgroup.pid[:,uzp0pid] = top.pgroup.uzp

# --- Adjustments implemented particle advance: 
#       Species weights and particle scrape aperture   
@callfrombeforeloadrho
def adjustments_before_rho():
  for sp_name in sp.keys():
    s = sp[sp_name]
    # --- Adjust species weights 
    s.w0 = 1.-rho_neut_f(top.zbeam,sp_name)
    s.w[:] = s.w0*s.pid[:,uzp0pid]/s.uzp
    #  --- Scraping scraping aperture:
    #       Set for efficient removal of out of bounds particles. 
    #top.prwall = aperture_r(top.zbeam)


# Modify ion distribution at launch point based on different assumptions on their birth 
#  in the ECR to reflect a target value of beam canonical angular momentum.  

execfile("frib-front-xy-load.py")

if birth_mode == 1 or birth_mode == 2:
	diag_plt_krot_launch() 
	diag_plt_krot_v()

# Carry out explicit fieldsolve with adjusted rho consistent with neutralization 
adjustments_before_rho()
loadrho() 
fieldsolve()

# Fix intitial history diagnostics to account for species weight changes
top.jhist = top.jhist-1   # needed to be minus 1 to reset save in right postion 
from warp.diagnostics.getzmom import *
zmmnt() 
savehist(0.) 

# Read in diagnostics for slice run 
execfile("frib-front-xy-diag.py") 

# Call initial history diagnostics accumulation to save initial load value 
#  if flagged for accumulation.   
diag_hist_hl()   

# Plot initial Brho by species 
plt_diag_bro(label = "Initial Rigidity by Species") 

# Plot initial beam potential profile 
         
# --- neutralized 
#       Using setup fields including neutralization 

diag_plt_phi_ax(label="Initial Neutralized Beam Potential at y,x = 0 b,r")
fma()

diag_plt_phi_ax(label="Initial Neutralized Beam Potential at y,x = 0 b,r",xmax=1.5*r_x)
fma()

# *** Code below did not work.  Appear to get same plots as unneutralized case.  Error?  Must be triggering weight 
#     adjustment if previous plots are correct.  
## --- unneutralized 
##       Must carry out fieldsolve with bare, unneutralized beam, then reset consistent with any neutralization 
#
## --- --- generate initial unneutralized self-field 
#for s in sp.values():       
#  s.w0 = 1.
#  s.w[:] = s.w0*s.pid[:,uzp0pid]/s.uzp
#
#fieldsolve() 
#
## --- --- unneutralized potential plots 
#diag_plt_phi_ax(label="Initial Unneutralized Beam Potential at y,x = 0 b,r") 
#fma()
#
#diag_plt_phi_ax(label="Initial Unneutralized Beam Potential at y,x = 0 b,r",xmax=1.5*r_x)
#fma()
#
## --- --- regenerate initial neutralized self-field 
#adjustments_before_rho() 
#loadrho() 
#fieldsolve() 
#

# Install diagnostic calls after simulation step
installafterstep(diag_calls)

# Step 0 diagnostics (if any) of the initial distribution loaded 
diag_calls() 

#
# Advance simulation specified steps 
#

#raise Exception("to here")
n_step = nint((z_adv-z_launch)/wxy.ds) + 1   # add one step in case of roundoff 
step_rigidities_plot = nint(((d5p1_zc + s4p2_zc)/2-z_launch)/wxy.ds)

step(step_rigidities_plot)

# Make plot of final Brho by species 
plt_diag_bro(label = "Final Rigidity by Species") # now takes place between 2nd S4 Solenoid and 1st D5 Dipole

step(n_step-step_rigidities_plot)


# Make additional history plots for final run if not already called 
#if not(top.it >= diag_hist_step.max()):
#  diag_hist()  # Add full arg list 

# Save restart dump of run.  By default the name of the dump file is
# top.runid (or script name if this is not set) with the step number (iii)
# and ".dump" appended to the name:
#       runidiii.pdb 
# To restart:
#   % python
#     >>> from warp import *
#     >>> restart("runidiii.dump") 
#
#dump() 



# Output data to auxillary file
output_data = false 
output_data_file = "frib-front-xy_data.txt" 

if output_data:
  fout = open(output_data_file,"a")
  #
  #fout.write("Run Number %s :\n"%irun) 
  fout.write(" Solenoid Excitations:\n")
  fout.write("  S41 = %s \n"%s4p1_str)
  fout.write("  S41 = %s \n"%s4p2_str)
  for ii in sp_target:
    s  = sp[ii] 
    js = s.js 
    rmsx    = top.hxrms[0,top.jhist,js]
    drmsxds = top.hxxpbar[0,top.jhist,js]/top.hxrms[0,top.jhist,js]
    # Final rms envelope width and angles 
    fout.write("   "+ii+": \n")
    fout.write("     sqrt(<x^2>) = %s mm, d sqrt(<x^2>)/ds = %s mr \n"%(rmsx/mm,drmsxds/mr))     
  #
  fout.close() 


# Print out timing statistics of run 
printtimers() 



# Multi-Species Envelope Model 

execfile("frib-front-env.py")

# Diagnostics to compare envelope model to Warp xy simulation 

execfile("frib-front-env-diag.py")


# Make sure that last plot is flushed from buffer
fma() 


