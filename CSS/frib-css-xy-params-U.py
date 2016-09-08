# Parameters for FRIB front-end simulations using the x-y slice model 
# described in frib-front-xy.py 
# Notes:
#  * include inputs describing beam AND numerical (mesh, advance) parameters  

 
#########################
# Begin Inputs 
#########################

# Define species
#   Operate = Operating species targeted to be transported (e.g. U)
#   Support = Support   species used in ECR (e.g. O in U operation) 
# 
#     Syntax of later use:
#
#       for sp in species_name:
#         sp.ppzx()
#
#     Dictionaries also setup. 

# Setup U operation with O support gas 

# --- Define charge states 
Operate_charge_states = [33,34,25,26,27,28,29,30,31,32,35,36,37,38,39,40]
Support_charge_states = [1,2,3,4]

# --- Create species consistent with charge states 
Operate_species = [Species(type=Uranium,charge_state=i,name='U%d'%i) for i in Operate_charge_states]
Support_species = [Species(type=Oxygen, charge_state=i,name='O%d'%i) for i in Support_charge_states]

# --- Make charge state dictionary 
Operate_charge_states_dict = {'U%d'%i: i for i in Operate_charge_states} 
Support_charge_states_dict = {'O%d'%i: i for i in Support_charge_states} 

charge_states = {} 
charge_states.update(Operate_charge_states_dict)
charge_states.update(Support_charge_states_dict) 

# --- Count species 
Operate_ns = len(Operate_species) 
Support_ns = len(Support_species) 
nspecies = Operate_ns + Support_ns 

# --- Make abbreviated dictionary of species "sp" for cleaner use in later diagnostics 
#       sp.keys()
sp_Operate = {Operate_species[i].name:Operate_species[i] for i in range(Operate_ns)}
sp_Support = {Support_species[i].name:Support_species[i] for i in range(Support_ns)}

sp = {}
sp.update(sp_Operate)
sp.update(sp_Support)

# --- Define "target" species list for later use in diagnostics 
#       Target species are those desired to transport downstream. 
sp_target = ['U33','U34']






#
# Reference particle specification 
# 
# The reference particle dictates how the dipole bend strength(s) are set to bend the required angles. 
#   * Set consistently with kinetic energy gain of reference charge state when accelerated through 
#     the source extractor and downstream ES gap to achieve the desired KeV/u energy to feed the 
#     downstream RFQ. 
#   * Works by specifying mass and charge of reference particle: A_ref and Q_ref .
#   * Values of A_ref and Q_ref below may be reset depending on set option.      
#   
#
# ref_mode = 0             : Use average species mass and charge state from target species list sp_target  
#            "Species Name": Use specified species name: e.g. "U33" 
#                              (typically should be in target species list sp_target, but can be any)
#            1             : Use spefified values of A_ref and Q_ref
#               A_ref = mass   reference species [amu]  (reset if ref_mode not equal 1)  
#               Q_ref = charge reference species [e]    (reset if ref_mode not equal 1)
#            2             : Use current-weighted average of species mass and charge state from target species list sp_target

ref_mode = 1 
A_ref = 238.02891    # reset when ref_mode not 1 
Q_ref = 34.0         # reset when ref_mode not 1 




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
elif ref_mode == 1:
  pass  
# --- error trap 
else:
  raise Exception("Error: ref_mode not set properly")   

m_ref = A_ref*amu
q_ref = Q_ref*echarge






# --- Assign species colors to help distinguish on plots: 
#        
#        O        species: blue 
#        High   U species: cyan 
#        Target U species: magenta/green 
#        Low    U species: red 

sp['O1'].color = "blue"
sp['O2'].color = "blue"
sp['O3'].color = "blue"
sp['O4'].color = "blue"

sp['U25'].color = "red"
sp['U26'].color = "red"
sp['U27'].color = "red"
sp['U28'].color = "red"
sp['U29'].color = "red"
sp['U30'].color = "red"
sp['U31'].color = "red"
sp['U32'].color = "red"

sp['U33'].color = "magenta"
sp['U34'].color = "green"

sp['U35'].color = "cyan"
sp['U36'].color = "cyan"
sp['U37'].color = "cyan"
sp['U38'].color = "cyan"
sp['U39'].color = "cyan"
sp['U40'].color = "cyan"

# --- Species unneutralized electric currents
#      * Set in a species named dictionary ibeam. 
#      * Data from ECR measurement data.    

Operate_ibeam = array([0.210,0.205,0.035,0.051,0.068,0.088,0.115,0.150,0.175,0.192,0.178,0.142,0.110,0.072,0.043,0.031])*mA 
Support_ibeam = array([0.300,0.300,0.300,0.200])*mA 

ibeam = {
'U25': 0.035*mA,
'U26': 0.051*mA,
'U27': 0.068*mA,
'U28': 0.088*mA,
'U29': 0.115*mA,
'U30': 0.150*mA,
'U31': 0.175*mA,
'U32': 0.192*mA,
'U33': 0.210*mA,  
'U34': 0.205*mA,
'U35': 0.178*mA,
'U36': 0.142*mA,
'U37': 0.110*mA,
'U38': 0.072*mA,
'U39': 0.043*mA,
'U40': 0.031*mA,
'O1': 0.300*mA,
'O2': 0.300*mA,
'O3': 0.300*mA,
'O4': 0.200*mA
         }

# --- Species kinetic energy at injection
#      * Set in a species named dictionary ekin  
#      * Set as Q*SourceBias where Q is the species charge state  

SourceBias = 35.*keV  # source voltage: set for Q_ref*SourceBias/A =>  4.9264706 keV/u for U 

ekin_birth = {key: SourceBias*charge_states[key] for key in sp.keys()} 


ekin_per_u = 12.*keV                             # target kinetic energy/u for LEBT post ES gap  

ekin_launch = {key: ekin_per_u*A_ref/Q_ref*charge_states[key] for key in sp.keys()} 



### Phase space volume setup:
# 1) Define the phase space volume contribution at launch before CSS
# 2) Set up initial P_theta and thermal emit as if the beam is launched at source
# 3) Redefine P_theta or thermal emit in accordance with launch mode
# 4) Inject rotation using frib-css-xy-load.py if P_theta does not vanish

# suffix "birth" denotes value at source extraction
# suffix "launch" or no suffice denotes value at launch in front of CSS

# Phase space volume contribution:

phase_space_vol_mode = 0    # 0: initial; 1: fully thermal; 2: fully P-theta 


# --- Species initial canonical angular momentum.
# 
#     Note: <P_theta>/(m*c) for each species behaves like a normalized emittance contributing to beam 
#     phase-space area. 
#
# Each ion species is assumed to be born from neutral to its full charge state instantly at birth. 
# birth_mode specifications below allow effective tuning of <P_theta> to cover a range of possibilities.     
#
# birth_mode = 0: Ions born at launch point B of applied field.
#                   * Option does nothing to adjust load to reflect initial distribution correlations 
#                     from being born at different magnetic field values within the ECR ion source.  
#              1: Ions born with statistical canonical angular momentum specified by species 
#                 with values set by a species dictionary: 
#                        ptheta = {key: value} value = <P_theta>/(m*c) [m-rad] 
#                   * Allows detailed tuning of distribution.
#              2: Ions in each species are born consistent with a specified field value and rms beam size.  
#                   * Allows adjusting <P_theta> consistent with beam being born a specified size 
#                     in a specified magnetic field.  
#                        rms_birth = sqrt(<r^2>) at birth [m]
#                        bz0_birth = B_z(r=0) at birth [Tesla]                  
#              3: Ions are born randomly between two peaks in the ECR ion source B-field. 
#                   * Each ion is born at a different B-field using a uniform spatial
#                     distribution between peaks of ECR field. 
#
#  * Angular velocities are added to the ions to conserve the linear approx 
#    statistical canonical angular momentum in birth_mode = 1,2 and 3  
#  * birth_mode = 2 only take into account the linear approx solenoidal (B_z) field component of 
#    source field. 
#  * birth_mode = 2 can be set with the average field value between the two peak values of the ECR 
#    ion source field and the likely size of the beam in the ECR to give an approx guess 
#    consistent with a particular ECR. 
#  * For birth_mode = 2, rms_launch = 2.8 mm, bz0_launch = 2.15 T
#    gives <P_theta>/(m*c) = 0.384 mm-mrad for U33 with the present setup.  


birth_mode = 2   # reminder: birth_mode 1 & 2 requries user input below


## -- Input for birth_mode == 1

ptheta_birth = {  # Values below are rounded numbers obtained from default values of birth_mode = 2 
'U25': 0.23*mm*mr,
'U26': 0.23*mm*mr,
'U27': 0.24*mm*mr,
'U28': 0.25*mm*mr,
'U29': 0.26*mm*mr,
'U30': 0.27*mm*mr,
'U31': 0.28*mm*mr,
'U32': 0.29*mm*mr,
'U33': 0.30*mm*mr,  
'U34': 0.31*mm*mr,
'U35': 0.32*mm*mr,
'U36': 0.32*mm*mr,
'U37': 0.33*mm*mr,
'U38': 0.34*mm*mr,
'U39': 0.35*mm*mr,
'U40': 0.36*mm*mr,
'O1': 0.13*mm*mr,
'O2': 0.27*mm*mr,
'O3': 0.40*mm*mr,
'O4': 0.54*mm*mr
         }  


## -- Input for birth_mode == 2

rms_birth_i = 2.82*mm  # Same birth radius as source aperture size
                       #  2.82 mm = sqrt(2.)*(r_extractor/2.)
                       # (factor 2. due to rms x equiv, sqrt(2.) due to rms r not rms x)

rms_birth = {  
'U25': rms_birth_i,
'U26': rms_birth_i,
'U27': rms_birth_i,
'U28': rms_birth_i,
'U29': rms_birth_i,
'U30': rms_birth_i,
'U31': rms_birth_i,
'U32': rms_birth_i,
'U33': rms_birth_i,
'U34': rms_birth_i,
'U35': rms_birth_i,
'U36': rms_birth_i,
'U37': rms_birth_i,
'U38': rms_birth_i,
'U39': rms_birth_i,
'U40': rms_birth_i,
'O1': rms_birth_i,
'O2': rms_birth_i,
'O3': rms_birth_i,
'O4': rms_birth_i
         }
         
bz0_birth_i = 1.67   # [Tesla], default value 1.67 T is the average field between Venus ECR peaks 

bz0_birth = {  
'U25': bz0_birth_i,
'U26': bz0_birth_i,
'U27': bz0_birth_i,
'U28': bz0_birth_i,
'U29': bz0_birth_i,
'U30': bz0_birth_i,
'U31': bz0_birth_i,
'U32': bz0_birth_i,
'U33': bz0_birth_i,
'U34': bz0_birth_i,
'U35': bz0_birth_i,
'U36': bz0_birth_i,
'U37': bz0_birth_i,
'U38': bz0_birth_i,
'U39': bz0_birth_i,
'U40': bz0_birth_i,
'O1': bz0_birth_i,
'O2': bz0_birth_i,
'O3': bz0_birth_i,
'O4': bz0_birth_i
         }


if birth_mode == 2:
	ptheta_birth={}
	for ii in sp.keys():
		ptheta_birth[ii] = 0.5*sp[ii].charge*bz0_birth[ii]*rms_birth[ii]**2/(sp[ii].mass*clight)






# --- Species initial thermal phase-space area measures  
#      * Set in species-named dictionaries:
#          temp_p = transverse   temperature [eV]  (after extractor accel value)
#          temp_z = longitudinal temperature [eV]  (after extractor accel value)
#          emitnx = normalized x-plane thermal rms emittance ** NOT edge measure **
#          emitny = normalized y-plane thermal rns emittance ** NOT edge measure **
#      * Can set individual dictionary values as in ibeam to tune as necessary.   
# 
#  Approach: Allow setting of transverse thermal phase-space area by either a 
#            temperature specification 
#            OR a normalized rms emittance measure. Convert input values to a normalized rms 
#            edge emittance to work with the Warp loader.  
#
#  init_ps_spec = "emitn"    => set by emitnx and emitny 
#                 or "temp"  => set by init_temp 
# 
#  Notes:
#  * Guilliaume:  Ions likely 2-3 eV and electrons few to 100s of keV.  
#                 Ions not equilibrated with electrons.
#  * Previous simulations used 0.4*mm*mr thermal phase-space area 
#    launched outside of magnetic field.      
#  * We set only the *thermal* component of the emittance and the distribution without 
#    correlation terms. Later, correlation terms are set to produce the value of canonical angular 
#    momentum desired based on axial magnetic field values at launch or birth. 
#  * Syntax is used to allow different phase space area measures for each species. Initial runs 
#    may employ a common value for all species.    

init_ps_spec = "temp" 

temp_perp_ion  = 3.          # initial transverse ion temp [eV] 
temp_z_ion     = 3.          # initial z          ion temp [eV] 
emitn          = 0.4*mm*mr   # initial rms normalized emittance [m-rad]  ** NOT edge measure **

temp_p_birth = {key: temp_perp_ion for key in sp.keys()} 
temp_z_birth = {key: temp_z_ion    for key in sp.keys()} 

emitnx_birth = {key: emitn for key in sp.keys()}
emitny_birth = {key: emitn for key in sp.keys()}






# --- Species initial beam centroid 
# 
#   xc  = <x>   Initial centroid coordinate [m]
#   xcp = <x'>  Initial centroid angle [rad] 
# 
#   + Analogous y-plane measures 
#

xc  = {key: 0. for key in sp.keys()}
yc  = {key: 0. for key in sp.keys()}

xcp = {key: 0. for key in sp.keys()}
ycp = {key: 0. for key in sp.keys()}




r_extractor = 4.*mm   # ECR extraction aperture radius [m]
                      #   Value both for Venus and Artemis A,B ECR Sources 

rx_birth  = {key: r_extractor for key in sp.keys()}
ry_birth  = {key: r_extractor for key in sp.keys()}

rxp_birth = {key: 0. for key in sp.keys()}
ryp_birth = {key: 0. for key in sp.keys()}




emit_thermal_birth = {}

for s in sp.keys():
  ekin_i  = ekin_birth[s]
  betab_i = sqrt(2.*jperev*ekin_i/sp[s].sm)/clight   # NR beta associated the KE 
  #
  rb_i    = sqrt(rx_birth[s]*ry_birth[s])          # take mean measure in case inhomogeneous
  #
  emitn_i = sqrt(emitnx_birth[s]*emitny_birth[s])  # take mean measure in case inhomogeneous  
  temp_i  = temp_p_birth[s]  
  if init_ps_spec == "emitn":
    emit_thermal_birth[s] = 4.*emitn_i 
  elif init_ps_spec == "temp": 
    emit_thermal_birth[s] = betab_i*((2.*rb_i)/sqrt(2.))*sqrt(temp_i/ekin_i)
  else:
    raise Exception("Error: init_emit_spec not set properly") 


phase_space_vol_birth = {key: sqrt(ptheta_birth[key]**2 + emit_thermal_birth[key]**2) for key in sp.keys()}

ptheta_launch = {}
emit_thermal_launch = {}

for ii in sp.keys():
	if phase_space_vol_mode == 0:
		ptheta_launch[ii] = ptheta_birth[ii]
		emit_thermal_launch[ii] = emit_thermal_birth[ii]
	if phase_space_vol_mode == 1:
		ptheta_launch[ii] = 0.
		emit_thermal_launch[ii] = phase_space_vol_birth[ii]
	if phase_space_vol_mode == 2:
		ptheta_launch[ii] = phase_space_vol_birth[ii]
		emit_thermal_launch[ii] = 0.










 #--- Species initial beam size: elliptical cross-section beam 
      #* Set as species-named dictionaries with:
   
   #rx  = 2*sqrt(<x^2>)           Initial x-plane beam rms edge radius [m] 
   #rxp = 2*<x*x'>/sqrt(<x^2>)    Initial rms envelope angle   [rad] 
 
   #ry, ryp  are analogous y-plane measures 

 #Notes:
 #* Now set rx etc by fraction of ECR aperture size.  
 #* Previous simulations used betatron functions with a specific emittance value to get 
   #desired beam size. Code used was (y-plane similar):   
alpha_x_launch = -0.3796434382
beta_x_launch  = 3.971339547
gamma_x_launch = (1. + alpha_x_launch**2)/beta_x_launch

#emitn_edge = 0.4*mm*mr   # norm rms edge emittance used to set beam size 

rx = {}
ry = {}
rxp = {}
ryp = {}

for ii in sp.keys():
	gamma = 1. + ekin_launch[ii]*jperev/(sp[ii].mass*clight**2)
	v_z = clight*sqrt(1. - 1./gamma**2)
	emit_edge = phase_space_vol_birth[ii]/(gamma*v_z/clight)
	rx[ii] = sqrt(emit_edge*beta_x_launch)
	ry[ii] = sqrt(emit_edge*beta_x_launch)
	rxp[ii] = -sqrt(emit_edge/beta_x_launch)*alpha_x_launch
	ryp[ii] = -sqrt(emit_edge/beta_x_launch)*alpha_x_launch

#ref_gamma_post_gap = 1. + ekin_per_u*jperev/(amu*clight**2)
#ref_vel_post_gap = clight*sqrt(1. - 1./ref_gamma_post_gap**2)
#ref_brho_post_gap = ref_gamma_post_gap*ref_vel_post_gap*A_ref*amu/(Q_ref*jperev)

#v_z_ref   = sqrt(2.*jperev*Q_ref*SourceBias/m_ref)  # nonrel approx ref z-velocity 
#gamma_ref = 1./sqrt(1.-(v_z_ref/clight)**2)         # ref axial gamma factor (approx 1) 
#emit_edge = emitn_edge/(gamma_ref*v_z_ref/clight)   # unnormalized rms edge emittance 

#rx  = sqrt(emit_edge*beta_x)             # envelope x-edge 
#rxp = -sqrt(emit_edge/beta_x)*alpha_x    # envelope x-angle 









#
# Define transverse simulation grid and properties of self-field solver 
#

# --- Single particle run
#       single_particle = flag to control whether to advance as a single particle run 
#                         True  => advance as single particles in applied field
#                         False => Advance with self-field and any neutralization factors 
single_particle = False   

# --- Symmetries.  Set for increased statistical efficiency.  These should
#     only be used in cases where lattice symmetries and initial beam
#     conditions and loads allow.  For no symmetry, set both options to false.
w3d.l2symtry = false     # 2-fold symmetry
w3d.l4symtry = false     # 4-fold symmetry

# -- Grid extent and increments
#     Chosen without consideration of symmetry. Consistently reset later 
#     if symmetry operations are set.  

r_grid = 8.*cm    # radial extent of grid 
n_grid = 400      # number grid cells (no symmetries) 


#
# Particle moving
#

z_launch  = 0.  # Axial position in lattice to launch  beam 
                       #   = ecr_z_extr defined in lattice file.
#z_adv     = 69.2       # Axial position in lattice to advance to (69.2 was end point before bend was introduced)  
#z_adv     = 69.581900 + 0.75       # Center of D5 dipole + half of axial length of 3d dipole field data
z_adv     = 6.       # End of CSS


ds = 2.*mm             # Axial advance step 

# 
# Particle loading 
#

nppg = 100.    # number of particles per grid cell


# Distribution loads type 
#  Comment: 
#   * See stptcl routine for how loading works. 
#   * Neutralization is NOT taken into account in the loading.    
#   * At present, all species are loaded with the same value of distrbtn.  
#   * rms equivalent beam loaded with the specified distribution from:  
#     KV => KV Distribution 
#
#     SG => semi-Gaussian distribution 
#             (KV density and local Gaussian angle spread about KV flutter) 
#
#     TE => Pseudoequilibrium with Thermal  Equilibrium form 
#     WB => Pseudoequilibrium with Waterbag Equilibrium form
#             The Pseudoequilibrium distributions use continuous focused 
#             equilibrium forms which are canoically transformed to AG 
#             symmetry of the lattice. 
#
#     For more info on loads, see review paper:
#       Lund, Kikuchi, and Davidson, PRSTAB 12, 114801 (2009) 


#w3d.distrbtn = "KV"          # initial KV distribution
#w3d.distrbtn = "TE"          # initial thermal distribution
w3d.distrbtn = "WB"           # initial waterbag distribution
#w3d.distrbtn = "SG"          # initial semi-Gaussian distribution 
