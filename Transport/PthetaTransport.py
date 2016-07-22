"""
Python input script for a Warp xy slice simulation of 
the FRIB front end.  For help and/or additional information contact:

     Steve Lund     lund@frib.msu.edu    (517) 908-7291

For documentation see the Warp web-site:

     https://warp.lbl.gov

All code inputs are mks with the exception of particle kinetic energy (eV).   

To run this Warp script:
  Interactive (return to interpreter after executing):

    % python -i frib-front-xy.py 

  Non-interactive (exit interpreter after executing):

    % python frib-front-xy.py
"""



## from warp.particles.singleparticle  import TraceParticle
## TraceParticle(x,y,z,..., enforceinitboundaries = False)








# Load Warp and various script packages 
from warp        import *               # Warp code 
from warp.utils.errorcheck  import checksymmetry   # Check for input errors
from warp.particles.singleparticle  import TraceParticle
#from runcounters import *               # Counter for parametric runs 

# Set informational labels included on all output cgm plots.   
top.pline2   = "xy Slice Simulation: FRIB Front End" 
top.pline1   = " "   # Add more info, if desired.  

# Invoke setup routine for graphics and output files (THIS IS MANDATORY)
setup()

# Set runmaker - included in informational labels on output plots
top.runmaker = "SMLund"

# Beam parameters for simulation
#
#   Other than numerical parameters, the physical simulation is specified
#   by the numbers input immediately below.  

# --- Define species
# 
#     Syntax of later use:
#
#       for us in U_species:
#         us.ppzx()
#
#     Dictionaries also setup later. 
















# Solenoid Field Strength (in Telsa)

SolenoidStrength = 0.45

# Lattice Length

LatticeLength = 0.8 # (times 5)

# Number of Periods

NumberofPeriods = 1

# Initial Radius

InitialRadius = 0.002

# set linear or nonlinear fields

abc = "lin"

# set neutralization factor

NeutralizationFactor = 0.75

#set velocity correction method: 0 - no correction, 1 - dBdz only, 2 - dBdz + d2Edz2

CorrectionMode = 1

# Number of particles per grid (default = 400, adjust for quick checks)

ParticleNumber = 400

# Envelope Plot Limits:
plotxmin = "e"
plotxmax = "e"
plotymin = 0.
plotymax = "e"


CorrectionFactor = 0.5

s4p2_zc  = 0.3



















# Setup handling for species weights 
# Via Dave Grote: need to setup scaling by hand for multi-species 
#  top.wpid = nextpid()     # setup/allocate pid array on generate 
#
#  species.getw()  = s.w   gets weights of species  [equivalent to species.w] (top.wpid must be set to work)  
#  species.getweights() returns product of species.sw*species.w 
#
#  put in adjustment using beforeloadrho() 
# 
#  beam.w * beam.sw = total weight = beam.getweights() 
# 

U_charge_states = [33,34,25,26,27,28,29,30,31,32,35,36,37,38,39,40]
O_charge_states = [1,2,3,4]

U_species = [Species(type=Uranium,charge_state=i,name='U%d'%i) for i in U_charge_states]
O_species = [Species(type=Oxygen, charge_state=i,name='O%d'%i) for i in O_charge_states]

U_ns = len(U_species) 
O_ns = len(O_species) 

# --- --- make abbreviated dictionary of species for later diagnostics 
#         sp.keys()
sp_U = {U_species[i].name:U_species[i] for i in range(U_ns)}
sp_O = {O_species[i].name:O_species[i] for i in range(O_ns)}

sp = {}
sp.update(sp_U)
sp.update(sp_O)

# --- Setup for variable particle weights on generate.
#       Note: pid array elements hold particle properties and are consistently mirrored when particles scraped  
#             nextpid() gets next pid array index 

top.wpid = nextpid()       # set pid index for variable weights: will initialize on generate when set  
uzp0pid  = nextpid() - 1   # set pid index to hold particle initial uz to scale weights (set after generate)  


# --- --- assign species colors using values which may help distinguish on plots: target magenta/green
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

# --- Target species tuple for later use 
sp_target = ['U33','U34']

# --- Reference species as average of target 
A_ref = 0. 
Q_ref = 0. 
for ii in sp_target:
  s = sp[ii]
  A_ref += s.mass/amu 
  Q_ref += s.charge/echarge 
A_ref = A_ref/len(sp_target) 
Q_ref = Q_ref/len(sp_target)  

m_ref = A_ref*amu 

# --- Set species properties for injection 
#        The beam kinetic energy (ekin) and axial velocity (vbeam) should not both
#        be set unless done so consistently.  If one is zero the code sets from
#        the other on generation:
#           vbeam = 0    => set vbeam from ekin
#           ekin  = 0    => set ekin  from vbeam
#

# --- --- unneutralized electric currents ... array elements corresponding to charge state arrays  
mA = 1.e-3
U_ibeam = array([0.210,0.205,0.035,0.051,0.068,0.088,0.115,0.150,0.175,0.192,0.178,0.142,0.110,0.072,0.043,0.031])*mA/10000. 
O_ibeam = array([0.300,0.300,0.300,0.200])*mA/10000. 

# --- --- kinetic energy 
SourceBias = 35.*keV  # source voltage: set for Q_ref*SourceBias/A =>  4.9264706 keV/u for U 

U_ekin = array(U_charge_states)*SourceBias
O_ekin = array(O_charge_states)*SourceBias 

# --- --- ion temp
#  Guilliaume:  ions likely 2-3 eV and electrons few - 100s of keV.  Ions not equilibrated with electrons.   

# --- --- beam size via edge reference emittance and Twiss parameters:

alpha_x = 0.
#beta_x  = 12.9*cm
beta_x  = 424.*cm
gamma_x = (1. + alpha_x**2)/beta_x 

alpha_y = 0.
#beta_y  = 12.9*cm
beta_y  = 424.*cm
gamma_y = (1. + alpha_y**2)/beta_y

mr = 0.001
emitn_edge = 0.4*mm*mr 

v_z_ref   = sqrt(2.*jperev*Q_ref*SourceBias/m_ref)
gamma_ref = 1./sqrt(1.-(v_z_ref/clight)**2)  
emit_edge = emitn_edge/(gamma_ref*v_z_ref/clight) 

top.lrelativ   = 1                # turn on relativity 

#r_x  = sqrt(emit_edge*beta_x)
#r_y  = sqrt(emit_edge*beta_y) 
#rp_x = -sqrt(emit_edge/beta_x)*alpha_x 
#rp_y = -sqrt(emit_edge/beta_y)*alpha_y 

r_x  = InitialRadius*sqrt(2)
r_y  = InitialRadius*sqrt(2)
rp_x = 0.
rp_y = 0.

# --- thermal velocity and energy (eV) of ref particle from emittance 
vt = v_z_ref*emit_edge/(2.*r_x) 
Et = 0.5*m_ref*vt**2/jperev 

# --- intrinsic thermal emittance scale 
Et_therm = 3.   # Guilliaume's estimated ion temp scale (eV) 
vt_therm = sqrt(2.*(jperev*Et_therm)/m_ref)
emit_therm  = 2.*r_x*vt_therm/v_z_ref
emitn_therm = (gamma_ref*v_z_ref/clight)*emit_therm


# Ratio of thermal to edge emittance suggests value of P_theta contributing to effective emittance 
# emit_therm/emit_edge = 0.10  => most beam PS area from P_theta  

for i in range(U_ns):
  Usp = U_species[i]
  Usp.ekin   = U_ekin[i]           # kinetic energy of beam particle [eV]
  Usp.vbeam  = 0.                  # beam axial velocity [m/sec] (set from ekin if 0) 
  Usp.ibeam  = U_ibeam[i]          # beam current [amps] 
  Usp.emitnx = emitn_therm         # beam x-emittance, rms edge [m-rad] 
  Usp.emitny = emitn_therm         # beam y-emittance, rms edge [m-rad]
  #Usp.emitnx = emitn_edge          
  #Usp.emitny = emitn_edge          
  Usp.vthz   = 0.                  # axial velocity spread [m/sec] 

for i in range(O_ns):
  Osp = O_species[i]
  Osp.ekin   = O_ekin[i] 
  Osp.vbeam  = 0.
  Osp.ibeam  = O_ibeam[i]
  Osp.emitnx = emitn_therm
  Osp.emitny = emitn_therm
  #Osp.emitnx = emitn_edge 
  #Osp.emitny = emitn_edge 
  Osp.vthz   = 0.

# Calculate vbeam and other species quantities 
derivqty()

# ---  Calculate and printout Q/M by species and store in a dictionary 
sp_qovm = {}
print("Species Charge to Mass Ratios:")
for ii in sort(sp.keys()):
  s = sp[ii]
  qovm = (s.charge/echarge)/(s.mass/amu)
  sp_qovm.update({ii:qovm})
  print("   Species: "+ii+" Q/A = %s"%qovm)

# ---  Calculate and printout rigidity by species and store in a dictionary 
sp_brho = {}
print("Species Rigidity:")
for ii in sort(sp.keys()):
  s = sp[ii]
  gamma = 1./sqrt(1.-(s.vbeam/clight)**2)
  brho  = gamma*s.mass*s.vbeam/s.charge
  sp_brho.update({ii:brho})
  print("   Species: "+ii+" [B rho] = %s T-m"%brho)

# --- Diagnostic plot of [B rho] vs Q/A for species 


def plt_diag_bro(label=None):
  if label == None: label = " "
  brho_min = largepos 
  brho_max = 0.  
  for ii in sp.keys():
    s = sp[ii]
    js = s.js 
    #
    weight = sum(s.sw*s.w) 
    #
    vbeam = sum( (s.sw*s.w)*s.getvz() )/weight
    gammabeam = 1./sqrt(1.-(vbeam/clight)**2)      
    brho  = s.mass*gammabeam*vbeam/s.charge
    #
    brho_min = min(brho,brho_min) 
    brho_max = max(brho,brho_max) 
    #
    plt(ii,sp_qovm[ii],brho,tosys=1,color=s.color) 
    #
  [qovm_min,qovm_max] = [minnd(sp_qovm.values()),maxnd(sp_qovm.values())]
  qovm_pad = 0.1*(qovm_max - qovm_min)
  brho_pad = 0.1*(brho_max - brho_min)
  #
  limits(qovm_min-qovm_pad,qovm_max+qovm_pad,brho_min-brho_pad,brho_max+brho_pad) 
  ptitles(label,"Q/A","[B rho] [Tesla-m]",)
  fma()


#
# Beam centroid and rms envelope initial conditions at s=0      
#    
#   x0:   initial x-centroid xc = <x> [m]
#   y0:   initial y-centroid yc = <y> [m]
#   xp0:  initial x-centroid angle xc' = <x'> = d<x>/ds [rad]
#   yp0:  initial y-centroid angle yc' = <y'> = d<y>/ds [rad]
#
#   a0:   initial x-envelope edge a = 2*sqrt(<(x-xc)^2>) [m]
#   b0:   initial y-envelope edge b = 2*sqrt(<(y-yc)^2>) [m]
#   ap0:  initial x-envelope angle ap = a' = d a/ds [rad]
#   bp0:  initial y-envelope angle bp = b' = d b/ds [rad]

for i in range(U_ns):
  Usp = U_species[i]
  # --- centroid 
  Usp.x0  = 0.
  Usp.y0  = 0.
  Usp.xp0 = 0.
  Usp.yp0 = 0.
  # --- envelope 
  Usp.a0   = r_x
  Usp.b0   = r_y
  Usp.ap0  = rp_x
  Usp.bp0  = rp_y

for i in range(O_ns):
  Osp = O_species[i]
  # --- centroid 
  Osp.x0  = 0.   
  Osp.y0  = 0.   
  Osp.xp0 = 0.   
  Osp.yp0 = 0.   
  # --- envelope 
  Osp.a0   = r_x              
  Osp.b0   = r_y           
  Osp.ap0  = rp_x  
  Osp.bp0  = rp_y  


#
# Setup Lattice  
#
ekin_per_u = 12.*keV                             # target kinetic energy/u for LEBT 
StandBias = A_ref*ekin_per_u/Q_ref - SourceBias  # Bias of Injector Column  

Bias = StandBias + SourceBias


 #--- Venus ECR Source 
     #Comment: Must have same z-grids for linear and nonlinear forms.  Minimal error checking to enforce this. 

 #--- --- element specification 

ecr_shift  = 11.*cm                 # shift of ecr from lattice file spec to make room for s4p1 
ecr_z_extr = 0.01 - ecr_shift  # z-location of beam extraction aperture in simulation coordinates     
ecr_sc     = 1.0                    # scale factor to muliply field data by 
ecr_typ    = "lin"                  # type: "lin" = linear optics fields or "nl" = nonlinear r-z field  

# --- --- linear element data  
fi = PRpickle.PR("lat_ecr_venus.lin.20140602.pkl")
ecr_bz_extr = fi.ecr_venus_bz_extr
ecr_dz = fi.ecr_venus_dz 
ecr_nz = fi.ecr_venus_nz  
ecr_z_m     = fi.ecr_venus_z_m
ecr_zm_extr = fi.ecr_venus_z_extr  # extraction location on z_m mesh field    
ecr_bz0_m   = fi.ecr_venus_bz0_m
ecr_bz0p_m  = fi.ecr_venus_bz0p_m
fi.close() 

ecr_zlen  = ecr_z_m.max() - ecr_z_m.min()                 # length ecr field mesh  
ecr_zmmin = ecr_z_extr - (ecr_zm_extr - ecr_z_m.min())    # start point of ecr field mesh in sim coordinates 
ecr_zmmax = ecr_z_extr + (ecr_z_m.max() - ecr_zm_extr)    # end   point of ecr field mesh in sim coordinates

ecr_lin_id = addnewmmltdataset(zlen=ecr_zlen,ms=ecr_bz0_m,msp=ecr_bz0p_m,nn=0,vv=0)

# --- --- define venus ecr fields  
if ecr_typ == "lin":
  ecr = addnewmmlt(zs=ecr_zmmin,ze=ecr_zmmax,id=ecr_lin_id,sc=ecr_sc) 
elif ecr_typ == "nl":
  #addnewbgrd(xs=0.,zs=s41_zc-s4_zlen/2.,ze=s41_zc+s4_zlen/2.,id=s4_nl_id,func=s41_scale)
  raise Exception("No ECR Venus Nonlinear Applied Fields Defined") 
  ecr = None
else:
  print("Warning: No ECR Applied Fields Defined") 
  ecr = None


# --- S4 solenoids 
#     Comment: linear and nonlinear variants must have same z-grid.  Minimal error checking only for input 
#              consistency.   

# --- --- element specification 

#s4p1_zc  = 66.956900   # S4 1: z-center  
#s4p1_str = 0.6 # 0.754 # S4 1: peak on-axis B_z field strength [Tesla]
#s4p1_typ = "nl"        # S4 1: type: "lin" = linear optics fields or "nl" = nonlinear r-z field  

#s4p2_zc  = LatticeLength / 2   # S4 2: z-center 
#s4p2_str = 2.0 # 0.617 # s4 2: peak on-axis B_z field strength [Tesla]
s4p2_str = SolenoidStrength # 0.617 # s4 2: peak on-axis B_z field strength [Tesla]
s4p2_typ = abc        # S4 1: type: "lin" = linear optics fields or "nl" = nonlinear r-z field  

# --- --- linear element data  
#fi = PRpickle.PR("lat_s4.lin.20140603.pkl")
fi = PRpickle.PR("lat_s4.lin.20141031.pkl")
s4_dz  = fi.s4_dz 
s4_nz  = fi.s4_nz  
s4_z_m = fi.s4_z_m 
s4_bz0_m   = fi.s4_bz0_m
s4_bz0p_m  = fi.s4_bz0p_m
fi.close() 

s4_zlen = s4_z_m.max() - s4_z_m.min() 
s4_lin_id = addnewmmltdataset(zlen=s4_zlen,ms=s4_bz0_m,msp=s4_bz0p_m,nn=0,vv=0)

# --- --- nonlinear element field data 
#fi = PRpickle.PR('lat_s4.rz.20140603.pkl') 
fi = PRpickle.PR('lat_s4.rz.20141031.pkl') 
#
s4_len_coil   = fi.s4_len_coil 
s4_len_magnet = fi.s4_len_magnet 
s4_r_coil_i   = fi.s4_r_coil_i 
s4_r_coil_o   = fi.s4_r_coil_o
#
if fi.s4_nz != s4_nz: raise Exception("S4: Nonlinear field model nz not equal to linear field model nz") 
s4_dr   = fi.s4_dr
s4_nr   = fi.s4_nr 
s4_r_m  = fi.s4_r_m 
s4_br_m_in = fi.s4_br_m
s4_bz_m_in = fi.s4_bz_m
fi.close() 

# --- --- nonlinear element vector potential data 
#fi = PRpickle.PR('lat_s4.at.20140603.pkl') 
fi = PRpickle.PR('lat_s4.at.20141031.pkl') 
#
if fi.s4_nz != s4_nz: raise Exception("S4: Nonlin Vector potential model nz not equal to nonlinear/linear model nz")
if fi.s4_nr != s4_nr: raise Exception("S4: Nonlin Vector potential model nr not equal to nonlinear model nr")
s4_at_m  = fi.s4_at_m
fi.close() 

# Axisymmetric b-field arrays must be 3d shape (nr+1,arb,nz+1) to load into Warp  
s4_br_m = fzeros((s4_nr+1,1,s4_nz+1))  
s4_br_m[:,0,:] = s4_br_m_in
s4_bz_m = fzeros((s4_nr+1,1,s4_nz+1))
s4_bz_m[:,0,:] = s4_bz_m_in

s4_nl_id = addnewbgrddataset(dx=s4_dr,dy=1.,zlength=s4_zlen,bx=s4_br_m,bz=s4_bz_m,rz = true)  # pass arb dy to avoid error trap  

s4_aspect = s4_r_coil_i/s4_len_coil 

## --- --- define solenoid s4 1 
#if s4p1_typ == "lin":
  #s4p1 = addnewmmlt(zs=s4p1_zc-s4_zlen/2.,ze=s4p1_zc+s4_zlen/2.,id=s4_lin_id,sc=s4p1_str) 
#elif s4p1_typ == "nl":
  #s4p1 = addnewbgrd(xs=0.,zs=s4p1_zc-s4_zlen/2.,ze=s4p1_zc+s4_zlen/2.,id=s4_nl_id,sc=s4p1_str)
#else:
  #print("Warning: No S4 1st Solenoid Applied Fields Defined") 
  #s4p1 = None

# --- --- define solenoid s4 2 
if s4p2_typ == "lin":
  s4p2 = addnewmmlt(zs=s4p2_zc-s4_zlen/2.,ze=s4p2_zc+s4_zlen/2.,id=s4_lin_id,sc=s4p2_str) 
elif s4p2_typ == "nl":
  s4p2 = addnewbgrd(xs=0.,zs=s4p2_zc-s4_zlen/2.,ze=s4p2_zc+s4_zlen/2.,id=s4_nl_id,sc=s4p2_str)
else:
  print("Warning: No S4 2nd Solenoid Applied Fields Defined") 
  s4p2 = None
  
# --- --- define solenoid s4 2 
if s4p2_typ == "lin":
  s4p2 = addnewmmlt(zs=s4p2_zc+ LatticeLength -s4_zlen/2.,ze=s4p2_zc+LatticeLength+s4_zlen/2.,id=s4_lin_id,sc=s4p2_str) 
elif s4p2_typ == "nl":
  s4p2 = addnewbgrd(xs=0.,zs=s4p2_zc-s4_zlen/2.,ze=s4p2_zc+s4_zlen/2.,id=s4_nl_id,sc=s4p2_str)
else:
  print("Warning: No S4 2nd Solenoid Applied Fields Defined") 
  s4p2 = None

# --- --- define solenoid s4 2 
if s4p2_typ == "lin":
  s4p2 = addnewmmlt(zs=s4p2_zc+ LatticeLength*2 -s4_zlen/2.,ze=s4p2_zc+LatticeLength*2+s4_zlen/2.,id=s4_lin_id,sc=s4p2_str) 
elif s4p2_typ == "nl":
  s4p2 = addnewbgrd(xs=0.,zs=s4p2_zc-s4_zlen/2.,ze=s4p2_zc+s4_zlen/2.,id=s4_nl_id,sc=s4p2_str)
else:
  print("Warning: No S4 2nd Solenoid Applied Fields Defined") 
  s4p2 = None

# --- --- define solenoid s4 2 
if s4p2_typ == "lin":
  s4p2 = addnewmmlt(zs=s4p2_zc+ LatticeLength*3 -s4_zlen/2.,ze=s4p2_zc+LatticeLength*3+s4_zlen/2.,id=s4_lin_id,sc=s4p2_str) 
elif s4p2_typ == "nl":
  s4p2 = addnewbgrd(xs=0.,zs=s4p2_zc-s4_zlen/2.,ze=s4p2_zc+s4_zlen/2.,id=s4_nl_id,sc=s4p2_str)
else:
  print("Warning: No S4 2nd Solenoid Applied Fields Defined") 
  s4p2 = None

# --- --- define solenoid s4 2 
if s4p2_typ == "lin":
  s4p2 = addnewmmlt(zs=s4p2_zc+ LatticeLength*4 -s4_zlen/2.,ze=s4p2_zc+LatticeLength*4+s4_zlen/2.,id=s4_lin_id,sc=s4p2_str) 
elif s4p2_typ == "nl":
  s4p2 = addnewbgrd(xs=0.,zs=s4p2_zc-s4_zlen/2.,ze=s4p2_zc+s4_zlen/2.,id=s4_nl_id,sc=s4p2_str)
else:
  print("Warning: No S4 2nd Solenoid Applied Fields Defined") 
  s4p2 = None

# Define vector potential function for both linear and nonlinear solenoid magnetic fields  
def getatheta(r):
  # --- gather vector potential 
  n = len(r) 
  at = zeros(n)
  at_scratch = zeros(n) 
  #z  = top.zbeam*ones(n)
  
  zzz = top.zbeam
  
  if zzz > LatticeLength:
	  zzz = zzz % LatticeLength
  
  z = zzz*ones(n)
  
  if   top.zbeam >= ecr_zmmin and top.zbeam <= ecr_zmmax:
    # --- contribution in venus 
    if ecr_typ == "lin":
      getgrid1d(n,z,at_scratch,ecr_nz,ecr_sc*ecr_bz0_m,ecr_zmmin,ecr_zmmax)
      at_scratch = at_scratch*r/2.
    elif ecr_typ == "nl":
       raise Exception("Vector Potential: ECR Nonlinear not defined")  
    else:
       raise Exception("Vector Potential: ECR not defined") 
    at += at_scratch
  #if top.zbeam >= s4p1_zc-s4_zlen/2. and top.zbeam <= s4p1_zc+s4_zlen/2.:
    ## --- contribution from 1st s4 
    #if s4p1_typ == "lin": 
      #getgrid1d(n,z,at_scratch,s4_nz,s4p1_str*s4_bz0_m,s4p1_zc-s4_zlen/2.,s4p1_zc+s4_zlen/2.)
      #at_scratch = at_scratch*r/2.
    #elif s4p1_typ == "nl":
      #getgrid2d(n,r,z,at_scratch,s4_nr,s4_nz,s4p1_str*s4_at_m,s4_r_m.min(),s4_r_m.max(), 
                #s4p1_zc-s4_zlen/2.,s4p1_zc+s4_zlen/2.)
    #else:
      #raise Exception("Vector Potential: S4.1 not defined")
    #at += at_scratch
  if zzz >= s4p2_zc-s4_zlen/2. and zzz <= s4p2_zc+s4_zlen/2.:
    # --- contribution from 2nd s4
    if s4p2_typ == "lin": 
      getgrid1d(n,z,at_scratch,s4_nz,s4p2_str*s4_bz0_m,s4p2_zc-s4_zlen/2.,s4p2_zc+s4_zlen/2.)
      at_scratch = at_scratch*r/2.
    elif s4p2_typ == "nl": 
      getgrid2d(n,r,z,at_scratch,s4_nr,s4_nz,s4p2_str*s4_at_m,s4_r_m.min(),s4_r_m.max(), 
                s4p2_zc-s4_zlen/2.,s4p2_zc+s4_zlen/2.)
    else:
      raise Exception("Vector Potential: S4.2 not defined")
    at += at_scratch
  return at 


## --- Grated Acceleration Gap
##   Note: for ideal zero-length gap:  top.lacclzl=true for zero length gap.  Accel given given by acclez*(accelze-acclzs) 
##   see dave grote email on caution on setting top.acclsw for gaps.   
##   Comment: Linear and nonlinear forms must have same axial grid.  Miminal error checking only for this.  

## --- --- element specification 
#gag_zc  = 67.811564  # Grated Accel Gap: z-center  
#gag_typ = "nl"       # Grated Accel Gap: type: "ideal" = Short gap kick, "lin" = linear r-z field imported, "nl" = nonlinear r-z field imported   

## --- --- linear element data  
## fi = PRpickle.PR("lat_gag.lin.20140624.pkl")  # Original Warp model with simplified geometry  
#fi = PRpickle.PR("lat_gag.lin.20141029.pkl")    # Poisson model with high detail 
#gag_dz = fi.gag_dz0 
#gag_nz = fi.gag_nz0  
#gag_z_m     = fi.gag_z0_m  
#gag_ez0_m   = fi.gag_ez0_m
#gag_ez0p_m  = fi.gag_ez0p_m
#fi.close() 

#gag_zlen = gag_z_m.max() - gag_z_m.min() 

#gag_lin_id = addnewemltdataset(zlen=gag_zlen,es=gag_ez0_m,esp=gag_ez0p_m,nn=0,vv=0)

## --- --- nonlinear element data
##fi = PRpickle.PR('lat_gag.rz.20140624.pkl')  # Original Warp model with simplified geometry  
#fi = PRpickle.PR('lat_gag.rz.20141029.pkl')   # Poisson model with high detail 
#if fi.gag_nz != gag_nz: raise Exception("GAG: Nonlinear and linear field model nz not equal") 
#gag_nr = fi.gag_nr
#gag_dr = fi.gag_dr
#gag_r_m = fi.gag_r_m
#gag_z_m_cen = fi.gag_z_m_cen
#gag_phi_m    = fi.gag_phi_m 
#gag_er_m_in  = fi.gag_er_m
#gag_ez_m_in  = fi.gag_ez_m
#fi.close() 

#gag_zlen = gag_z_m.max()-gag_z_m.min()          # axial length nonlin/lin structure on mesh 

#gag_zs   = gag_zc - (gag_z_m_cen-gag_z_m.min()) # z_start of nonlin/lin mesh structure  
#gag_ze   = gag_zc + (gag_z_m.max()-gag_z_m_cen) # z_end   of nonlin/lin mesh structure 

## Geometry parameters ?? Read these in grated gap file to be safe .. ** must be consistent with input data ** ?? 
#gag_rp = 7.3*cm                  # pipe radius of inner extent of rings in grated gap 
#gag_col_zs = gag_zc - 11.989*cm  # z-start (at end   of biased upstream pipe)     of grated gap mechanical structure 
#gag_col_ze = gag_zc + 15.611*cm  # z-end   (at start of grounded downstream pipe) of grated gap mechanical structure  

#gag_er_m = fzeros((gag_nr+1,1,gag_nz+1)) # Axisymmetric e-field arrays must be 3d shape (nr+1,arb,nz+1) to load into Warp  
#gag_er_m[:,0,:] = gag_er_m_in
#gag_ez_m = fzeros((gag_nr+1,1,gag_nz+1))
#gag_ez_m[:,0,:] = gag_ez_m_in

#gag_nl_id = addnewegrddataset(dx=gag_dr,dy=1.,zlength=gag_zlen,ex=gag_er_m,ez =gag_ez_m,rz = true) 

## --- --- define grated acceleration gap  
#if gag_typ == "ideal":
  #print("Wanning: No Ideal Acceleration Gap model yet implemented")  
  #gag = None 
#elif gag_typ == "lin":
  #gag = addnewemlt(zs=gag_zs,ze=gag_ze,id=gag_lin_id,sc=StandBias) 
#elif gag_typ == "nl":
  #gag = addnewegrd(xs=0.,zs=gag_zs,ze=gag_ze,id=gag_nl_id,sc=StandBias)
#else:
  #print("Warning: No Grated Acceleration Gap Applied Fields Defined") 
  #gag = None 


## --- D5 Bending Dipole 

## --- --- element specification 

#d5p1_zc  = 69.581900   # D5 1: z-center  
#d5p1_str = 1.0         # D5 1: Input field scale factor
#d5p1_typ = "3d"        # D5 1: type: "lin" = linear optics fields or "3d" = 3d field  

## --- --- nonlinear element data 
#fi = PRpickle.PR('lat_d5.3d.20140527.pkl') 
#d5_3d_nx = fi.d5_nx
#d5_3d_ny = fi.d5_ny
#d5_3d_nz = fi.d5_nz
#d5_3d_dx = fi.d5_dx
#d5_3d_dy = fi.d5_dy
#d5_3d_dz = fi.d5_dz
#d5_3d_x_m = fi.d5_x_m
#d5_3d_y_m = fi.d5_y_m
#d5_3d_z_m = fi.d5_z_m
#d5_3d_z_m_cen = fi.d5_z_m_cen
#d5_3d_bx_m = fi.d5_bx_m
#d5_3d_by_m = fi.d5_by_m
#d5_3d_bz_m = fi.d5_bz_m
#fi.close() 
#d5_3d_zlen = d5_3d_z_m.max() - d5_3d_z_m.min()

#d5_3d_id = addnewbgrddataset(dx=d5_3d_dx,dy=d5_3d_dy,zlength=d5_3d_zlen,bx=d5_3d_bx_m,by=d5_3d_by_m,bz =d5_3d_bz_m) 

## --- --- define dipole d5 
#if d5p1_typ == "lin":
  #print("Warning: No D5 1st Dipole Linear Applied Fields Defined")
  #d5p1 = None
#elif d5p1_typ == "nl":
  #d5p1 = addnewbgrd(dx=d5_3d_dx,dy=d5_3d_dy,xs=d5_3d_x_m.min(),ys=d5_3d_y_m.min(),
    #zs=d5p1_zc-d5_3d_zlen/2.,ze=d5p1_zc+d5_3d_zlen/2.,id=d5_3d_id,sc=d5p1_str)
#else:
  #print("Warning: No D5 1st Dipole Applied Fields Defined") 
  #d5p1 = None


# --- Neutralization specifications 
#     Break points z1 and z2 correspond to z values where neutralization is turned off and then back on so the beam is 
#     unneutralized in the grated acceleration gap.  

#neut_z1 = gag_zc - 20.90*cm    # z of neutralization stop post injector before grated gap: set via 1% of gap E_z field reached  
#neut_z2 = gag_zc + 22.28*cm  

neut_f1 = NeutralizationFactor                 # corresponding electron neutralization factors 
#neut_f2 = 0.01 
neut_ode = 1 - neut_f1

# --- Aperture specfications 
#     load scraper after generation of pic code  

r_p = 8.00*cm

#r_p_up   = 8.00*cm 
#r_p_down = 7.62*cm 

#r_ap   = [r_p_up,gag_rp,r_p_down]
#v_ap   = [SourceBias+StandBias,StandBias/2.,0.] 
#z_ap_l = [ecr_z_extr,  gag_col_zs, gag_col_ze]
#z_ap_u = [gag_col_zs,  gag_col_ze, d5p1_zc   ]

#r_p = max(r_ap)   # Max aperture in simulations 


#aperture = [] 
#for i in range(len(r_ap)):
  #rp = r_ap[i] 
  #v  = v_ap[i] 
  #zl = z_ap_l[i] 
  #zu = z_ap_u[i]
  ##
  #aperture.append( ZCylinderOut(radius=rp,zlower=zl,zupper=zu,condid="next") )


# --- Add a circular aperture particle scraper at pipe radius aperture. 
#       Could also use ParticleScraper(conductors=aperture) but 
#       setting prwall is faster for a simple cylinder. 
top.prwall = r_p    # reset this later consistent with actual aperture in range simulated in advances 

L_p = LatticeLength*5

aperture = ZCylinderOut(radius=r_p, length=largepos, zcent=L_p/2., condid="next")

# --- Lattice periodicity  
#top.zlatperi  = largepos # periodicity length [m]  
top.zlatperi  = L_p # periodicity length [m]  
top.zlatstrt  = 0.       # z of lattice start; added to element z's [m] 
                         #   (can use to change lattice phase) 

nadvance = NumberofPeriods
nstep = 500































# Define transverse simulation grid

# --- Symmetries.  Set for increased statistical efficiency.  These should
#     only be used in cases where lattice symmetries and initial beam
#     conditions and loads allow.  For no symmetry, set both options to false.
w3d.l2symtry = false     # 2-fold symmetry
w3d.l4symtry = false     # 4-fold symmetry

# -- Grid increments
#      First choose number grid cells without symmetry and reset 
#      consistent with symmetry options

n_grid = ParticleNumber # 200 previous # number grid cells (no symmetries) 

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
l_grid = 2.*r_p               # length edge of simulation grid [m]      
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

nppg = 100.    # number of particles per grid cell
top.npmax = int(nppg*pi*(r_x*r_y)/dx**2*sym_x*sym_y) # max initial particles loaded (each species) 

# Distribution loads type 
#  Comment: 
#     See stptcl routine for how loading works.   
#     At present, all species are loaded with the same value of distrbtn.   zbeam can be 
#     set before the generate for the beam location. 

z_launch  = 0. 
top.zbeam = z_launch

# rms equivalent beam loaded with the specified distribution form 
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
w3d.distrbtn = "WB"          # initial waterbag distribution
#w3d.distrbtn = "SG"          # initial semi-Gaussian distribution 

# --- random number options to use in loading 
w3d.xrandom  = "pseudo" # "digitrev"    # load x,y,z  with digitreverse random numbers 
w3d.vtrandom = "pseudo" # "digitrev"    # load vx, vy with digitreverse random numbers
w3d.vzrandom = "pseudo" # "digitrev"    # load vz     with digitreverse random numbers 
w3d.cylinder = true          # load a cylinder

# Particle moving
#

z_adv = nadvance*top.zlatperi   # Axial position in lattice to advance to 

top.lrelativ   =  false    # turn off relativistic kinematics
top.relativity = 0         # turn off relativistic self-field correction
                           #   to account for approx diamagnetic B-field of beam

wxy.ds = L_p/float(nstep)             # ds for part adv [m] 
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
#top.depos = "none"


# Turn on x-window plots, if desired; use winkill() to close interactively.  
#winon()


# Potential profile diagnostic primarily for initial beam 
def diag_plt_phi_ax(xmax=None,label=None):
  if xmax == None: xmax = max(w3d.xmesh.max(),w3d.ymesh.max())
  ixmax = sum(where(w3d.xmesh < xmax, 1, 0))
  iymax = sum(where(w3d.ymesh < xmax, 1, 0)) 
  if label == None: label = "Beam Potential at y,x = 0 b,r"
  #
  ix_cen = sum(where(w3d.xmesh < 0., 1, 0))
  iy_cen = sum(where(w3d.ymesh < 0., 1, 0))
  phix = getphi(iy=iy_cen)
  phiy = getphi(ix=ix_cen)
  phimin = min(phix[ixmax],phiy[iymax]) 
  #
  plg(phix,w3d.xmesh/mm)
  plg(phiy,w3d.ymesh/mm,color="red") 
  ptitles(label,"x,y [mm]","phi [V]", )
  limits(-xmax/mm,xmax/mm,phimin,'e') 


################################
# Particle simulation
################################

# Generate the xy PIC code.  In the generate, particles are allocated and
# loaded consistent with initial conditions and load parameters
# set previously.  Particles are advanced with the step() command later
# after various diagnostics are setup.
package("wxy"); generate()

# Read in diagnostics for applied lattice fields 
#execfile("diag_lattice.py") 

# Install conducting aperture on mesh 
installconductors(aperture,dfill=largepos)

# Check that inputs are consistent with symmetries (errorcheck package function)
checksymmetry()

# Carry out an initial unneutralized field solve with conducting pipe after generate 

for s in sp.values():       
  s.w   = 1.      # Need full charge: set relative weight to unity 

loadrho() 
fieldsolve() 

# Make plot of initial unneutralized beam potential profile 
         
diag_plt_phi_ax(label="Initial Unneutralized Beam Potential at y,x = 0 b,r") 
fma()

diag_plt_phi_ax(label="Initial Unneutralized Beam Potential at y,x = 0 b,r",xmax=1.5*r_x)
fma()


# Setup variable weight species needs for neutralization and acceleration 

# --- set initial weight factors consistent with neutralization factor 
#       w0 = initial weight (same for all particles in species) 
#       species.w = array of variable weight factors 
for s in sp.values():       
  s.w0  = 1.-neut_f1
  #s.w   = 1.-neut_f1
  s.sw0    = s.sw       # save initial sw    (in case later changed)  
  s.vbeam0 = s.vbeam    # save initial vbeam (in case later changed)  

# --- save initial uzp for all species at once 
top.pgroup.pid[:,uzp0pid] = top.pgroup.uzp

# --- adjust weights  
@callfrombeforeloadrho
def adjustweights():
  for s in sp.values():
    s.w[:] = s.w0*s.pid[:,uzp0pid]/s.uzp

# Carry out explicit fieldsolve with adjusted rho consistent with neutralization 
loadrho() 
fieldsolve()


# Make plot of initial neutralized beam potential profile 
         
diag_plt_phi_ax(label="Initial f = %s Neutralized Beam Potential at y,x = 0 b,r"%(neut_f1))
fma()

diag_plt_phi_ax(label="Initial f = %s Neutralized Beam Potential at y,x = 0 b,r"%(neut_f1),xmax=1.5*r_x)
fma()


# Modify intital distribution loaded on generate to include canonical angular momentum
# 
bz0_extr   = getappliedfields(x=0.,y=0.,z=ecr_z_extr)[5]    # B_z on-axis at ECR extraction plane
bz0_launch = getappliedfields(x=0.,y=0.,z=z_launch)[5]      # B_z on-axis at simulation launch location 

inj_ang_mom = true

sp_krot_launch = {}
sp_krot_v      = {} 
for ii in sp.keys():
  s = sp[ii]
  # --- rigidity 
  gamma = 1./sqrt(1.-(s.vbeam/clight)**2)
  brho  = gamma*s.mass*s.vbeam/s.charge
  # --- rms calculation
  rms_launch = sqrt(average( (s.xp)**2 + (s.yp)**2 ))
  rms_extr   = sqrt( (s.a0/2.)**2 + (s.b0/2.)**2 )     # *** Replace with beam size at extraction ! ****
  # --- rot wavenumbers at launch and in vacuum v 
  krot_launch = (bz0_extr*rms_extr**2/rms_launch**2 - bz0_launch)/(2.*brho)
  krot_v      = bz0_extr/(2.*brho)
  # 
  sp_krot_launch.update({ii:krot_launch})
  sp_krot_v.update({ii:krot_v}) 
  #
  if inj_ang_mom: 
    s.uxp -= krot_launch*s.yp*s.uzp
    s.uyp += krot_launch*s.xp*s.uzp

#raise Exception("to here")














electrons_tracked = Species(type=Uranium, charge_state = 33)
ntrack = 2
electron_tracker = TraceParticle(js=electrons_tracked.jslist[0],
                     x=[1*mm, 0.],
                     y=[1*mm, 0.],
                     z=0.,
                     vx=[0., 1*10**-6*967654.00],
                     vy=[1*10**-6*967654.00, 0.],
                     savedata=1,
                     enforceinitboundaries = False)



















# --- make plots of initial rotation by species at launch point and if transported in vacuuo    

def diag_plt_krot_launch(): 
  for ii in sp.keys():
    plt(ii,sp_qovm[ii],sp_krot_launch[ii],tosys=1,color=sp[ii].color) 

  [qovm_min,qovm_max] = [minnd(sp_qovm.values()),maxnd(sp_qovm.values())]
  [krot_launch_min,krot_launch_max] = [minnd(sp_krot_launch.values()),maxnd(sp_krot_launch.values())]
  qovm_pad = 0.1*(qovm_max - qovm_min)
  krot_launch_pad = 0.1*(krot_launch_max - krot_launch_min)

  limits(qovm_min-qovm_pad,qovm_max+qovm_pad,krot_launch_min-krot_launch_pad,krot_launch_max+krot_launch_pad) 
  ptitles("Angular Phase Advance Wavenumber: Beam Launch","Q/A","Wavenumber [Rad/m]",)
  fma() 

diag_plt_krot_launch() 

def diag_plt_krot_v():
  for ii in sp.keys():
    plt(ii,sp_qovm[ii],sp_krot_v[ii],tosys=1,color=sp[ii].color)

  [qovm_min,qovm_max] = [minnd(sp_qovm.values()),maxnd(sp_qovm.values())]
  [krot_v_min,krot_v_max] = [minnd(sp_krot_v.values()),maxnd(sp_krot_v.values())]

  krot_v_pad = 0.1*(krot_v_max - krot_v_min)
  qovm_pad = 0.1*(qovm_max - qovm_min)
  limits(qovm_min-qovm_pad,qovm_max+qovm_pad,krot_v_min-krot_v_pad,krot_v_max+krot_v_pad) 
  ptitles("Angular Phase Advance Wavenumber: Beam Launch in Bz=0","Q/A","Wavenumber [Rad/m]",)
  fma()

diag_plt_krot_v() 


# Make plot of initial Brho by species 
plt_diag_bro(label = "Initial Rigidity by Species") 

#raise exception("to here")


# Setup diagnostics
# Diagnostics are grouped into several classes:
#   - Particle:  Snapshot plots of distribution function projections 
#   - Field:     Snapshot plots of self fields 
#   - History:   History plots on the evolution of moments and particle counts 
#                   accumulated as the simulation advances.   

# --- set max simulation step for diagnostic setup 
max_diag_step = 1.e10

# --- set history diagnostic and moment accumulations 
ds_diag = 1.*cm 
top.nhist = max(1,nint(ds_diag/wxy.ds))           # step interval for histories 
top.itmomnts[0:3] = [0,max_diag_step,top.nhist]   # do loop ranges for moments 
                                                  #   and status writes to tty

#raise Exception("to here")

# Fix intitial history diagnostics to account for species weight changes
#top.jhist = top.jhist-1   # needed to be minus 1 to reset save in right postion 
#from getzmom import *
#zmmnt() 
#savehist(0.) 


# ---- local diagnostic history arrays 
hl_lenhist_max = 10000
hl_ekin     = fzeros([hl_lenhist_max,top.ns])  # axial beam kinetic energy 
hl_spnum    = fzeros([hl_lenhist_max,top.ns])  # number simulation particles
hl_spnumt   = fzeros([hl_lenhist_max])         # number simulation particles (all species) 
hl_ibeam_p  = fzeros([hl_lenhist_max,top.ns])  # beam current (particle)   
hl_ibeam_e  = fzeros([hl_lenhist_max,top.ns])  # beam current (electrical) 
hl_ibeam_pt = fzeros([hl_lenhist_max])         # total beam current (particle)   
hl_ibeam_et = fzeros([hl_lenhist_max])         # total beam current (electrical)   
hl_lambda_p = fzeros([hl_lenhist_max,top.ns])  # line charge (particle) 
hl_lambda_e = fzeros([hl_lenhist_max,top.ns])  # line charge (electrical) 
hl_pth_bar  = fzeros([hl_lenhist_max,top.ns])  # canonical angular momentum (tilde{P_theta}) 
hl_pthn_bar = fzeros([hl_lenhist_max,top.ns])  # normalized canonical angular momentum (P_theta/(m*c))  
hl_lz_bar   = fzeros([hl_lenhist_max,top.ns])  # mechanical angular momentum
hl_lz_bar2   = fzeros([hl_lenhist_max,top.ns])  # mechanical angular momentum
hl_krot     = fzeros([hl_lenhist_max,top.ns])  # rotation wavenumber  
hl_lang     = fzeros([hl_lenhist_max,top.ns])  # Larmor rotation angle  
hl_epspv    = fzeros([hl_lenhist_max,top.ns])  # rms total phase volume emittance 
hl_epspvn   = fzeros([hl_lenhist_max,top.ns])  # rms normalized total phase volume emittance ** warning save scaled mm-mrad **
#
hl_Rt       = fzeros([hl_lenhist_max])         # rms rad   (all species) 
hl_Qt       = fzeros([hl_lenhist_max])         # Perveance (all species) 
hl_emitt    = fzeros([hl_lenhist_max])         # emittance (all species) 
hl_sovs0    = fzeros([hl_lenhist_max])         # effective SC depression  

hl_dz = top.nhist*wxy.ds

@callfromafterstep
def diag_hist_hl():
  # check step in history accumulation cycle 
  if top.it%top.nhist != 0: return
  # accumulate history diagnostics by species 
  for ii in sp.keys():
    s = sp[ii]
    js = s.js 
    #
    weight = sum(s.sw*s.w) 
    #
    vbeam = sum( (s.sw*s.w)*s.getvz() )/weight
    gammabeam = 1./sqrt(1.-(vbeam/clight)**2)      
    brho  = s.mass*gammabeam*vbeam/s.charge
    #
    #rsq = (s.getx())**2 + (s.gety())**2 # species radii squared for ptheta calculation
    #r   = sqrt(rsq)
    r   = s.getr() 
    rsq = r*r  
    #
    #avg_xsq = sum( (s.sw*s.w)*(s.xp)**2 )/weight
    #avg_ysq = sum( (s.sw*s.w)*(s.yp)**2 )/weight
    avg_rsq = sum( (s.sw*s.w)*rsq )/weight
    #
    #avg_xyp = sum( (s.sw*s.w)*s.getx()*s.getyp() )/weight
    #avg_yxp = sum( (s.sw*s.w)*s.gety()*s.getxp() )/weight
    #
    avg_xpy = s.mass*sum( (s.sw*s.w)*s.getx()*s.getuy() )/weight
    avg_ypx = s.mass*sum( (s.sw*s.w)*s.gety()*s.getux() )/weight
    #
    bz0  = getappliedfields(x=0.,y=0.,z=top.zbeam)[5]
    # 
    # --- Axial kinetic energy 
    hl_ekin[top.jhist,js] = s.mass*clight**2*(gammabeam - 1.)/jperev     
    # --- Simulation Particle Number 
    hl_spnum[top.jhist,js] = s.getn()  
    # --- Current, particle (approx here) 
    hl_ibeam_p[top.jhist,js] = s.charge*s.sw*(s.vbeam0/vbeam)*sum( s.getvz() )  # Fix! use weights correctly  ?? 
    # --- Current, electrical  
    hl_ibeam_e[top.jhist,js] = s.charge*sum( (s.sw*s.w)*s.getvz() )             # Fix! use weights correctly  ??
    # --- line charge 
    hl_lambda_p[top.jhist,js] = hl_ibeam_p[top.jhist,js]/vbeam 
    hl_lambda_e[top.jhist,js] = hl_ibeam_e[top.jhist,js]/vbeam 
    # --- Mechanical angular momentum: <x*y'> - <y*x'>  
    hl_lz_bar2[top.jhist,js] = avg_xpy - avg_ypx
    hl_lz_bar[top.jhist,js] = (avg_xpy - avg_ypx)/(s.mass*gammabeam*vbeam)  
    # --- Normalized canonical angular momentum 
    #hl_pthn_bar[top.jhist,js] = gammabeam*(vbeam/clight)*hl_pth_bar[top.jhist,js]
    #hl_pthn_bar[top.jhist,js] = ( avg_xpy - avg_ypx + (s.charge*bz0/2.)*avg_rsq )/(s.mass*clight)  
    hl_pthn_bar[top.jhist,js] = ( avg_xpy - avg_ypx + sum( (s.sw*s.w)*s.charge*r*getatheta(r) )/weight )/(s.mass*clight)
    # --- Canonical angular momentum (scaled by gamma_b*beta_b*m*c) 
    #hl_pth_bar[top.jhist,js] = hl_lz_bar[top.jhist,js] + bz0/(2.*brho)*avg_rsq 
    hl_pth_bar[top.jhist,js] = hl_pthn_bar[top.jhist,js]/(gammabeam*vbeam/clight) 
    # --- rms total phase volume emittance
    hl_epspv[top.jhist,js] = sqrt( (top.hepsr[0,top.jhist,js])**2 + 4.*(hl_pth_bar[top.jhist,js])**2 ) 
    # --- rms normalized total phase volume emittance ** warning norm emittance scaled mm-mrad to keep to Warp pattern ** 
    hl_epspvn[top.jhist,js] = (gammabeam*(vbeam/clight))*hl_epspv[top.jhist,js]*1.e6 
    # --- Rotation wavenumber 
    hl_krot[top.jhist,js] = hl_lz_bar[top.jhist,js]/avg_rsq
    # --- Larmor Rotation angle: integrate from previous step  
    if top.jhist == 0:
      hl_lang[0,js] = 0.  
    else:
      hl_lang[top.jhist,js] = hl_lang[top.jhist-1,js] + 0.5*hl_dz*(hl_krot[top.jhist-1,js]+hl_krot[top.jhist,js])  
  # --- total number of simulation particles 
  hl_spnumt[top.jhist] = float(sum(hl_spnum[top.jhist,:]))
  # --- total currents 
  hl_ibeam_pt[top.jhist] = sum(hl_ibeam_p[top.jhist,:]) 
  hl_ibeam_et[top.jhist] = sum(hl_ibeam_e[top.jhist,:]) 
  # --- total perveance 
  #hl_Rt[top.jhist] = sqrt(average( getx(jslist=-1)**2 + gety(jslist=-1))) 
 

diag_hist_hl()   # make sure initial diagnostic saved before any steps 


# --- Plot limits for particle phase space plots. If lframe = true (default
#     false) diagnostics such as ppxxp for x-x' particle phase space will
#     use these ranges.  
#      max/min x,y   plot coordinates (m) 
#      max/min x',y' plot coordinates (rad)
l_diag = r_p
top.xplmax =  l_diag  
top.xplmin = -l_diag
top.yplmax =  l_diag
top.yplmin = -l_diag         
top.xpplmax = 75.*mr
top.xpplmin = -top.xpplmax    
top.ypplmax =  top.xpplmax 
top.ypplmin = -top.xpplmax

# --- Color palette for phase-space plots (comment for default)
#     Search for .gp suffix files in the Warp scripts directory for possible
#     choices.  Some useful ones include:
#       earth.gp   (default)        heat.gp     (heat) 
#       gray.gp    (gray scale)     rainbow.gp  (rainbow) 
#palette("heat.gp")

# --- Set a chop factor for particle phase space plots to avoid plotting
#     too many particles (large storage and features will obscure).  Set
#     for approx 10 K particles per species plotted.  
chop_fraction = 10.e3/float(top.npmax) 


# --- Particle phase space diagnostics.
#     The list diag_step_part contains all steps where diagnostics in
#     diag_part() are made.  The list can contain repeated elements
#     and need not be ordered

diag_part_z = array([
  z_launch,
  #s4p1_zc,
  #s4p1_zc-20.*cm,
  #s4p1_zc+20.*cm, 
  s4p2_zc,
  s4p2_zc-20.*cm,
  s4p2_zc+20.*cm, 
  #gag_zc, 
  #gag_col_zs,
  #gag_col_zs-5.*cm,
  #gag_col_ze,
  #gag_col_ze+5.*cm,
  z_adv,
  #(s4p1_zc+gag_zc)/2.,
  (s4p2_zc+z_adv)/2. ]) 

diag_part_z_name = [ 
  "Initial Launch", 
  "S4 Solenoid #1: z-Center", 
  "S4 Solenoid #1: z-Center - 20 cm", 
  "S4 Solenoid #1: z-Center + 20 cm",
  "S4 Solenoid #2: z-Center", 
  "S4 Solenoid #2: z-Center - 20 cm", 
  "S4 Solenoid #2: z-Center + 20 cm",
  "Grated Gap: z-Center",
  #"Grated Gap: z-Start", 
  #"Grated Gap: z-Center - 5 cm", 
  #"Grated Gap: z-End",
  "Grated Gap: z-End + 5 cm", 
  "Final: Before D2 Dipole",
  "Between S4 Solenoid #1 and Grated Gap",
  "Between S4 Solenoid #2 and Final (Before D2 Dipole)" 
                   ]

diag_part_step = nint((diag_part_z-z_launch)/wxy.ds)

diag_part_z_names = {diag_part_step[i]:diag_part_z_name[i] for i in range(len(diag_part_step))}



def diag_part(plt_xy=False,plt_xxp=False,plt_yyp=False,plt_xpyp=False,
              plt_trace=False, plt_denxy=False, plt_denr=False):
  print "Making particle diagnostic plots"
  #
  try:  
    z_label = diag_part_z_names[top.it]
  except:
    z_label = ""
  #
  # --- x-y projection
  if plt_xy:
    # --- All Species 
    #  Caution:  js=-1 with density plot will just overlay species contour plots 
    #ppxy(js=-1,lframe=true,chopped=chop_fraction,color='density',ncolor=25,
    #     titles=false,yscale=1./mm,xscale=1./mm)
    ppxy(js=-1,lframe=true,chopped=chop_fraction,titles=false,yscale=1./mm,xscale=1./mm)
    ptitles("x-y Phase Space: All Species, z = %5.2f m"%(top.zbeam),
            "x [mm]","y [mm]",z_label)
    fma()
    # --- Target Species 
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      co = s.color
      lab+= ii + "("+co+"), "
      s.ppxy(lframe=true,chopped=chop_fraction,titles=false,yscale=1./mm,xscale=1./mm)
    ptitles("x-y Phase Space: "+lab+" z = %5.2f m"%(top.zbeam),"x [mm]","y [mm]",z_label)
    fma()
  # --- x-x' projection
  if plt_xxp: 
    # --- All Species
    #   Caution:  js = -1 with density plot will overlay species contour plots  
    #ppxxp(js = -1,lframe=true,chopped=chop_fraction,slope='auto',color='density',ncolor=25,
    #      titles=false,yscale=1./mr,xscale=1./mm)
    ppxxp(js = -1,lframe=true,chopped=chop_fraction,slope='auto',titles=false,yscale=1./mr,xscale=1./mm)
    ptitles("x-x' Phase Space: All Species, z = %5.2f m"%(top.zbeam),"x [mm]","x' [mrad]",z_label)
    fma()
    # --- Target Species 
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      co = s.color
      lab+= ii + "("+co+"), "
      s.ppxxp(lframe=true,chopped=chop_fraction,slope='auto',titles=false,yscale=1./mr,xscale=1./mm)
    ptitles("x-x' Phase Space: "+lab+" z = %5.2f m"%(top.zbeam),"x [mm]","x' [mrad]",z_label)
    fma()
  # --- y-y' projection
  if plt_yyp:
    # --- All Species 
    #   Caution: js=-1 with denisty plot will overlay species contour plots 
    #ppyyp(js=-1,lframe=true,chopped=chop_fraction,slope='auto',color='density',ncolor=25,
    #      titles=false,yscale=1./mr,xscale=1./mm)
    ppyyp(js=-1,lframe=true,chopped=chop_fraction,slope='auto',color='density',ncolor=25,
          titles=false,yscale=1./mr,xscale=1./mm)
    ptitles("y-y' Phase Space: All Species, z = %5.2f m"%(top.zbeam),
            "y [mm]","y' [mrad]",z_label)
    fma()
    # --- Target Species 
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      co = s.color
      lab+= ii + "("+co+"), "
      s.ppyyp(lframe=true,chopped=chop_fraction,slope='auto',titles=false,yscale=1./mr,xscale=1./mm)
    ptitles("y-y' Phase Space: "+lab+" z = %5.2f m"%(top.zbeam),"y [mm]","y' [mrad]",z_label)
    fma()
  # --- x'-y' projection
  if plt_xpyp:
    # --- All Species 
    #   Caution:  js=-1 with density plot will overlay species countours 
    #ppxpyp(js=-1,lframe=true,chopped=chop_fraction,slope='auto',color='density',ncolor=25,
    #       titles=false,yscale=1./mr,xscale=1./mr)
    ppxpyp(js=-1,lframe=true,chopped=chop_fraction,slope='auto',titles=false,yscale=1./mr,xscale=1./mr)
    ptitles("x'-y' Phase Space: All Species, z = %5.2f m"%(top.zbeam),"x' [mrad]","y' [mrad]",z_label)
    fma()
    # --- Target Species 
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      co = s.color
      lab+= ii + "("+co+"), "
      s.ppxpyp(lframe=true,chopped=chop_fraction,slope='auto',titles=false,yscale=1./mr,xscale=1./mm)
    ptitles("x'-y' Phase Space: "+lab+" z = %5.2f m"%(top.zbeam),"x' [mrad]","y' [mrad]",z_label)
    fma()
  # --- x-y, x-x', y-y', x'-y' projections, 4 to a page (trace-space)
  if plt_trace:
    # --- All Species 
    pptrace(lframe=true,chopped=chop_fraction,slope='auto',color='density',ncolor=25)
    fma()
  # --- charge density on x and y axes
  if plt_denxy:
    rho_sc = 1.
    ix_cen = sum(where(w3d.xmesh < 0.,1,0))
    iy_cen = sum(where(w3d.ymesh < 0.,1,0))
    # --- All Species 
    rho_x = getrho(iy=iy_cen)
    rho_y = getrho(ix=ix_cen) 
    # 
    plg(rho_x/rho_sc,w3d.xmesh/mm)
    if w3d.l4symtry: plg(rho_x/rho_sc,-w3d.xmesh/mm) 
    plg(rho_y/rho_sc,w3d.ymesh/mm,color="red")
    if w3d.l4symtry or w3d.l2symtry: 
      plg(rho_y/rho_sc,-w3d.ymesh/mm,color="red")
    ptitles("Charge Density: All Species, on x[b], y[r] Axes: z = %5.2f m"%(top.zbeam),
            "x,y [mm]","Density [arb units]",z_label)
    fma()
    # --- Target Species: species.get_density() returns density     
    for ii in sp_target:
      s = sp[ii]
      co = s.color
      den = s.get_density()/cm**3
      plg(den[:,iy_cen],w3d.xmesh/mm)
      if w3d.l4symtry: plg(den[:,iy_cen],-w3d.xmesh/mm) 
      plg(den[ix_cen,:],w3d.ymesh/mm,color="red")
      if w3d.l4symtry or w3d.l2symtry: plg(den[ix_cen,:],-w3d.ymesh/mm,color="red")
      ptitles("Density: "+ii+" on x[b], y[r] Axes: z = %5.2f m"%(top.zbeam),
              "x,y [mm]","Density [#/cm^3]",z_label)
      fma()
  # --- charge density on radial mesh 
  if plt_denr:
    # --- radial mesh reflecting x-y grid structure to illustrate simulation noise
    nr    = nint(sqrt(w3d.nx/(2.*sym_x)*w3d.ny/(2.*sym_y)))
    rmax  = sqrt(w3d.xmmax*w3d.ymmax)
    dr    = rmax/nr 
    rmesh = linspace(0.,rmax,num=nr+1)
    #
    sp_list = sp_target #+ ["All"] 
    ns   = len(sp_list) 
    # --- density as a function or r on mesh array 
    den  = zeros(nr+1)
    #   
    weightr = zeros(nr+1)   
    count   = zeros(nr+1)   
    # --- for all species on mesh 
    for ii in sp.keys():
       s  = sp[ii]
       #
       np = s.getn() 
       rp = s.getr() 
       wp = s.getweights()
       #
       deposgrid1d(1,np,rp,wp,nr,weightr,count,0.,rmax)
    #
    den[1:nr+1] = weightr[1:nr+1]/(2.*pi*dr*rmesh[1:nr+1])
    den[0]      = den[1]   # set origin by next grid up to remove distraction
    # 
    plg(den/cm**3, rmesh/mm)
    plg(den/cm**3,-rmesh/mm) 
    ptitles("Radial Number Density: All Species, z = %5.2f m"%(top.zbeam),"radius r [mm]","rho [particles/cm**3]",z_label)
    ir = min(nr,sum(where(den>0,1,0)))      # index farthest radial extent of rho in radial mesh assuming no halo  
    rmmax = max(1.2*rmesh[ir],0.01) # set curoff to contain radial density  
    rmmax = cm*nint(rmmax/cm + 0.5) # round up to nearest cm to contain plot 
    denmax = 1.2*maxnd(den) 
    limits(-rmmax/mm,rmmax/mm,0.,denmax/cm**3)
    fma() 
    # --- for target species on mesh 
    for ii in sp_target:
       s  = sp[ii]
       co = s.color 
       lab = ii + "("+co+"), "
       #
       np = s.getn() 
       rp = s.getr() 
       wp = s.getweights()
       #
       weightr = zeros(nr+1)   # reset for clean accumulation/count with itask = 1 
       count   = zeros(nr+1)   
       deposgrid1d(1,np,rp,wp,nr,weightr,count,0.,rmax)
       # 
       den[1:nr+1] = weightr[1:nr+1]/(2.*pi*dr*rmesh[1:nr+1])
       den[0]      = den[1]   # set origin by next grid up to remove distraction
       # 
       plg(den/cm**3, rmesh/mm,color=co)
       plg(den/cm**3,-rmesh/mm,color=co) 
       ptitles("Radial Number Density: "+lab+" z = %5.2f m"%(top.zbeam),"radius r [mm]","rho [particles/cm**3]",z_label)
       ir = sum(where(den>0,1,0))      # index farthest radial extent of rho in radial mesh assuming no halo  
       rmmax = max(1.2*rmesh[ir],0.01) # set curoff to contain radial density  
       rmmax = cm*nint(rmmax/cm + 0.5) # round up to nearest cm to contain plot 
       denmax = 1.2*maxnd(den) 
       limits(-rmmax/mm,rmmax/mm,0.,denmax/cm**3)
       fma() 
  

# --- Field diagnostics.  
#     The list diag_step_field containins all steps where
#     diagnostics in diag_field() are made. The list can contain repeated
#     elements and need not be ordered.   

diag_field_z = array([
  z_launch,
  #s4p1_zc,gag_zc,
  s4p2_zc,
  z_adv 
                    ]) 

diag_field_z_name = [ 
  "Initial Launch", 
  "S4 Solenoid #1: z-Center", 
  "Grated Gap: z-Center",
  "S4 Solenoid #1: z-Center", 
  "Final: Before D2 Dipole"
                     ]

diag_field_step = nint((diag_field_z-z_launch)/wxy.ds)

diag_field_z_names = {diag_field_step[i]:diag_field_z_name[i] for i in range(len(diag_field_step))}



def diag_field(plt_pa=False,plt_pc=False,plt_pc_xy=False):
  print "Making field diagnostic plots"
  #
  try:  
    z_label = diag_field_z_names[top.it]
  except:
    z_label = ""
  # --- self-field electrostatic potential
  if plt_pc:
    pfxy(cond=true,titles=false,yscale=1./mm,xscale=1./mm,iz = 0)
    ptitles("Self-Field Potential: z = %5.2f"%(top.zbeam),
            "x [mm]","y [mm]",z_label)
    fma()
  # --- self-field electrostatic potential and particles together
  if plt_pc_xy:
    # --- All particle species included 
    pfxy(cond=true,titles=false,yscale=1./mm,xscale=1./mm)
    #   Caution: js=-1 with density plot will superimpose species contours 
    #ppxy(js=-1,lframe=true,chopped=chop_fraction,color='density',ncolor=25,
    #     titles=false,yscale=1./mm,xscale=1./mm)
    ppxy(js=-1,lframe=true,chopped=chop_fraction,titles=false,yscale=1./mm,xscale=1./mm)
    ptitles("Self-Field Potential: z = %5.2f"%(top.zbeam),
            "x [mm]","y [mm]",z_label)
    fma()
    # --- Target particle species 
    lab = ""
    pfxy(cond=true,titles=false,yscale=1./mm,xscale=1./mm)
    for ii in sp_target:
      s = sp[ii]
      co = s.color
      lab+= ii + "("+co+"), "
      s.ppxy(lframe=true,chopped=chop_fraction,titles=false,yscale=1./mm,xscale=1./mm)
      s.ppxy(lframe=true,chopped=chop_fraction,titles=false,yscale=1./mm,xscale=1./mm)
    ptitles("Self-Field Potential: + "+lab+" Particles, z = %5.2f"%(top.zbeam),"x [mm]","y [mm]",z_label)
    fma()
  # --- Electrostatic potential on principal axes 
  if plt_pa:
    diag_plt_phi_ax(label="Beam Potential along y,x = 0 [b,r] at z = %5.2f"%(top.zbeam))
    fma()
    # 
    xrms = max(top.xrms[0,sp['U33'].js],top.xrms[0,sp['U34'].js]) 
    diag_plt_phi_ax(label="Beam Potential along y,x = 0 [b,r] at z = %5.2f"%(top.zbeam),xmax=2.*xrms) 
    fma() 


# --- History diagnostics.  These can be made at intermediate stages of the
#     run as well as at the end.  The list diag_step_hist contains all
#     steps where diagnostics in diag_hsit() are made. The list can
#     contain repeated elements and need not be ordered.
#     Notes:
#      * Many additional history diagnostics can be added by looking for
#        relevant moments accumulated in the Warp (see the variable group
#        "Hist" in top.v for an extensive list of variables that can be
#         used) and using gist commands to make relevant plots

diag_hist_z    = array([z_adv]) #array([gag_col_zs,z_adv])
diag_hist_step = nint((diag_hist_z-z_launch)/wxy.ds)

def diag_hist(plt_ekin=False,plt_spnum=False,plt_curr_p=False,plt_curr_e=False,plt_lam_p=False,plt_lam_e=False,
              plt_lz=False,plt_pth=False,plt_pthn=False,plt_krot=False,plt_lang=False,
              plt_cen=False,plt_envrms=False,plt_envmax=False,plt_envrmsp=False, 
              plt_emit=False,plt_emitn=False,plt_emitg=False,plt_emitng=False,plt_emitr=False,plt_emitnr=False,
              plt_emitpv=False,emitpvn=False, 
              plt_temp=False,plt_temp_flow=False):
  print "Making history diagnostic plots"
  #
  # --- kinetic energy 
  if plt_ekin:
    # --- All Species Combined, MeV
    #hpekin(titles=false,yscale=1.,lhzbeam=true)
    #ptitles("History: All Species Kinetic Energy","z [m]","MeV", )
    #fma()
    # --- All Species, in keV/u 
    for ii in sort(sp.keys()):
      s = sp[ii]
      js = s.js
      co = s.color
      A  = s.mass/amu
      plg(hl_ekin[0:top.jhist+1,js]/(A*kV),top.hzbeam[0:top.jhist+1],color=co)        
      #hpekin(js=js,color=co,titles=false,yscale=1./A,lhzbeam=true)    
    ptitles("History: Kinetic Energy","z [m]","KeV/u", )
    fma()
    # --- U species, in MeV/u
    for ii in sort(sp_U.keys()):
      s = sp[ii]
      js = s.js
      co = s.color
      A  = s.mass/amu
      #hpekin(js=js,color=co,titles=false,yscale=1./A,lhzbeam=true)
      plg(hl_ekin[0:top.jhist+1,js]/(A*kV),top.hzbeam[0:top.jhist+1],color=co)        
    ptitles("History: U Species Kinetic Energy","z [m]","KeV/u", )
    fma()
    # --- O species, in MeV/u 
    for ii in sort(sp_O.keys()):
      s = sp[ii]
      js = s.js
      co = s.color
      A  = s.mass/amu
      plg(hl_ekin[0:top.jhist+1,js]/(A*kV),top.hzbeam[0:top.jhist+1],color=co)        
      #hpekin(js=js,color=co,titles=false,yscale=1./A,lhzbeam=true) # Was getting wrong answer !!
    ptitles("History: O Species Kinetic Energy","z [m]","KeV/u", )
    fma()
    # --- By Target Species, in kV/Q
    #     Plot by KV/Q so you can see total potential gain falling through 
    #     full bias to check system tuning  
    zi = top.hzbeam[0]
    zf = top.hzbeam[top.jhist]
    ekin_t = Bias/kV
    lab = ""
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      Q  = s.charge_state
      lab+= ii + "("+co+"), "
      plg(hl_ekin[0:top.jhist+1,js]/(Q*kV),top.hzbeam[0:top.jhist+1],color=co)        
      #hpekin(js=js,color=co,titles=false,yscale=1./Q,lhzbeam=true)
    plg(array([ekin_t,ekin_t]),array([zi,zf]),type="dash") 
    ptitles("History: "+lab+"Kinetic Energy","z [m]","KeV/Q", )
    limits(zi,zf,0.,1.2*ekin_t) 
    fma()
  # --- simulation particle number (to check for lost particles)
  #     Comment: tried using hppnum() but was unclear what was being plotted 
  if plt_spnum:
    # --- All Species Combined  
    plg(hl_spnumt[0:top.jhist+1],top.hzbeam[0:top.jhist+1])    
    ptitles("History: Live Sim Particle Number (all species)", "z [m]","Particle Number (simulation)", )
    fma()
    # --- All Species Individually 
    for ii in sort(sp.keys()):
      s = sp[ii]
      js = s.js
      co = s.color
      plg(hl_spnum[0:top.jhist+1,js],top.hzbeam[0:top.jhist+1],color=co)        
    ptitles("History: Live Sim Particle Number (by species)","z [m]","Particle Number (simulation)", )
    fma() 
    # --- Target Species
    lab = "" 
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg(hl_spnum[0:top.jhist+1,js],top.hzbeam[0:top.jhist+1],color=co)        
    ptitles("History: "+lab+" Live Sim Particle Number","z [m]","Particle Number (simulation)", )
    fma()             
  # --- current (particle)  
  if plt_curr_p:
    # --- All Species Combined  
    for ii in sort(sp.keys()):
      s = sp[ii]        
      js = s.js
      co = s.color
      plg(hl_ibeam_p[0:top.jhist+1,js]*1.e6,top.hzbeam[0:top.jhist+1],color=co)    
    ptitles("History: Species Particle Current (approx)", "z [m]","Current (microA)", )
    fma() 
    # --- Target Species 
    lab = "" 
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg(hl_ibeam_p[0:top.jhist+1,js]*1.e6,top.hzbeam[0:top.jhist+1],color=co)    
    ptitles("History: "+lab+" Particle Current (approx)","z [m]","Current (microA)", )
    fma() 
    # --- Total
    plg(hl_ibeam_pt[0:top.jhist+1]*1.e3,top.hzbeam[0:top.jhist+1])    
    ptitles("History: Total Particle Current (approx)","z [m]","Current (mA)", )
    fma()             
  # --- current (electrical)  
  if plt_curr_e:
    # --- All Species Combined  
    for ii in sort(sp.keys()):
      s = sp[ii]        
      js = s.js
      co = s.color
      plg(hl_ibeam_e[0:top.jhist+1,js]*1.e6,top.hzbeam[0:top.jhist+1],color=co)    
    ptitles("History: Species Electrical Current", "z [m]","Current (microA)", )
    fma()
    # --- Target Species
    lab = "" 
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg(hl_ibeam_e[0:top.jhist+1,js]*1.e6,top.hzbeam[0:top.jhist+1],color=co)    
    ptitles("History: "+lab+" Electrical Current","z [m]","Current (microA)", )
    fma()
    # --- Total
    plg(hl_ibeam_et[0:top.jhist+1]*1.e3,top.hzbeam[0:top.jhist+1])    
    ptitles("History: Total Electrical Current","z [m]","Current (mA)", )
    fma()                          
  # --- line charge (particle)  
  if plt_lam_p:
    # --- All Species Combined  
    for ii in sort(sp.keys()):
      s = sp[ii]        
      js = s.js
      co = s.color
      plg(hl_lambda_p[0:top.jhist+1,js]*10**9,top.hzbeam[0:top.jhist+1],color=co)    
    ptitles("History: Species Particle Line Charge", "z [m]","Line Charge (nC/m)", )
    fma()
    # --- Target Species
    lab = "" 
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg(hl_lambda_p[0:top.jhist+1,js]*10**9,top.hzbeam[0:top.jhist+1],color=co)    
    ptitles("History: "+lab+" Particle Line Charge","z [m]","Line Charge (nC/m)", )
    fma()             
  # --- line charge (electrical)  
  if plt_lam_e:
    # --- All Species Combined  
    for ii in sort(sp.keys()):
      s = sp[ii]        
      js = s.js
      co = s.color
      plg(hl_lambda_e[0:top.jhist+1,js]*10**9,top.hzbeam[0:top.jhist+1],color=co)    
    ptitles("History: Species Electrical Line Charge", "z [m]","Line Charge (nC/m)", )
    fma()
    # --- Target Species
    lab = "" 
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg(hl_lambda_e[0:top.jhist+1,js]*10**9,top.hzbeam[0:top.jhist+1],color=co)    
    ptitles("History: "+lab+" Electrical Line Charge","z [m]","Line Charge (nC/m)", )
    fma()             
  # --- lz mechanical angular momentum  
  if plt_lz:
    # --- All Species Combined  
    for ii in sort(sp.keys()):
      s = sp[ii]        
      js = s.js
      co = s.color
      plg(hl_lz_bar[0:top.jhist+1,js]*10**6,top.hzbeam[0:top.jhist+1],color=co)    
    ptitles("History: Species Mechanical Angular Mom", "z [m]","<xy'>-<yx'>  [mm-mrad]", )
    fma()
    # --- Target Species
    lab = "" 
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg(hl_lz_bar[0:top.jhist+1,js]*10**6,top.hzbeam[0:top.jhist+1],color=co)    
    ptitles("History: "+lab+" Mechanical Angular Mom","z [m]","<xy'>-<yx'>  [mm-mrad]", )
    fma()             
  # --- canonical angular momentum  
  if plt_pth:
    # --- All Species Combined  
    for ii in sort(sp.keys()):
      s = sp[ii]        
      js = s.js
      co = s.color
      plg(hl_pth_bar[0:top.jhist+1,js]*10**6,top.hzbeam[0:top.jhist+1],color=co)    
    ptitles("History: Species Canonical Angular Mom", "z [m]","Ptheta  [mm-mrad]", )
    fma()
    # --- Target Species
    lab = "" 
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg(hl_pth_bar[0:top.jhist+1,js]*10**6,top.hzbeam[0:top.jhist+1],color=co)    
    ptitles("History: "+lab+" Canonical Angular Mom","z [m]","Ptheta  [mm-mrad]", )
    fma()             
  # --- canonical angular momentum (normalized)   
  if plt_pthn:
    # --- All Species Combined  
    for ii in sort(sp.keys()):
      s = sp[ii]        
      js = s.js
      co = s.color
      plg(hl_pthn_bar[0:top.jhist+1,js]*10**6,top.hzbeam[0:top.jhist+1],color=co)    
    ptitles("History: Species Norm Canonical Angular Mom", "z [m]","Ptheta  [mm-mrad]", )
    fma()
    # --- Target Species
    lab = "" 
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg(hl_pthn_bar[0:top.jhist+1,js]*10**6,top.hzbeam[0:top.jhist+1],color=co)    
    ptitles("History: "+lab+" Norm Canonical Angular Mom","z [m]","Ptheta  [mm-mrad]", )
    fma()             
  # --- effective rotation wavenumber 
  if plt_krot:
    # --- All Species Combined  
    for ii in sort(sp.keys()):
      s = sp[ii]        
      js = s.js
      co = s.color
      plg(hl_krot[0:top.jhist+1,js],top.hzbeam[0:top.jhist+1],color=co)    
    ptitles("History: Species Effective Rot Wavenumber", "z [m]","krot  [rad/m]", )
    fma()
    # --- Target Species
    lab = "" 
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg(hl_krot[0:top.jhist+1,js],top.hzbeam[0:top.jhist+1],color=co)    
    ptitles("History: "+lab+" Effective Rot Wavenumber","z [m]","krot  [rad/m]", )
    fma()             
  # --- larmor rotation angle   
  if plt_lang:
    # --- All Species Combined  
    for ii in sort(sp.keys()):
      s = sp[ii]        
      js = s.js
      co = s.color
      plg((180./pi)*hl_lang[0:top.jhist+1,js],top.hzbeam[0:top.jhist+1],color=co)    
    ptitles("History: Species Larmor Rot Angle", "z [m]","Rotation [deg]", )
    fma()
    # --- Target Species
    lab = "" 
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg((180./pi)*hl_lang[0:top.jhist+1,js],top.hzbeam[0:top.jhist+1],color=co)    
    ptitles("History: "+lab+" Larmor Rot Angle","z [m]","Rotation  [deg]", )
    fma()             
  # --- centroid
  if plt_cen:
    # All Species Combined, x- and y-plane 
    hpxbar(titles=false,yscale=1./mm,lhzbeam=true)
    hpybar(titles=false,yscale=1./mm,lhzbeam=true,color="red")
    ptitles("History: All Species x-,y-Centroid: x[b], y[r]","z [m]","<x>, <y> Centroids [mm]", )
    fma()
    # --- By Target Species, x-plane 
    hpxbar(titles=false,yscale=1./(sqrt(2.)*mm),lhzbeam=true)
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpxbar(js=js,color=co,titles=false,yscale=1./(sqrt(2.)*mm),lhzbeam=true)    
    ptitles("History: "+lab+"x-Centroid","z [m]","<x> [mm]", )
    fma()
    #
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpxbar(js=js,color=co,titles=false,yscale=1./(sqrt(2.)*mm),lhzbeam=true)    
    ptitles("History: "+lab+"x-Centroid","z [m]","<x> [mm]", )
    fma()
    # --- By Target Species, y-plane 
    hpybar(titles=false,yscale=1./(sqrt(2.)*mm),lhzbeam=true)
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpybar(js=js,color=co,titles=false,yscale=1./(sqrt(2.)*mm),lhzbeam=true)    
    ptitles("History: "+lab+"y-Centroid","z [m]","<y> [mm]", )
    fma()
    #
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpybar(js=js,color=co,titles=false,yscale=1./(sqrt(2.)*mm),lhzbeam=true)    
    ptitles("History: "+lab+"y-Centroid","z [m]","<y> [mm]", )
    fma()
  # --- rms envelope width 
  if plt_envrms:
    # --- All Species Combined, x- and y-plane  
    hpenvx(titles=false,yscale=1./(2.*mm),lhzbeam=true)    
    hpenvy(titles=false,yscale=1./(2.*mm),lhzbeam=true,color="red")
    ptitles("History: All Species RMS Envelope: x[b], y[r]","z [m]","RMS Width [mm]", )
    fma()
    # --- Target Species, x-plane 
    hpenvx(titles=false,yscale=1./(2.*mm),lhzbeam=true)
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpenvx(js=js,color=co,titles=false,yscale=1./(2.*mm),lhzbeam=true)    
    ptitles("History: "+lab+"RMS x-Envelope","z [m]","RMS Width [mm]", )
    fma()
    #
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpenvx(js=js,color=co,titles=false,yscale=1./(2.*mm),lhzbeam=true)    
    ptitles("History: "+lab+"RMS x-Envelope","z [m]","RMS Width [mm]", )
    fma()
    # --- Target Species, y-plane 
    hpenvy(titles=false,yscale=1./(2.*mm),lhzbeam=true)
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpenvy(js=js,color=co,titles=false,yscale=1./(2.*mm),lhzbeam=true)    
    ptitles("History: "+lab+"RMS y-Envelope","z [m]","RMS Width [mm]", )
    fma()
    #
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpenvy(js=js,color=co,titles=false,yscale=1./(2.*mm),lhzbeam=true)    
    ptitles("History: "+lab+"RMS y-Envelope","z [m]","RMS Width [mm]", )
    fma()
  # --- max particle envelopes 
  if plt_envmax:
    # --- x-plane, All Species Combined  
    for ii in sort(sp.keys()):
      s = sp[ii]        
      js = s.js
      co = s.color
      plg(top.hxmaxp[0:top.jhist+1,js]/mm,top.hzbeam[0:top.jhist+1],color=co)    
    ptitles("History: Species max particle x", "z [m]","Max x [mm]", )
    fma()
    # --- x-plane, Target Species
    lab = "" 
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg(top.hxmaxp[0:top.jhist+1,js]/mm,top.hzbeam[0:top.jhist+1],color=co)    
    ptitles("History: "+lab+" max particle x","z [m]","Max x [mm]", )
    fma()             
    # --- y-plane, All Species Combined  
    for ii in sort(sp.keys()):
      s = sp[ii]        
      js = s.js
      co = s.color
      plg(top.hymaxp[0:top.jhist+1,js]/mm,top.hzbeam[0:top.jhist+1],color=co)    
    ptitles("History: Species max particle y", "z [m]","Max y [mm]", )
    fma()
    # --- y-plane, Target Species
    lab = "" 
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg(top.hymaxp[0:top.jhist+1,js]/mm,top.hzbeam[0:top.jhist+1],color=co)    
    ptitles("History: "+lab+" max particle y","z [m]","Max y [mm]", )
    fma()             
  # --- rms envelope angle  
  if plt_envrmsp:
    # --- Target Species, x-plane 
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg(top.hxxpbar[0,0:top.jhist+1,js]/(top.hxrms[0,0:top.jhist+1,js]*mr),top.hzbeam[0:top.jhist+1],color=co)    
    ptitles("History: "+lab+"RMS x-Envelope Angle","z [m]","RMS Angle [mr]", )
    fma()
    # --- Target Species, y-plane 
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg(top.hyypbar[0,0:top.jhist+1,js]/(top.hyrms[0,0:top.jhist+1,js]*mr),top.hzbeam[0:top.jhist+1],color=co)    
    ptitles("History: "+lab+"RMS y-Envelope Angle","z [m]","RMS Angle [mr]", )
    fma()
  # --- emittance, unnormalized 
  if plt_emit:
    # --- All Species Combined, x- and y-plane 
    hpepsx(titles=false,yscale=1./(mm*mr),lhzbeam=true)
    hpepsy(titles=false,yscale=1./(mm*mr),lhzbeam=true,color="red")
    ptitles("History: All Species RMS Edge x-, y-Emittance: x[b],y[r]","z [m]","Emittance [mm-mr]", )
    fma()
    # --- Target Species, x-plane 
    hpepsx(titles=false,yscale=1./(mm*mr),lhzbeam=true)
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsx(js=js,color=co,titles=false,yscale=1./(mm*mr),lhzbeam=true)    
    ptitles("History: "+lab+"RMS Edge x-Emittance","z [m]","Emittanace [mm-mr]", )
    fma()
    #
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsx(js=js,color=co,titles=false,yscale=1./(mm*mr),lhzbeam=true)    
    ptitles("History: "+lab+"RMS Edge x-Emittance","z [m]","Emittanace [mm-mr]", )
    fma()
    # --- Target Species, y-plane 
    hpepsy(titles=false,yscale=1./(mm*mr),lhzbeam=true)
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsy(js=js,color=co,titles=false,yscale=1./(mm*mr),lhzbeam=true)    
    ptitles("History: "+lab+"RMS Edge y-Emittance","z [m]","Emittance [mm-mr]", )
    fma()
    #
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsy(js=js,color=co,titles=false,yscale=1./(mm*mr),lhzbeam=true)    
    ptitles("History: "+lab+"RMS Edge y-Emittance","z [m]","Emittance [mm-mr]", )
    fma()
  # --- emittance, normalized ** warning norm emittance scaled mm-mrad by default **
  if plt_emitn:
    # --- All Species Combined, x- and y-plane 
    hpepsnx(titles=false,yscale=1.,lhzbeam=true)
    hpepsny(titles=false,yscale=1.,lhzbeam=true,color="red")
    ptitles("History: All Species Norm RMS Edge x-, y-Emittance: x[b],y[r]","z [m]","Norm Emittance [mm-mr]", )
    fma()
    # --- By Target Species, x-plane 
    hpepsnx(titles=false,yscale=1.,lhzbeam=true)
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsnx(js=js,color=co,titles=false,yscale=1.,lhzbeam=true)    
    ptitles("History: "+lab+"Norm RMS Edge x-Emittance","z [m]","Norm Emittanace [mm-mr]", )
    fma()
    #
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsnx(js=js,color=co,titles=false,yscale=1.,lhzbeam=true)    
    ptitles("History: "+lab+"Norm RMS Edge x-Emittance","z [m]","Norm Emittanace [mm-mr]", )
    fma()
    # --- By Target Species, y-plane 
    hpepsny(titles=false,yscale=1.,lhzbeam=true)
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsny(js=js,color=co,titles=false,yscale=1.,lhzbeam=true)    
    ptitles("History: "+lab+"Norm RMS Edge y-Emittance","z [m]","Norm Emittance [mm-mr]", )
    fma()
    #
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsny(js=js,color=co,titles=false,yscale=1.,lhzbeam=true)    
    ptitles("History: "+lab+"Norm RMS Edge y-Emittance","z [m]","Emittance [mm-mr]", )
    fma()
  # --- emittance, generalized unnormalized 
  if plt_emitg:
    # --- All Species Combined, g- and h-plane 
    hpepsg(titles=false,yscale=1./(mm*mr),lhzbeam=true)
    hpepsh(titles=false,yscale=1./(mm*mr),lhzbeam=true,color="red")
    ptitles("History: All Species RMS g-, h-Emittance: g[b],h[r]","z [m]","Emittance [mm-mr]", )
    fma()
    # --- By Target Species, g-plane 
    hpepsg(titles=false,yscale=1./(mm*mr),lhzbeam=true)
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsg(js=js,color=co,titles=false,yscale=1./(mm*mr),lhzbeam=true)    
    ptitles("History: "+lab+"RMS g-Emittance","z [m]","Emittanace [mm-mr]", )
    fma()
    #
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsg(js=js,color=co,titles=false,yscale=1./(mm*mr),lhzbeam=true)    
    ptitles("History: "+lab+"RMS g-Emittance","z [m]","Emittanace [mm-mr]", )
    fma()
    # --- By Target Species, h-plane 
    hpepsh(titles=false,yscale=1./(mm*mr),lhzbeam=true)
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsh(js=js,color=co,titles=false,yscale=1./(mm*mr),lhzbeam=true)    
    ptitles("History: "+lab+"RMS h-Emittance","z [m]","Emittance [mm-mr]", )
    fma()
    #
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsh(js=js,color=co,titles=false,yscale=1./(mm*mr),lhzbeam=true)    
    ptitles("History: "+lab+"RMS h-Emittance","z [m]","Emittance [mm-mr]", )
    fma()
  # --- emittance, generalized normalized 
  if plt_emitng:
    # --- All Species Combined, g- and h-plane 
    hpepsng(titles=false,yscale=1.,lhzbeam=true)
    hpepsnh(titles=false,yscale=1.,lhzbeam=true,color="red")
    ptitles("History: All Species RMS Norm g-, h-Emittance: g[b],h[r]","z [m]","Norm Emittance [mm-mr]", )
    fma()
    # --- By Target Species, g-plane  
    hpepsng(titles=false,yscale=1.,lhzbeam=true)
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsng(js=js,color=co,titles=false,yscale=1.,lhzbeam=true)    
    ptitles("History: "+lab+"RMS Norm g-Emittance","z [m]","Norm Emittanace [mm-mr]", )
    fma()
    #
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsng(js=js,color=co,titles=false,yscale=1.,lhzbeam=true)    
    ptitles("History: "+lab+"RMS Norm g-Emittance","z [m]","Norm Emittanace [mm-mr]", )
    fma()
    # --- By Target Species, h-plane 
    hpepsnh(titles=false,yscale=1.,lhzbeam=true)
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsnh(js=js,color=co,titles=false,yscale=1.,lhzbeam=true)    
    ptitles("History: "+lab+"RMS Norm h-Emittance","z [m]","Norm Emittance [mm-mr]", )
    fma()
    #
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsnh(js=js,color=co,titles=false,yscale=1.,lhzbeam=true)    
    ptitles("History: "+lab+"RMS Norm h-Emittance","z [m]","Norm Emittance [mm-mr]", )
    fma()
  # --- emittance, generalized radial unnormalized 
  if plt_emitr:
    # --- All Species Combined
    hpepsr(titles=false,yscale=1./(mm*mr),lhzbeam=true)
    ptitles("History: All Species RMS r-Emittance","z [m]","Emittance [mm-mr]", )
    fma()
    # --- By Target Species  
    hpepsr(titles=false,yscale=1./(mm*mr),lhzbeam=true)
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsr(js=js,color=co,titles=false,yscale=1./(mm*mr),lhzbeam=true)    
    ptitles("History: "+lab+"RMS r-Emittance","z [m]","Emittanace [mm-mr]", )
    fma()
    #
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsr(js=js,color=co,titles=false,yscale=1./(mm*mr),lhzbeam=true)    
    ptitles("History: "+lab+"RMS r-Emittance","z [m]","Emittanace [mm-mr]", )
    fma()
  # --- emittance, generalized radial normalized ** warning norm emittance scaled mm-mrad by default **
  if plt_emitnr:
    # --- All Species Combined
    hpepsnr(titles=false,yscale=1.,lhzbeam=true)
    ptitles("History: All Species Norm RMS r-Emittance","z [m]","Norm Emittance [mm-mr]", )
    fma()
    # --- By Target Species  
    hpepsnr(titles=false,yscale=1.,lhzbeam=true)
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsnr(js=js,color=co,titles=false,yscale=1.,lhzbeam=true)    
    ptitles("History: "+lab+"RMS Norm r-Emittance","z [m]","Norm Emittanace [mm-mr]", )
    fma()
    #
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsnr(js=js,color=co,titles=false,yscale=1.,lhzbeam=true)    
    ptitles("History: "+lab+"RMS Norm r-Emittance","z [m]","Norm Emittanace [mm-mr]", )
    fma()
  # --- emittance, total phase volume, unnormalized 
  if plt_emitpv:
    # --- All Species Combined  
    for ii in sort(sp.keys()):
      s = sp[ii]        
      js = s.js
      co = s.color
      plg(hl_epspv[0:top.jhist+1,js]/(mm*mr),top.hzbeam[0:top.jhist+1],color=co)    
    ptitles("History: Species Total Phase Volume Emittance", "z [m]","Emittance [mm-mrad]", )
    fma()
    # --- Target Species
    lab = "" 
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg(hl_epspv[0:top.jhist+1,js]/(mm*mr),top.hzbeam[0:top.jhist+1],color=co)    
    ptitles("History: "+lab+" Total Phase Volume Emittance","z [m]","Emittance [mm-mrad]", )
    fma()             
  # --- emittance, total phase volume, normalized  ** warning norm emittance scaled mm-mrad by default **
  if plt_emitpv:
    # --- All Species Combined  
    for ii in sort(sp.keys()):
      s = sp[ii]        
      js = s.js
      co = s.color
      plg(hl_epspvn[0:top.jhist+1,js],top.hzbeam[0:top.jhist+1],color=co)    
    ptitles("History: Species Total Phase Volume Norm Emittance", "z [m]","Norm Emittance [mm-mrad]", )
    fma()
    # --- Target Species
    lab = "" 
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg(hl_epspvn[0:top.jhist+1,js],top.hzbeam[0:top.jhist+1],color=co)    
    ptitles("History: "+lab+" Total Phase Volume Norm Emittance","z [m]","Norm Emittance [mm-mrad]", )
    fma()             
  # --- temperature using <x'2> and corelation of x and x' via <x*x'> 
  if plt_temp and top.jhist > 0:
    xrms  = top.hxrms[0,1:top.jhist+1,0]
    emitx = top.hepsx[0,1:top.jhist+1,0]
    plg(beam.mass*top.vbeam**2*emitx**2/(jperev*16.*xrms**2),
        top.hzbeam[1:top.jhist+1] )
    yrms  = top.hyrms[0,1:top.jhist+1,0]
    emity = top.hepsy[0,1:top.jhist+1,0]
    plg(beam.mass*top.vbeam**2*emity**2/(jperev*16.*yrms**2),
        top.hzbeam[1:top.jhist+1],color="red")
    ptitles("History: Spatial Avg Temp: x [b], y [r]",
            "z [m]","Spatial Avg Temp (eV)", )
    fma()

#  -- Install diagnostics at appropriate intervals after steps
#       Add options to generate plots desired 
#  -- Install diagnostics at appropriate intervals after steps
#       Add options to generate plots desired 

# Function to call diagnostics at a timestep in step control lists 
def diag_calls():
  if top.it in diag_part_step:
    #diag_part(plt_xy=true,plt_xxp=true,plt_yyp=false,plt_xpyp=true,
              #plt_trace=false,plt_denxy=true,plt_denr=true)
    diag_part(plt_xy=true,plt_xxp=true,plt_yyp=false,plt_xpyp=false,
              plt_trace=false,plt_denxy=false,plt_denr=false)
  if top.it in diag_field_step: 
    diag_field(plt_pc=true,plt_pc_xy=true,plt_pa=true)
  if top.it in diag_hist_step:
    diag_hist(plt_ekin=true,plt_spnum=true,plt_curr_e=true,plt_curr_p=true,plt_lam_p=true,plt_lam_e=true,
              plt_lz=true,plt_pth=true,plt_pthn=true,plt_krot=true,plt_lang=true, 
              plt_cen=true,plt_envrms=true,plt_envmax=true,plt_envrmsp=true,  
              plt_emit=true,plt_emitn=true,plt_emitg=true,plt_emitng=true,plt_emitr=true,plt_emitnr=true, 
              plt_emitpv=true,emitpvn=true)

# Install diagnostic calls after simulation step
installafterstep(diag_calls)

# Step 0 diagnostics (if any) of the initial distribution loaded 
diag_calls() 

# Advance simulation specified steps 

#raise "to here"

# ---- to grated accel gap

n_step = nstep
  
#n_step = nint((neut_z1-z_launch)/wxy.ds) 
top.prwall = r_p    # consistent aperture 
step(n_step*int(nadvance))

## --- reset species weights to turn off neutralization  
#for s in sp.values():
   #s.w0 = 1. 

#loadrho()     # applies adjusted species weights  
#fieldsolve()  # make field consistent with turned off neutralization 

## --- unneutralized advance in acceleration column  
#n_step = nint((neut_z2-top.zbeam)/wxy.ds)
#top.prwall = gag_rp   # consistent aperture 
#step(n_step)

## --- reset species weights to turn on post accel gap neutralization  
#for s in sp.values():
   #s.w0 = 1.-neut_f2 

## --- neutralized advance after acceleration column to start of dipole   
#n_step = nint((z_adv-top.zbeam)/wxy.ds) + 2   # add two extra steps in case of roundoff accumulation 
#top.prwall = r_p_down    # consistent aperture 
#step(n_step)



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

# Make plot of final Brho by species 
plt_diag_bro(label = "Final Rigidity by Species") 

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

# Make sure that last plot is flushed from buffer
fma() 










from scipy.integrate import odeint
import matplotlib.pyplot as plt


# Make time array for solution

stepnum = 2000

sss = linspace(z_launch, z_adv, stepnum)

stepsize = (z_adv - z_launch)/(stepnum - 1)




# Data needed in ODE

speciesq = U_charge_states
speciesq.extend(O_charge_states)
speciesq[:] = [x*jperev for x in speciesq]


speciesI = append(U_ibeam, O_ibeam)

#speciesI = U_ibeam
#speciesI.extend(O_ibeam)

speciesCombined = U_species
speciesCombined.extend(O_species)

# state vector, for a total of n species
# ordering: same as rest of the code, U then O
# first lot of n elements are KE of the respective species
# second lot of n elements are the sigma
# third lot of elements are the dsigmadz

nstotal = U_ns + O_ns

xyz = [0]*3*nstotal

dct = {}

deltaz = stepsize/2.

callcounter = 0

termdata = []

def f(xyz, rrr):
	
	efieldz = getappliedfields(0, 0, rrr)[2][0]
	bfieldz = getappliedfields(0, 0, rrr)[5][0]
	
	dEdz = (getappliedfields(0, 0, rrr + deltaz/2)[2][0] - getappliedfields(0, 0, rrr - deltaz/2)[2][0])/deltaz
	dBdz = (getappliedfields(0, 0, rrr + deltaz/2)[5][0] - getappliedfields(0, 0, rrr - deltaz/2)[5][0])/deltaz
	
	d2Edz2 = (getappliedfields(0, 0, rrr + deltaz/2)[2][0] - getappliedfields(0, 0, rrr)[2][0]*2 + getappliedfields(0, 0, rrr - deltaz/2)[2][0])/(deltaz**2)*4
	#dVdz = 0
	
	derivs = []
	
	#for i in range(U_ns):
		#derivs.append(U_species[i].charge/jperev*efieldz) #build first lot in deriv output
	
	#for i in range(O_ns):
		#derivs.append(O_species[i].charge/jperev*efieldz) #build first lot in deriv output
	
	
	speciesbeta = []
	
	for i in range(nstotal):
		speciesbeta.append(sqrt((2*xyz[i]*jperev)/(speciesCombined[i].mass*clight**2)))
		
	
	### Modified Longitudinal Dynamics ###
	
	for j in range(nstotal):
		
		if CorrectionMode == 0:
			derivs.append(speciesq[j]*(efieldz)/jperev )
			
		if CorrectionMode == 1:
			derivs.append(speciesq[j]/jperev*(efieldz + (speciesbeta[j]*clight*hl_pth_bar[0,j]/2 - speciesq[j]*xyz[j+nstotal]**2*bfieldz/4/speciesCombined[j].mass)*dBdz) )#/speciesCombined[j].mass*dBdz) )
		
		if CorrectionMode == 2:
			derivs.append(speciesq[j]/jperev*(efieldz - xyz[j+nstotal]**2/4*d2Edz2 + (speciesbeta[j]*clight*hl_pth_bar[0,j]/2 - speciesq[j]*xyz[j+nstotal]**2*bfieldz/4/speciesCombined[j].mass)*dBdz) )#/speciesCombined[j].mass*dBdz) )
		
		#derivs.append(speciesq[j]*(efieldz)/jperev )
		#derivs.append(speciesq[j]/jperev*(efieldz + (speciesbeta[j]*clight*hl_pth_bar[0,j]/2 - speciesq[j]*xyz[j+nstotal]**2*bfieldz/4/speciesCombined[j].mass)*dBdz) )#/speciesCombined[j].mass*dBdz) )
		#derivs.append(speciesq[j]/jperev*(efieldz - xyz[j+nstotal]**2/4*d2Edz2 + (speciesbeta[j]*clight*hl_pth_bar[0,j]/2 - speciesq[j]*xyz[j+nstotal]**2*bfieldz/4/speciesCombined[j].mass)*dBdz) )#/speciesCombined[j].mass*dBdz) )
	
	for i in range(nstotal):
		dct["sigma%s" %(i)] = xyz[i+nstotal]
		dct["dsigmadz%s" %(i)] = xyz[i+2*nstotal]
		derivs.append(xyz[i+2*nstotal]) #build second lot in deriv output
	
	
	#if neut_z1 < rrr < neut_z2:
		#neut_ode = 1 - 0
	#else:
		#neut_ode = 1 - neut_f1
	
	for j in range(nstotal):
		
		scterm = 0
		
		for s in range(nstotal):
			if True:
				QQQ = (speciesq[j]*speciesI[s])/(2*pi*eps0*speciesCombined[j].mass*speciesbeta[j]**2*speciesbeta[s]*clight**3)
				scterm += QQQ*neut_ode*xyz[j+nstotal]/(xyz[j+nstotal]**2 + xyz[s+nstotal]**2)
		
		term1 = (speciesq[j]*-efieldz)/(2*xyz[j]*jperev) * xyz[j+2*nstotal]
		
		term2 = (speciesq[j]*-dEdz)/(4*xyz[j]*jperev) * xyz[j+nstotal]
		
		term3 = ((speciesq[j]*bfieldz)/(2*speciesCombined[j].mass*speciesbeta[j]*clight))**2*xyz[j+nstotal]
		
		#emitterm = (hl_epspv[0,j]**2 + (hl_pth_bar[0,j]**2) / (speciesCombined[j].mass*speciesbeta[j]*clight)**2) / xyz[j+nstotal]**3
		#emitterm = ((top.hepsr[0,0,j]/speciesbeta[j]*10**-6)**2 + (hl_pthn_bar[0,j] /speciesbeta[j])**2) / xyz[j+nstotal]**3
		#emitterm = ((0.015/speciesbeta[j]*10**-6)**2 + (hl_pthn_bar[0,j] /speciesbeta[j])**2) / xyz[j+nstotal]**3
		emitterm = ((top.hepsr[0,0,j]*CorrectionFactor)**2 + (hl_pthn_bar[0,j] /speciesbeta[j])**2) / xyz[j+nstotal]**3
		
		d2sigmadz2 = term1 + term2 - term3 + scterm + emitterm
		
		derivs.append(d2sigmadz2)
		
		global callcounter
		
		callcounter += 1
		
		termdata.append([j, rrr, term1, term2, -term3, scterm, emitterm, xyz[j]])
		
	return derivs
		
		
		
	
	

	

	

# sp["U34"].js = 1, sp["O1"].js = 16
# top.ns = 20 = total number of species
#top.jhist = 256
		
	# species ke = xyz[i]
	








# Initial States

zyx = []

for i in range(U_ns):
	zyx.append(U_ekin[i])
	
for i in range(O_ns):
	zyx.append(O_ekin[i])

#for i in range(U_ns):
	#zyx.append(sqrt(U_species[i].a0**2 + U_species[i].b0**2))
	
#for i in range(O_ns):
	#zyx.append(sqrt(O_species[i].a0**2 + O_species[i].b0**2))

#for i in range(U_ns):
	#zyx.append((U_species[i].a0*U_species[i].ap0 + U_species[i].b0*U_species[i].bp0) / zyx[i+nstotal])
	
#for i in range(O_ns):
	#zyx.append((O_species[i].a0*O_species[i].ap0 + O_species[i].b0*O_species[i].bp0) / zyx[i+U_ns+nstotal])


for i in range(U_ns):
	zyx.append(sqrt((U_species[i].a0/2)**2 + (U_species[i].b0/2)**2))
	
for i in range(O_ns):
	zyx.append(sqrt((O_species[i].a0/2)**2 + (O_species[i].b0/2)**2))

for i in range(U_ns):
	zyx.append((U_species[i].a0*U_species[i].ap0 + U_species[i].b0*U_species[i].bp0) / zyx[i+nstotal])
	
for i in range(O_ns):
	zyx.append((O_species[i].a0*O_species[i].ap0 + O_species[i].b0*O_species[i].bp0) / zyx[i+U_ns+nstotal])


psoln = odeint (f, zyx, sss, hmax = stepsize, mxstep=5000)


print (psoln[399])

warpke = []

for i in range(U_ns):
	warpke.append(mean(U_species[i].getke()))

print (warpke)

#print (array(psoln[399]-warpke)/array(warpke)










termdatawarp = []

for iii in range(0, top.jhist+1):
	
	xyz = 0
	
	zzz = top.hzbeam[iii]
	deltaz = diff(top.hzbeam)[0]/2
	
	efieldz = getappliedfields(0, 0, zzz)[2][0]
	bfieldz = getappliedfields(0, 0, zzz)[5][0]
	dEdz = (getappliedfields(0, 0, zzz + deltaz/2)[2][0] - getappliedfields(0, 0, zzz - deltaz/2)[2][0])/deltaz
	
	kelist = array([hl_ekin[iii,j] for j in range(nstotal)])
	speciesbeta = array([sqrt((2*kelist[j]*jperev)/(speciesCombined[j].mass*clight**2)) for j in range(nstotal)])
	anglelist = array([top.hxxpbar[0,iii,j]/(top.hxrms[0,iii,j])*sqrt(2) for j in range(nstotal)])
	envelopelist = array([top.hxrms[0,iii,j]*sqrt(2) for j in range(nstotal)])
	pthetalist = array([hl_pthn_bar[iii,j] for j in range(nstotal)])
	
	xyz = concatenate((kelist, envelopelist, anglelist))
	
	#if neut_z1 < zzz < neut_z2:
		#neut_ode = 1 - 0
	#else:
		#neut_ode = 1 - neut_f1
	
	for j in range(nstotal):
				
		scterm = 0
		
		for s in range(nstotal):
			if True:
				QQQ = (speciesq[j]*speciesI[s])/(2*pi*eps0*speciesCombined[j].mass*speciesbeta[j]**2*speciesbeta[s]*clight**3)
				scterm += QQQ*neut_ode*xyz[j+nstotal]/(xyz[j+nstotal]**2 + xyz[s+nstotal]**2)
		
		term1 = (speciesq[j]*-efieldz)/(2*xyz[j]*jperev) * xyz[j+2*nstotal]
		
		term2 = (speciesq[j]*-dEdz)/(4*xyz[j]*jperev) * xyz[j+nstotal]
		
		term3 = ((speciesq[j]*bfieldz)/(2*speciesCombined[j].mass*speciesbeta[j]*clight))**2*xyz[j+nstotal]
		
		#emitterm = (hl_epspv[0,j]**2 + (hl_pth_bar[0,j]**2) / (speciesCombined[j].mass*speciesbeta[j]*clight)**2) / xyz[j+nstotal]**3
		#emitterm = ((top.hepsr[0,0,j]/speciesbeta[j]*10**-6)**2 + (hl_pthn_bar[0,j] /speciesbeta[j])**2) / xyz[j+nstotal]**3
		#emitterm = ((0.015/speciesbeta[j]*10**-6)**2 + (pthetalist[j] /speciesbeta[j])**2) / xyz[j+nstotal]**3
		emitterm = ((top.hepsr[0,iii,j]*CorrectionFactor)**2 + (pthetalist[j] /speciesbeta[j])**2) / xyz[j+nstotal]**3
		
		#d2sigmadz2 = term1 + term2 - term3 + scterm + emitterm
		
		#derivs.append(d2sigmadz2)
		
		global callcounter
		
		callcounter += 1
		
		termdatawarp.append([j, zzz, term1, term2, -term3, scterm, emitterm, xyz[j]])























from scipy import interpolate

xyzxyz = [0]*3*nstotal

dct = {}

deltaz = stepsize/2.

callcounter = 0



zlist = []
		
zlist = array([top.hzbeam[kkk] for kkk in range(0, top.jhist+1)])

def fwarp(xyzxyz, rrr):
	
	efieldz = getappliedfields(0, 0, rrr)[2][0]
	bfieldz = getappliedfields(0, 0, rrr)[5][0]
	
	dEdz = (getappliedfields(0, 0, rrr + deltaz/2)[2][0] - getappliedfields(0, 0, rrr - deltaz/2)[2][0])/deltaz
	#dVdz = 0
	
	derivs = []
	
	for i in range(U_ns):
		derivs.append(U_species[i].charge/jperev*efieldz) #build first lot in deriv output
	
	for i in range(O_ns):
		derivs.append(O_species[i].charge/jperev*efieldz) #build first lot in deriv output
	
	for i in range(nstotal):
		derivs.append(xyzxyz[i+2*nstotal]) #build second lot in deriv output
	
	speciesbeta = []
	
	for i in range(nstotal):
		speciesbeta.append(sqrt((2*xyzxyz[i]*jperev)/(speciesCombined[i].mass*clight**2)))

	#if neut_z1 < rrr < neut_z2:
		#neut_ode = 1 - 0
	#else:
		#neut_ode = 1 - neut_f1

	for j in range(nstotal):
		
		scterm = 0
		
		
		
		emittancelist = []
		
		#for kkk in range(len(top.hepsr)):
			#emittancelist = emittancelist.append(top.hepsr[0,kkk,j])
		
		emittancelist = array([top.hepsr[0,kkk,j]*CorrectionFactor for kkk in range(0, top.jhist+1)])
		
		pthetaLIST = []
		
		pthetaLIST = array([hl_pthn_bar[kkk,j] for kkk in range(0, top.jhist+1)])
		
		kineticenergylist = []
		
		kineticenergylist = array([hl_ekin[kkk,j] for kkk in range(0, top.jhist+1)])
		

		emitinter = interpolate.interp1d(zlist, emittancelist, kind='slinear')
		
		pthetainter = interpolate.interp1d(zlist, pthetaLIST, kind='slinear')
		
		keinter = interpolate.interp1d(zlist, kineticenergylist, kind='slinear')
		
		for s in range(nstotal):
			if True:
				QQQ = (speciesq[j]*speciesI[s])/(2*pi*eps0*speciesCombined[j].mass*speciesbeta[j]**2*speciesbeta[s]*clight**3)
				scterm += QQQ*neut_ode*xyzxyz[j+nstotal]/(xyzxyz[j+nstotal]**2 + xyzxyz[s+nstotal]**2)
		
		term1 = (speciesq[j]*-efieldz)/(2*xyzxyz[j]*jperev) * xyzxyz[j+2*nstotal]
		
		term2 = (speciesq[j]*-dEdz)/(4*xyzxyz[j]*jperev) * xyzxyz[j+nstotal]
		
		term3 = ((speciesq[j]*bfieldz)/(2*speciesCombined[j].mass*speciesbeta[j]*clight))**2*xyzxyz[j+nstotal]
		
		#emitterm = (hl_epspv[0,j]**2 + (hl_pth_bar[0,j]**2) / (speciesCombined[j].mass*speciesbeta[j]*clight)**2) / xyzxyz[j+nstotal]**3
		#emitterm = ((top.hepsr[0,0,j]/speciesbeta[j]*10**-6)**2 + (hl_pthn_bar[0,j] /speciesbeta[j])**2) / xyzxyz[j+nstotal]**3
		#emitterm = ((0.015/speciesbeta[j]*10**-6)**2 + (hl_pthn_bar[0,j] /speciesbeta[j])**2) / xyzxyz[j+nstotal]**3
		emitterm = ((emitinter(rrr))**2 + (pthetainter(rrr) /speciesbeta[j])**2) / xyzxyz[j+nstotal]**3
		
		d2sigmadz2 = term1 + term2 - term3 + scterm + emitterm
		
		derivs.append(d2sigmadz2)
		
	return derivs
		
		


#psolnwarp = odeint (fwarp, zyx, sss, hmax = stepsize, mxstep=5000)

















#plt.figure(1)
#plt.plot(sss,psoln[:,0])

#plt.figure(2)
#plt.plot(sss,psoln[:,16])

plt.figure(3)
plt.plot(sss,psoln[:,20])

plt.figure(4)
plt.plot(sss,psoln[:,36])

termdata2 = []
termdatawarp2=[]

for i in range(0, len(termdata), 20):
	termdata2.append(termdata[i])

termdata2 = array(termdata2)

termdatawarp2 = array([termdatawarp[i] for i in  range(0, len(termdatawarp), 20)])

#plt.figure(5)
#plt.title("Term 1")
#plt.plot(termdata2[:,1],termdata2[:,2])
#plt.plot(termdatawarp2[:,1],termdatawarp2[:,2], ls='dashed')

#plt.figure(6)
#plt.title("Term 2")
#plt.plot(termdata2[:,1],termdata2[:,3])
#plt.plot(termdatawarp2[:,1],termdatawarp2[:,3], ls='dashed')

#plt.figure(7)
#plt.title("Term 3")
#plt.plot(termdata2[:,1],termdata2[:,4])
#plt.plot(termdatawarp2[:,1],termdatawarp2[:,4], ls='dashed')

#plt.figure(8)
#plt.title("Space Charge Term")
#plt.plot(termdata2[:,1],termdata2[:,5])
#plt.plot(termdatawarp2[:,1],termdatawarp2[:,5], ls='dashed')

#plt.figure(9)
#plt.title("Emittance Term")
#plt.plot(termdata2[:,1],termdata2[:,6])
#plt.plot(termdatawarp2[:,1],termdatawarp2[:,6], ls='dashed')

#plt.show()

print (psoln[:,20])



# --- Target Species, x-plane 
hpenvx(titles=false,yscale=1./(2.*mm),lhzbeam=true)
lab = ""    
for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpenvx(js=js,color=co,titles=false,yscale=1./(2.*mm),lhzbeam=true)    
ptitles("ODE (dashed), Warp (solid)","z [m]","RMS x-Envelope [mm]", )

plg(psoln[:,20]*1000/sqrt(2),sss, color="magenta", type="dash")
plg(psoln[:,21]*1000/sqrt(2),sss, color="green", type="dash")
limits(plotxmin,plotxmax, plotymin, plotymax)


combinedenvelop = [0]*stepnum

TotalWeight=0.
for j in speciesCombined:
	TotalWeight += j.sw

for i in range(stepnum):
	
	SumofR2overN = 0.
	
	for j in range(nstotal):
		#combinedenvelop[i] += psoln[i, j+nstotal]*1000/sqrt(2)*speciesI[j]/sum(speciesI)
		SumofR2overN += psoln[i, j+nstotal]**2*(speciesCombined[j].sw/TotalWeight)
	
	combinedenvelop[i] = SumofR2overN**0.5*1000/sqrt(2)

plg(combinedenvelop,sss, color="black", type="dash")

fma()	







SpeciesDictionary = {'U33':0, 'U34':1}

for kkk in range(2,10):
	abc = 'U%d' % (kkk+23)
	DummyDict = {abc:kkk}
	SpeciesDictionary.update(DummyDict)

for kkk in range(1,7):
	abc = 'U%d' % (kkk+34)
	DummyDict = {abc:kkk+9}
	SpeciesDictionary.update(DummyDict)
	
for kkk in range(1,5):
	abc = 'O%d' % kkk
	DummyDict = {abc:kkk+15}
	SpeciesDictionary.update(DummyDict)


def SpeciesEnvelope(iii):
    s = sp[iii]
    js = s.js
    co = s.color
    hpenvx(js=js,color=co,titles=false,yscale=1./(2.*mm),lhzbeam=true)

    
    
    RefNumber = SpeciesDictionary[iii]
    
    plg(psoln[:,RefNumber+20]*1000/sqrt(2),sss, color=co, type="dash")
    ptitles("ODE (dashed), Warp (solid), %s" % iii,"z [m]","RMS x-Envelope [mm]", )
    limits(plotxmin,plotxmax, plotymin, plotymax)
    
    fma()


s = sp['O1']
js = s.js
hpenvx(js=js,color='black',titles=false,yscale=1./(2.*mm),lhzbeam=true)
RefNumber = SpeciesDictionary['O1']
plg(psoln[:,RefNumber+20]*1000/sqrt(2),sss, color='black', type="dash")

s = sp['O4']
js = s.js
hpenvx(js=js,color='red',titles=false,yscale=1./(2.*mm),lhzbeam=true)
RefNumber = SpeciesDictionary['O4']
plg(psoln[:,RefNumber+20]*1000/sqrt(2),sss, color='red', type="dash")

ptitles("ODE (dashed), Warp (solid), O1 black, O4 red","z [m]","RMS x-Envelope [mm]", )
limits(plotxmin,plotxmax, plotymin, plotymax)

fma()


for jjj in sp:
	SpeciesEnvelope(jjj)










for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      A  = s.mass/amu
      #hpekin(js=js,color=co,titles=false,yscale=1./A,lhzbeam=true)
      plg(hl_ekin[0:top.jhist+1,js]/(A*kV),top.hzbeam[0:top.jhist+1],color=co)        
ptitles("Kinetic Energy, ODE (dashed), Warp (solid)","z [m]","KeV/u", )

plg(psoln[:,0]/1000/238.02891,sss, color="magenta", type="dash")
plg(psoln[:,1]/1000/238.02891,sss, color="green", type="dash")

fma()



# --- Target Species, x-plane 
lab = ""    
for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg(top.hepsr[0, 0:top.jhist+1,js]*CorrectionFactor,top.hzbeam[0:top.jhist+1],color=co)
      plg(sqrt(4*top.hepsx[0, 0:top.jhist+1,js]**2 - hl_lz_bar2[0:top.jhist+1,js]**2),top.hzbeam[0:top.jhist+1],color=co, type="dash")    
ptitles("History: "+lab+"RMS x-Envelope","z [m]","RMS Width [mm]", )

fma()






def plotallODEterms(jj):
	
	termdata4 = []
	termdatawarp4=[]

	
	termdata4 = array([termdata[i] for i in range(jj, len(termdata), nstotal)])

	termdatawarp4 = array([termdatawarp[i] for i in  range(jj, len(termdatawarp), nstotal)])
	
	plg(termdata4[:,2], termdata4[:,1], color="black")
	plg(termdata4[:,3], termdata4[:,1], color="green")
	plg(termdata4[:,4], termdata4[:,1], color="yellow")
	plg(termdata4[:,5], termdata4[:,1], color="red")
	plg(termdata4[:,6], termdata4[:,1], color="blue")
	
	ptitles("ODE Terms, black:E green:dEdz yellow:B red:SC blue:Emit","z [m]","sigma'' [rad/m]",)
	
	fma()


plotallODEterms(0)






def plotallWARPterms(jj):
	
	termdata4 = []
	termdatawarp4=[]

	
	termdata4 = array([termdata[i] for i in range(jj, len(termdata), nstotal)])

	termdatawarp4 = array([termdatawarp[i] for i in  range(jj, len(termdatawarp), nstotal)])
	
	plg(termdatawarp4[:,2], termdatawarp4[:,1], color="black")
	plg(termdatawarp4[:,3], termdatawarp4[:,1], color="green")
	plg(termdatawarp4[:,4], termdatawarp4[:,1], color="yellow")
	plg(termdatawarp4[:,5], termdatawarp4[:,1], color="red")
	plg(termdatawarp4[:,6], termdatawarp4[:,1], color="blue")
	
	ptitles("Terms with WARP input, black:E green:dEdz yellow:B red:SC blue:Emit","z [m]","sigma'' [rad/m]",)
	
	fma()


plotallWARPterms(0)







def termsodevswarp(jj):
	
	termdata3 = []
	termdatawarp3=[]

	for i in range(jj, len(termdata), nstotal):
		termdata3.append(termdata[i])

	termdata3 = array(termdata3)

	termdatawarp3 = array([termdatawarp[i] for i in  range(jj, len(termdatawarp), nstotal)])


	plg(termdata3[:,2],termdata3[:,1])
	plg(termdatawarp3[:,2],termdatawarp3[:,1], type="dash")
	ptitles("Term 1 (E-field), ODEinput (solid), WARPinput (dashed)","z [m]","sigma'' [rad/m]",)
	fma()
	
	plg(termdata3[:,3],termdata3[:,1])
	plg(termdatawarp3[:,3],termdatawarp3[:,1], type="dash")
	ptitles("Term 2 (dEdz), ODEinput (solid), WARPinput (dashed)","z [m]","sigma'' [rad/m]",)
	fma()
	
	plg(termdata3[:,4],termdata3[:,1])
	plg(termdatawarp3[:,4],termdatawarp3[:,1], type="dash")
	ptitles("Term 3 (B-field), ODEinput (solid), WARPinput (dashed)","z [m]","sigma'' [rad/m]",)
	fma()
	
	plg(termdata3[:,5],termdata3[:,1])
	plg(termdatawarp3[:,5],termdatawarp3[:,1], type="dash")
	ptitles("Space Charge Term, ODEinput (solid), WARPinput (dashed)","z [m]","sigma'' [rad/m]",)
	fma()
	
	plg(termdata3[:,6],termdata3[:,1])
	plg(termdatawarp3[:,6],termdatawarp3[:,1], type="dash")
	ptitles("Emittance Term, ODEinput (solid), WARPinput (dashed)","z [m]","sigma'' [rad/m]",)
	fma()



termsodevswarp(0)

#termsodevswarp(16)










def termsdifference(jj):
	
	termdata4 = []
	termdatawarp4=[]

	
	termdata4 = array([termdata[i] for i in range(jj, len(termdata), nstotal)])

	termdatawarp4 = array([termdatawarp[i] for i in  range(jj, len(termdatawarp), nstotal)])
	
	warpterm1inter = interpolate.interp1d(termdatawarp4[:,1], termdatawarp4[:,2], bounds_error = False, fill_value = 0., kind='linear')
	warpterm2inter = interpolate.interp1d(termdatawarp4[:,1], termdatawarp4[:,3], bounds_error = False, fill_value = 0., kind='linear')
	warpterm3inter = interpolate.interp1d(termdatawarp4[:,1], termdatawarp4[:,4], bounds_error = False, fill_value = 0., kind='linear')
	warpscterminter = interpolate.interp1d(termdatawarp4[:,1], termdatawarp4[:,5], bounds_error = False, fill_value = 0., kind='linear')
	warpemitterminter = interpolate.interp1d(termdatawarp4[:,1], termdatawarp4[:,6], bounds_error = False, fill_value = 0., kind='linear')
	
	warpterm1 = array([ warpterm1inter(kkk) for kkk in termdata4[:,1] ])
	warpterm2 = array([ warpterm2inter(kkk) for kkk in termdata4[:,1] ])
	warpterm3 = array([ warpterm3inter(kkk) for kkk in termdata4[:,1] ])
	warpscterm = array([ warpscterminter(kkk) for kkk in termdata4[:,1] ])
	warpemitterm = array([ warpemitterminter(kkk) for kkk in termdata4[:,1] ])
	
	if len(termdata4[:,2]) != len(warpterm1):
		print "termdata4[:,2] not same length as warpterm1"
	
	term1diff = termdata4[:,2] - warpterm1
	term2diff = termdata4[:,3] - warpterm2
	term3diff = termdata4[:,4] - warpterm3
	sctermdiff = termdata4[:,5] - warpscterm
	emittermdiff = termdata4[:,6] - warpemitterm
	
	print "calling termdiffernce"
	
	plg(term1diff, termdata4[:,1], color="black")
	plg(term2diff, termdata4[:,1], color="green")
	plg(term3diff, termdata4[:,1], color="yellow")
	plg(sctermdiff, termdata4[:,1], color="red")
	plg(emittermdiff, termdata4[:,1], color="blue")
	
	ptitles("ODEinput minus WARPinput for different terms","z [m]","sigma'' [rad/m]",)
	
	fma()




termsdifference(0)

#termsdifference(16)





def kedifference(jj):
	
	termdata4 = []
	termdatawarp4=[]

	
	termdata4 = array([termdata[i] for i in range(jj, len(termdata), nstotal)])

	termdatawarp4 = array([termdatawarp[i] for i in  range(jj, len(termdatawarp), nstotal)])
	
	warpkeinter = interpolate.interp1d(termdatawarp4[:,1], termdatawarp4[:,7], bounds_error = False, fill_value = 0., kind='linear')
	
	warpke = array([ warpkeinter(kkk) for kkk in termdata4[:,1] ])
	
	if len(termdata4[:,2]) != len(warpke):
		print "termdata4[:,2] not same length as warpterm1"
	
	kediff = termdata4[:,7] - warpke
	
	kediffcut = array(kediff[4:-4])
	warpkecut = array(warpke[4:-4])
	
	kepercentagediff = kediffcut/warpkecut*100
	
	print "calling kediffernce"
	
	if jj == 0:
		plg(kepercentagediff, termdata4[:,1][4:-4], color="magenta")
	elif jj == 1:
		plg(kepercentagediff, termdata4[:,1][4:-4], color="green")
	else:
		plg(kepercentagediff, termdata4[:,1][4:-4], color="black")
	

kedifference(0)
kedifference(1)

ptitles("KE Difference between ODE and WARP","z [m]","% difference",)
	
fma()




# --- Target Species, x-plane 
#hpenvx(titles=false,yscale=1./(2.*mm),lhzbeam=true)
#lab = ""    
#for ii in sp_target:
      #s = sp[ii]
      #js = s.js
      #co = s.color
      #lab+= ii + "("+co+"), "
      #hpenvx(js=js,color=co,titles=false,yscale=1./(2.*mm),lhzbeam=true)    
#ptitles("ODE with WARP input (dashed), WARP Simulations (solid)","z [m]","RMS x-Envelope [mm]", )

#plg(psolnwarp[:,20]*1000/sqrt(2),sss, color="magenta", type="dash")
#plg(psolnwarp[:,21]*1000/sqrt(2),sss, color="green", type="dash")


#combinedenvelop = [0]*stepnum

#TotalWeight=0.
#for j in speciesCombined:
	#TotalWeight += j.sw

#for i in range(stepnum):
	
	#SumofR2overN = 0.
	
	#for j in range(nstotal):
		##combinedenvelop[i] += psoln[i, j+nstotal]*1000/sqrt(2)*speciesI[j]/sum(speciesI)
		#SumofR2overN += psolnwarp[i, j+nstotal]**2*(speciesCombined[j].sw/TotalWeight)
	
	#combinedenvelop[i] = SumofR2overN**0.5*1000/sqrt(2)

#plg(combinedenvelop,sss, color="black", type="dash")	

#fma()




# --- Target Species, x-plane 

#ptitles("ODE (solid), ODE with WARP input (dashed)","z [m]","RMS x-Envelope [mm]", )

#plg(psoln[:,20]*1000/sqrt(2),sss, color="magenta")
#plg(psoln[:,21]*1000/sqrt(2),sss, color="green")
#plg(psolnwarp[:,20]*1000/sqrt(2),sss, color="magenta", type="dash")
#plg(psolnwarp[:,21]*1000/sqrt(2),sss, color="green", type="dash")

#combinedenvelop = [0]*stepnum

#TotalWeight=0.
#for j in speciesCombined:
	#TotalWeight += j.sw

#for i in range(stepnum):
	
	#SumofR2overN = 0.
	
	#for j in range(nstotal):
		##combinedenvelop[i] += psoln[i, j+nstotal]*1000/sqrt(2)*speciesI[j]/sum(speciesI)
		#SumofR2overN += psolnwarp[i, j+nstotal]**2*(speciesCombined[j].sw/TotalWeight)
	
	#combinedenvelop[i] = SumofR2overN**0.5*1000/sqrt(2)

#plg(combinedenvelop,sss, color="black", type="dash")	



#combinedenvelop = [0]*stepnum

#TotalWeight=0.
#for j in speciesCombined:
	#TotalWeight += j.sw

#for i in range(stepnum):
	
	#SumofR2overN = 0.
	
	#for j in range(nstotal):
		##combinedenvelop[i] += psoln[i, j+nstotal]*1000/sqrt(2)*speciesI[j]/sum(speciesI)
		#SumofR2overN += psoln[i, j+nstotal]**2*(speciesCombined[j].sw/TotalWeight)
	
	#combinedenvelop[i] = SumofR2overN**0.5*1000/sqrt(2)

#plg(combinedenvelop,sss, color="black")	

#fma()

ptitles("Solenoid field: %s T, Lattice Length: %s m, r0: %s m" %(SolenoidStrength, LatticeLength, InitialRadius))

fma()

ptitles("Fields: %s, f_neut: %s, Correction: %s, CentralDifference" %(abc, NeutralizationFactor, CorrectionMode))

fma()

electron_tracker.pzx(i=0,color=red)
electron_tracker.pzx(i=1,color=black)

fma()

electron_tracker.pzy(i=0,color=red)
electron_tracker.pzy(i=1,color=black)

fma()
