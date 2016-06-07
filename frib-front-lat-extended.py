# Lattice Description of FRIB Front End 
# Notes:
#  * Use coordinates in FRIB lattice file. These are set so the start of the RFQ is at coordinate 200 meters.  
#    This is large enough where source coordinates are all positive.  
#    - Use of this may create some confusion (simulations start at some odd z-value), but hopefully make it 
#      simpler to understand coordinates of elements consistently with the construction specification.  


# Lattice periodicity  
top.zlatperi  = largepos # periodicity length [m]  
top.zlatstrt  = 0.       # z of lattice start; added to element z's [m] 
                         #   (can use to change lattice phase) 


# Energy and Bias of Source High Voltage Platform 

ekin_per_u = 12.*keV                             # target kinetic energy/u for LEBT post ES gap  

StandBias  = A_ref*ekin_per_u/Q_ref - SourceBias # Conistent Bias of Injector Column
Bias       = StandBias + SourceBias              # Total bias to achieve ekin_per_u 


# Venus ECR Source Fields 
# Comment: Must have same z-grids for linear and nonlinear forms. 

# --- element specification 

ecr_shift  = 11.*cm                 # shift of ecr from lattice file spec to make room for s4p1 
ecr_z_extr = 66.650938 - ecr_shift  # z-location of beam extraction aperture in simulation coordinates     
ecr_sc     = 1.0                    # scale factor to muliply field data by 
ecr_typ    = "lin"                  # type: "lin" = linear optics fields or "nl" = nonlinear r-z field  

# --- linear element data  
#fi = PRpickle.PR("lat_ecr_venus.lin.20140602.pkl")  old data file which does not include fringe field extension
fi = PRpickle.PR("ecr_venus/lat_ecr_venus.lin.20160218.pkl")
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

# --- define venus ecr fields  
if ecr_typ == "lin":
  ecr = addnewmmlt(zs=ecr_zmmin,ze=ecr_zmmax,id=ecr_lin_id,sc=ecr_sc) 
elif ecr_typ == "nl":
  #addnewbgrd(xs=0.,zs=s41_zc-s4_zlen/2.,ze=s41_zc+s4_zlen/2.,id=s4_nl_id,func=s41_scale)
  raise Exception("No ECR Venus Nonlinear Applied Fields Defined") 
  ecr = None
else:
  print("Warning: No ECR Applied Fields Defined") 
  ecr = None


# S4 solenoids 
# Comment: linear and nonlinear variants must have same z-grid. 

# --- element specification 

s4p1_zc  = 66.956900   # S4 1: z-center  
s4p1_str = 0.6 # 0.754 # S4 1: peak on-axis B_z field strength [Tesla]
s4p1_typ = "nl"        # S4 1: type: "lin" = linear optics fields or "nl" = nonlinear r-z field  

s4p2_zc  = 68.306900   # S4 2: z-center 
s4p2_str = 0.5 # 0.617 # s4 2: peak on-axis B_z field strength [Tesla]
s4p2_typ = "nl"        # S4 1: type: "lin" = linear optics fields or "nl" = nonlinear r-z field  

# --- linear element data  
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

# --- nonlinear element field data 
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

# --- nonlinear element vector potential data 
#fi = PRpickle.PR('lat_s4.at.20140603.pkl') 
fi = PRpickle.PR('lat_s4.at.20141031.pkl') 
#
if fi.s4_nz != s4_nz: raise Exception("S4: Nonlin Vector potential model nz not equal to nonlinear/linear model nz")
if fi.s4_nr != s4_nr: raise Exception("S4: Nonlin Vector potential model nr not equal to nonlinear model nr")
s4_at_m  = fi.s4_at_m
fi.close() 

# --- Axisymmetric b-field arrays must be 3d shape (nr+1,arb,nz+1) to load into Warp  
s4_br_m = fzeros((s4_nr+1,1,s4_nz+1))  
s4_br_m[:,0,:] = s4_br_m_in
s4_bz_m = fzeros((s4_nr+1,1,s4_nz+1))
s4_bz_m[:,0,:] = s4_bz_m_in

s4_nl_id = addnewbgrddataset(dx=s4_dr,dy=1.,zlength=s4_zlen,bx=s4_br_m,bz=s4_bz_m,rz = true)  # pass arb dy to avoid error trap  

s4_aspect = s4_r_coil_i/s4_len_coil 

# --- define solenoid s4 1 
if s4p1_typ == "lin":
  s4p1 = addnewmmlt(zs=s4p1_zc-s4_zlen/2.,ze=s4p1_zc+s4_zlen/2.,id=s4_lin_id,sc=s4p1_str) 
elif s4p1_typ == "nl":
  s4p1 = addnewbgrd(xs=0.,zs=s4p1_zc-s4_zlen/2.,ze=s4p1_zc+s4_zlen/2.,id=s4_nl_id,sc=s4p1_str)
else:
  print("Warning: No S4 1st Solenoid Applied Fields Defined") 
  s4p1 = None

# --- define solenoid s4 2 
if s4p2_typ == "lin":
  s4p2 = addnewmmlt(zs=s4p2_zc-s4_zlen/2.,ze=s4p2_zc+s4_zlen/2.,id=s4_lin_id,sc=s4p2_str) 
elif s4p2_typ == "nl":
  s4p2 = addnewbgrd(xs=0.,zs=s4p2_zc-s4_zlen/2.,ze=s4p2_zc+s4_zlen/2.,id=s4_nl_id,sc=s4p2_str)
else:
  print("Warning: No S4 2nd Solenoid Applied Fields Defined") 
  s4p2 = None

# --- Define vector potential function for both linear and nonlinear solenoid magnetic fields  
def getatheta(r):
  # --- gather vector potential 
  n = len(r) 
  at = zeros(n)
  at_scratch = zeros(n) 
  z  = top.zbeam*ones(n)
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
  if top.zbeam >= s4p1_zc-s4_zlen/2. and top.zbeam <= s4p1_zc+s4_zlen/2.:
    # --- contribution from 1st s4 
    if s4p1_typ == "lin": 
      getgrid1d(n,z,at_scratch,s4_nz,s4p1_str*s4_bz0_m,s4p1_zc-s4_zlen/2.,s4p1_zc+s4_zlen/2.)
      at_scratch = at_scratch*r/2.
    elif s4p1_typ == "nl":
      getgrid2d(n,r,z,at_scratch,s4_nr,s4_nz,s4p1_str*s4_at_m,s4_r_m.min(),s4_r_m.max(), 
                s4p1_zc-s4_zlen/2.,s4p1_zc+s4_zlen/2.)
    else:
      raise Exception("Vector Potential: S4.1 not defined")
    at += at_scratch
  if top.zbeam >= s4p2_zc-s4_zlen/2. and top.zbeam <= s4p2_zc+s4_zlen/2.:
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


# Grated Acceleration Gap
#  Note: for ideal zero-length gap:  top.lacclzl=true for zero length gap.  Accel given given by acclez*(accelze-acclzs) 
#   see dave grote email on caution on setting top.acclsw for gaps.   
#   Comment: Linear and nonlinear forms must have same axial grid.  

# --- element specification 
gag_zc  = 67.811564  # Grated Accel Gap: z-center  
gag_typ = "nl"       # Grated Accel Gap: type: "ideal" = Short gap kick, "lin" = linear r-z field imported, "nl" = nonlinear r-z field imported   

# --- linear element data  
# fi = PRpickle.PR("lat_gag.lin.20140624.pkl")  # Original Warp model with simplified geometry  
fi = PRpickle.PR("lat_gag.lin.20141029.pkl")    # Poisson model with high detail 
gag_dz = fi.gag_dz0 
gag_nz = fi.gag_nz0  
gag_z_m     = fi.gag_z0_m  
gag_ez0_m   = fi.gag_ez0_m
gag_ez0p_m  = fi.gag_ez0p_m
fi.close() 

gag_zlen = gag_z_m.max() - gag_z_m.min() 

gag_lin_id = addnewemltdataset(zlen=gag_zlen,es=gag_ez0_m,esp=gag_ez0p_m,nn=0,vv=0)

# --- nonlinear element data
#fi = PRpickle.PR('lat_gag.rz.20140624.pkl')  # Original Warp model with simplified geometry  
fi = PRpickle.PR('lat_gag.rz.20141029.pkl')   # Poisson model with high detail 
if fi.gag_nz != gag_nz: raise Exception("GAG: Nonlinear and linear field model nz not equal") 
gag_nr = fi.gag_nr
gag_dr = fi.gag_dr
gag_r_m = fi.gag_r_m
gag_z_m_cen = fi.gag_z_m_cen
gag_phi_m    = fi.gag_phi_m 
gag_er_m_in  = fi.gag_er_m
gag_ez_m_in  = fi.gag_ez_m
fi.close() 

gag_zlen = gag_z_m.max()-gag_z_m.min()          # axial length nonlin/lin structure on mesh 

gag_zs   = gag_zc - (gag_z_m_cen-gag_z_m.min()) # z_start of nonlin/lin mesh structure  
gag_ze   = gag_zc + (gag_z_m.max()-gag_z_m_cen) # z_end   of nonlin/lin mesh structure 

# Geometry parameters ?? Read these in grated gap file to be safe .. ** must be consistent with input data ** ?? 
gag_rp = 7.3*cm                  # pipe radius of inner extent of rings in grated gap 
gag_col_zs = gag_zc - 11.989*cm  # z-start (at end   of biased upstream pipe)     of grated gap mechanical structure 
gag_col_ze = gag_zc + 15.611*cm  # z-end   (at start of grounded downstream pipe) of grated gap mechanical structure  

gag_er_m = fzeros((gag_nr+1,1,gag_nz+1)) # Axisymmetric e-field arrays must be 3d shape (nr+1,arb,nz+1) to load into Warp  
gag_er_m[:,0,:] = gag_er_m_in
gag_ez_m = fzeros((gag_nr+1,1,gag_nz+1))
gag_ez_m[:,0,:] = gag_ez_m_in

gag_nl_id = addnewegrddataset(dx=gag_dr,dy=1.,zlength=gag_zlen,ex=gag_er_m,ez =gag_ez_m,rz = true) 

# --- define grated acceleration gap  
if gag_typ == "ideal":
  print("Warning: No Ideal Acceleration Gap model yet implemented")  
  gag = None 
elif gag_typ == "lin":
  gag = addnewemlt(zs=gag_zs,ze=gag_ze,id=gag_lin_id,sc=StandBias) 
elif gag_typ == "nl":
  gag = addnewegrd(xs=0.,zs=gag_zs,ze=gag_ze,id=gag_nl_id,sc=StandBias)
else:
  print("Warning: No Grated Acceleration Gap Applied Fields Defined") 
  gag = None 


# D5 Bending Dipole 

# --- element specification 

d5p1_zc  = 69.581900   # D5 1: z-center  
d5p1_str = 1.0         # D5 1: Input field scale factor
d5p1_typ = "ideal"        # D5 1: type: "ideal" = uniform By, "lin" = linear optics fields, "3d" = 3d field  

# --- nonlinear element data 
fi = PRpickle.PR('lat_d5.3d.20140527.pkl') 
d5_3d_nx = fi.d5_nx
d5_3d_ny = fi.d5_ny
d5_3d_nz = fi.d5_nz
d5_3d_dx = fi.d5_dx
d5_3d_dy = fi.d5_dy
d5_3d_dz = fi.d5_dz
d5_3d_x_m = fi.d5_x_m
d5_3d_y_m = fi.d5_y_m
d5_3d_z_m = fi.d5_z_m
d5_3d_z_m_cen = fi.d5_z_m_cen
d5_3d_bx_m = fi.d5_bx_m
d5_3d_by_m = fi.d5_by_m
d5_3d_bz_m = fi.d5_bz_m
fi.close() 
d5_3d_zlen = d5_3d_z_m.max() - d5_3d_z_m.min()

d5_3d_id = addnewbgrddataset(dx=d5_3d_dx,dy=d5_3d_dy,zlength=d5_3d_zlen,bx=d5_3d_bx_m,by=d5_3d_by_m,bz =d5_3d_bz_m) 

# Starting and ending position of first ideal D5 dipole.
# Also used to define the lattice bend

d5p1_zs = 69.2 # may need to be revised upon obtaining lattice design data
d5p1_ze = d5p1_zc + (d5p1_zc - d5p1_zs)

# --- define dipole d5 

# --- ideal (uniform) field 
if d5p1_typ == "ideal": 
  bending_R = (d5p1_ze - d5p1_zs)/(pi/2.)
  bending_B = sqrt( A_ref*ekin_per_u*jperev*2.*A_ref*amu)/(Q_ref*jperev)/bending_R
  d5p1 = addnewdipo(zs = d5p1_zs, ze = d5p1_ze, by = bending_B)
# --- linear optic approximation field 
elif d5p1_typ == "lin":
  print("Warning: No D5 1st Dipole Linear Applied Fields Defined")
  d5p1 = None
# --- 3D field from magnet design code
elif d5p1_typ == "nl":
  d5p1 = addnewbgrd(dx=d5_3d_dx,dy=d5_3d_dy,xs=d5_3d_x_m.min(),ys=d5_3d_y_m.min(),
    zs=d5p1_zc-d5_3d_zlen/2.,ze=d5p1_zc+d5_3d_zlen/2.,id=d5_3d_id,sc=d5p1_str)
else:
  print("Warning: No D5 1st Dipole Applied Fields Defined") 
  d5p1 = None


#
# Lattice bends for D5 Bending Dipole 
#
#  Comments:
#   * Use ideal bends on lattice whether dipole field is ideal (uniform) or not.   

d5p1_bend = True  # True or False: Add ideal bend to lattice 


if d5p1_bend:
  top.diposet = False     # turn off By that automatically generated with addnewbend()
  addnewbend(zs = d5p1_zs, ze = d5p1_ze, rc = (d5p1_ze - d5p1_zs)/(pi/2.))





# Neutralization specifications 
#
#    Define neutralization fraction function 
#       rho_neut_f(z,s) 
#    that returns the neutralization fraction based on axial location z within the lattice for the named 
#    species "s". Values of rho_neut_f() are set by dictionary named arrays specified 
#    by species with:
#      
#       neut_z    = []   Break points in lattice with neutralization fraction [m]
#       neut_frac = []   Neutralization fraction to the right of corresponding break point: 
#                         0 = no neut, 1 = full neut 
#
#    * neut_z values must be in numerical order. 
#    * neut_frac should be the same dimension as neut_z 
#    * Comment: Should generalize this structure to allow for different neutralization factors 
#      by species.  
#


sp_neut_z = \
array(
[ecr_z_extr - 10.*cm, # z to the left of ECR extraction point ... beam should be launched to right  
 gag_zc - 20.90*cm,   # z of neut stop before grated gap, set where 1% of gap E_z field reached 
 gag_zc + 22.28*cm    # z of neut stop after  grated gap, set where 1% of gap E_z field reached
]    )

sp_neut_frac = \
array(
[0.75, 
 0., 
 0.75
]    )

neut_z    = {key: sp_neut_z    for key in sp.keys()}
neut_frac = {key: sp_neut_frac for key in sp.keys()}

def rho_neut_f(z,s):
  # --- Find index giving location in neutralization fraction array 
  index = sum(where(z/neut_z[s] >= 1., 1, 0)) - 1  
  if index < 0: index = 0 
  # --- Return neutralization fraction with no error checking to allow flexability
  f = neut_frac[s][index] 
  return(f) 


# Aperture specfications 
#   
#   r_ap   = array of aperture radii 
#   v_ap   = array of aperture bias voltages 
#   z_ap_l = array of aperture extent lower ranges 
#   z_ap_u = array of aperture extent upper ranges 
#  
#   aperture_r(z):  function returning aperture at value of z 
#
#   * Arrays should be same dimension. 
#   * Aperture is loaded into the xy simulation for fieldsolves.  
#   * Aperture and bias does not matter in xy simulations, but try to use the right values for 
#     consistent reference potential. 
#   

r_p_up   = 8.00*cm  # aperture   upstream of grated gap [m]
r_p_down = 7.62*cm  # aperture downstream of grated gap [m]

r_ap   = array([r_p_up,               gag_rp,       r_p_down  ])
v_ap   = array([SourceBias+StandBias, StandBias/2., 0.        ]) 
z_ap_l = array([ecr_z_extr,           gag_col_zs,   gag_col_ze])
z_ap_u = array([gag_col_zs,           gag_col_ze,   d5p1_zc   ])

r_p = max(r_ap)   # Max aperture in simulations 

def aperture_r(z):
  index = sum(where(z/z_ap_l >= 1.,1,0))-1 
  if index < 0: index = 0
  # 
  return(r_ap[index])  

aperture = [] 
for i in range(len(r_ap)):
  rp = r_ap[i] 
  v  = v_ap[i] 
  zl = z_ap_l[i] 
  zu = z_ap_u[i]
  #
  aperture.append( ZCylinderOut(radius=rp,zlower=zl,zupper=zu,condid="next") )

# --- Add a circular aperture particle scraper at pipe radius aperture. 
#       Could also use ParticleScraper(conductors=aperture) but 
#       setting prwall is faster for a simple cylinder. 
top.prwall = r_p    # reset this later consistent with actual aperture in range simulated in advances 

