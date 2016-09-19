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

ref_gamma_post_gap = 1. + ekin_per_u*jperev/(amu*clight**2)
ref_vel_post_gap = clight*sqrt(1. - 1./ref_gamma_post_gap**2)
ref_brho_post_gap = ref_gamma_post_gap*ref_vel_post_gap*A_ref*amu/(Q_ref*jperev)

# Venus ECR Source Fields 
# Comment: Must have same z-grids for linear and nonlinear forms. 

# --- element specification 

ecr_shift  = 11.*cm                 # shift of ecr from lattice file spec to make room for s4p1 
ecr_z_extr = 66.650938 - ecr_shift  # z-location of beam extraction aperture in simulation coordinates     
ecr_sc     = 1.0                    # scale factor to muliply field data by 
ecr_typ    = "lin"                  # type: "lin" = linear optics fields or "nl" = nonlinear r-z field  

# --- linear element data  
#fi = PRpickle.PR("lat_ecr_venus/lat_ecr_venus.lin.20140602.pkl")  old data file which does not include fringe field extension
fi = PRpickle.PR("lat_ecr_venus/lat_ecr_venus.lin.20160218.pkl")
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
#fi = PRpickle.PR("lat_s4/lat_s4.lin.20140603.pkl")
fi = PRpickle.PR("lat_s4/lat_s4.lin.20141031.pkl")
s4_dz  = fi.s4_dz 
s4_nz  = fi.s4_nz  
s4_z_m = fi.s4_z_m 
s4_bz0_m   = fi.s4_bz0_m
s4_bz0p_m  = fi.s4_bz0p_m
fi.close() 

s4_zlen = s4_z_m.max() - s4_z_m.min() 
s4_lin_id = addnewmmltdataset(zlen=s4_zlen,ms=s4_bz0_m,msp=s4_bz0p_m,nn=0,vv=0)

# --- nonlinear element field data 
#fi = PRpickle.PR('lat_s4/lat_s4.rz.20140603.pkl') 
fi = PRpickle.PR('lat_s4/lat_s4.rz.20141031.pkl') 
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
#fi = PRpickle.PR('lat_s4/lat_s4.at.20140603.pkl') 
fi = PRpickle.PR('lat_s4/lat_s4.at.20141031.pkl') 
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
# fi = PRpickle.PR("lat_gag/lat_gag.lin.20140624.pkl")  # Original Warp model with simplified geometry  
fi = PRpickle.PR("lat_gag/lat_gag.lin.20141029.pkl")    # Poisson model with high detail 
gag_dz = fi.gag_dz0 
gag_nz = fi.gag_nz0  
gag_z_m     = fi.gag_z0_m  
gag_ez0_m   = fi.gag_ez0_m
gag_ez0p_m  = fi.gag_ez0p_m
fi.close() 

gag_zlen = gag_z_m.max() - gag_z_m.min() 

gag_lin_id = addnewemltdataset(zlen=gag_zlen,es=gag_ez0_m,esp=gag_ez0p_m,nn=0,vv=0)

# --- nonlinear element data
#fi = PRpickle.PR('lat_gag/lat_gag.rz.20140624.pkl')  # Original Warp model with simplified geometry  
fi = PRpickle.PR('lat_gag/lat_gag.rz.20141029.pkl')   # Poisson model with high detail 
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










IdealCSS = True # if true, all positions in CSS are defined relative to 1st D5 in accordance with DIMAD design
                # idealCSS involves repositioning of the elements in additional to using ideal elements

# D5 Bending Dipole 

# --- element specification 

d5p1_zc  = 69.587759   # D5 1: z-center  
d5p1_str = 1.0         # D5 1: Input field scale factor
d5p1_typ = "nl"        # D5 1: type: "ideal" = uniform By, "lin" = linear optics fields, "3d" = 3d field  
d5p1_ideal_len = 1.0

d5p2_zc  = 73.781896   # D5 2: z-center  
if IdealCSS == True:
	d5p2_zc = d5p1_zc + 4.2
d5p2_str = 1.0         # D5 2: Input field scale factor
d5p2_typ = "nl"        # D5 2: type: "ideal" = uniform By, "lin" = linear optics fields, "3d" = 3d field
d5p2_ideal_len = 1.0

# --- nonlinear element data 
if False: # OLD VERSION
  fi = PRpickle.PR('lat_d5/lat_d5.3d.20140527.pkl') 
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

  d5_3d_id = addnewbgrddataset(dx=d5_3d_dx,dy=d5_3d_dy,zlength=d5_3d_zlen,
                               bx=d5_3d_bx_m,by=d5_3d_by_m,bz =d5_3d_bz_m)

if d5p1_typ == "nl":
  #Bend magnet grid data info
  d5_3d_cond_in = 43.0*cm                       # Conductor inner radius
  d5_3d_cond_out = 84.0*cm                      # Conductor outer radius
  #d5_3d_x_ap = (d5_3d_cond_out - d5__3dcond_in)/2.0  # x aperture size
  #d5_3d_y_ap = 5.0*cm                          # y aperture size
  d5_3d_rc = (d5_3d_cond_out + d5_3d_cond_in)/2.0     # Magnet center radius
  d5_3d_core_l = d5_3d_rc*0.5*pi                # Magnet core length
  d5_3d_s = 69.2                                # Magnet core start point
  d5_3d_e = d5_3d_s + d5_3d_core_l              # Magnet core end point

  d5_3d_frng_l = 53.0*cm                        # Fringe field length
  d5_3d_xw = 42.0*cm                            # B grid data X width
  d5_3d_yw = 10.0*cm                            # B grid data Y width
  d5_3d_zlen = d5_3d_frng_l*2.0 + d5_3d_core_l  # B grid data Z length
  d5_3d_str = d5_3d_s - d5_3d_frng_l            # B grid data start point
  d5_3d_end = d5_3d_str+d5_3d_zlen              # B grid data end point

  d5_3d_scl = 0.6015157305571277 * 1.00132179   # B grid scaling factor (center_line_By *correction_factor)

  d5_3d_nx = 105   # Number of X grid
  d5_3d_ny = 25    # Number of Y grid
  d5_3d_nz = 500   # Number of Z grid

  d5_3d_dx = d5_3d_xw/float64(d5_3d_nx)     # X grid size
  d5_3d_dy = d5_3d_yw/float64(d5_3d_ny)     # Y grid size
  d5_3d_dz = d5_3d_zlen/float64(d5_3d_nz)   # Z grid size

  bgdata = getdatafromtextfile("lat_d5/bend_trans.table",dims=[3,None],)
  if len(bgdata[0]) != (d5_3d_nx+1)*(d5_3d_ny+1)*(d5_3d_nz+1): raise Exception("bend grid data is invalid.")

  d5_3d_bx = resize(bgdata[0],(d5_3d_nx+1,d5_3d_ny+1,d5_3d_nz+1))*gauss
  d5_3d_by = resize(bgdata[1],(d5_3d_nx+1,d5_3d_ny+1,d5_3d_nz+1))*gauss
  d5_3d_bz = resize(bgdata[2],(d5_3d_nx+1,d5_3d_ny+1,d5_3d_nz+1))*gauss

  d5_3d_id = addnewbgrddataset(dx=d5_3d_dx ,dy=d5_3d_dx ,zlength=d5_3d_zlen ,bx=d5_3d_bx, by=d5_3d_by ,bz=d5_3d_bz)



# Starting and ending position of first ideal D5 dipole.
# Also used to define the lattice bend

#d5p1_zs = 69.2 # may need to be revised upon obtaining lattice design data
#d5p1_ze = d5p1_zc + (d5p1_zc - d5p1_zs)

# Define starting and ending position of 1st ideal D5 dipole using ideal length and centre position
# Makes the bend less tight than when the starting position is set to be 69.2 m

# --- define dipole d5 

# --- ideal (uniform) field 
if d5p1_typ == "ideal":
  d5p1_zs = d5p1_zc - d5p1_ideal_len / 2.
  d5p1_ze = d5p1_zc + d5p1_ideal_len / 2.
  bending_R = (d5p1_ze - d5p1_zs)/(pi/2.)
  bending_B = sqrt( A_ref*ekin_per_u*jperev*2.*A_ref*amu)/(Q_ref*jperev)/bending_R
  d5p1 = addnewdipo(zs = d5p1_zs, ze = d5p1_ze, by = bending_B)
  equivalent_G = 0.9/bending_R**2 * ref_brho_post_gap	#focusing effect from slanted poles
  d5p1_equiv_quad = addnewquad(zs = d5p1_zs, ze = d5p1_ze, db = -equivalent_G)
# --- linear optic approximation field 
elif d5p1_typ == "lin":
  print("Warning: No D5 1st Dipole Linear Applied Fields Defined")
  d5p1 = None
# --- 3D field from magnet design code
elif d5p1_typ == "nl":
  #d5p1 = addnewbgrd(dx=d5_3d_dx,dy=d5_3d_dy,xs=d5_3d_x_m.min(),ys=d5_3d_y_m.min(),
  #  zs=d5p1_zc-d5_3d_zlen/2.,ze=d5p1_zc+d5_3d_zlen/2.,id=d5_3d_id,sc=d5p1_str)

  d5p1_zs = d5p1_zc - d5_3d_core_l/2. # core starts
  d5p1_ze = d5p1_zc + d5_3d_core_l/2. # core ends
  d5p1_str = d5p1_zs - d5_3d_frng_l   # fringe starts
  d5p1_end = d5p1_ze + d5_3d_frng_l   # fringe ends

  d5p1 = addnewbgrd(dx=d5_3d_dx, dy=d5_3d_dy, xs=-d5_3d_xw/2., ys=-d5_3d_yw/2.,
    zs=d5p1_str, ze=d5p1_end, id=d5_3d_id, sc=d5_3d_scl)

else:
  print("Warning: No D5 1st Dipole Applied Fields Defined") 
  d5p1 = None

# --- ideal (uniform) field 
if d5p2_typ == "ideal": 
  d5p2_zs = d5p2_zc - d5p2_ideal_len / 2.
  d5p2_ze = d5p2_zc + d5p2_ideal_len / 2.
  bending_R = (d5p2_ze - d5p2_zs)/(pi/2.)
  bending_B = sqrt( A_ref*ekin_per_u*jperev*2.*A_ref*amu)/(Q_ref*jperev)/bending_R
  d5p2 = addnewdipo(zs = d5p2_zs, ze = d5p2_ze, by = bending_B)
  equivalent_G = 0.9/bending_R**2 * ref_brho_post_gap	#focusing effect from slanted poles
  d5p2_equiv_quad = addnewquad(zs = d5p2_zs, ze = d5p2_ze, db = -equivalent_G)
# --- linear optic approximation field 
elif d5p2_typ == "lin":
  print("Warning: No D5 2nd Dipole Linear Applied Fields Defined")
  d5p2 = None
# --- 3D field from magnet design code
elif d5p2_typ == "nl":
  #d5p1 = addnewbgrd(dx=d5_3d_dx,dy=d5_3d_dy,xs=d5_3d_x_m.min(),ys=d5_3d_y_m.min(),
  #  zs=d5p1_zc-d5_3d_zlen/2.,ze=d5p1_zc+d5_3d_zlen/2.,id=d5_3d_id,sc=d5p1_str)

  d5p2_zs = d5p2_zc - d5_3d_core_l/2.  # core starts
  d5p2_ze = d5p2_zc + d5_3d_core_l/2.  # core ends
  d5p2_str = d5p2_zs - d5_3d_frng_l    # fringe starts
  d5p2_end = d5p2_ze + d5_3d_frng_l    # fringe ends
  
  d5p2 = addnewbgrd(dx=d5_3d_dx, dy=d5_3d_dy, xs=-d5_3d_xw/2., ys=-d5_3d_yw/2.,
    zs=d5p2_str, ze=d5p2_end, id=d5_3d_id, sc=d5_3d_scl)

else:
  print("Warning: No D5 2nd Dipole Applied Fields Defined") 
  d5p2 = None

#
# Lattice bends for D5 Bending Dipole 
#
#  Comments:
#   * Use ideal bends on lattice whether dipole field is ideal (uniform) or not.   

d5p1_bend = True  # True or False: Add ideal bend to lattice 
d5p2_bend = True  # True or False: Add ideal bend to lattice 


if d5p1_bend:
  top.diposet = False     # turn off By that automatically generated with addnewbend()
  equivalent_ideal_R = (d5p1_ze - d5p1_zs)/(pi/2.)
  equivalent_ideal_B = sqrt( A_ref*ekin_per_u*jperev*2.*A_ref*amu)/(Q_ref*jperev)/equivalent_ideal_R
  addnewbend(zs = d5p1_zs, ze = d5p1_ze, rc = equivalent_ideal_R)

if d5p2_bend:
  top.diposet = False     # turn off By that automatically generated with addnewbend()
  equivalent_ideal_R = (d5p2_ze - d5p2_zs)/(pi/2.)
  equivalent_ideal_B = sqrt( A_ref*ekin_per_u*jperev*2.*A_ref*amu)/(Q_ref*jperev)/equivalent_ideal_R
  addnewbend(zs = d5p2_zs, ze = d5p2_ze, rc = equivalent_ideal_R)


dipole_exit = [[0.,0.],[0.,0.]]
cccounter = 0

for ii in sp_target:
	s = sp[ii]
	species_R = sqrt( s.charge*Bias*2.*s.mass)/(s.charge)/ equivalent_ideal_B
	offset = sqrt( species_R**2 - (species_R - equivalent_ideal_R)**2 ) - equivalent_ideal_R
	centroid_angle = arccos( 1 - equivalent_ideal_R / species_R ) - pi/2
	
	dipole_exit[cccounter][0] = offset
	dipole_exit[cccounter][1] = centroid_angle*sign(offset)
	
	cccounter += 1






















# Q7 Electrostatic Quads
# Comment: linear and nonlinear variants must have same z-grid. 

# --- calculate esq strength from k1 in DIMAD lattice design
# values used in Dr. Ren's model: 4.6kV, -8.3kV, 3.7kV

amu_eV = 931.4941e6

q7t1p1_k1 = 8.2957122
q7t1p2_k1 = -15.6040684
q7t1p3_k1 = 7.51275015

q7_design_len = 0.2068

q7_aper_r = 7.5*cm

#gamma_ref = (ekin_per_u + amu_eV) / amu_eV
#v_ref = sqrt(1. - 1./gamma_ref**2)*clight
#brho_ref = (A_ref * sqrt((ekin_per_u + amu_eV)**2 - (amu_eV)**2)/clight) / Q_ref

q7_str_mode = 1    # 0: electrode potential corresponds to kappa ; 1: equivalent focusing
                   # must be 0 for ideal q7

inter_quad_distance = 0.335 # centroid distance between two quads in a triplet

# --- element specification 

## 1st triplet
q7t1p1_zc = 70.537759 # (q7: Q7 device type; t1: 1st triplet; p1: part 1)
if IdealCSS == True:
	q7t1p1_zc = d5p1_zc + 0.95
#q7t1p1_str = 10000 # [V]
q7t1p1_sign = 1    # +1 for x_quad, -1 for y_quad
q7t1p1_typ = "nl"  # type: "lin" = linear optics fields or "nl" = nonlinear r-z field

q7t1p2_zc = q7t1p1_zc + inter_quad_distance
#q7t1p2_str = 10000 # [V]
q7t1p2_sign = -1   # +1 for x_quad, -1 for y_quad
q7t1p2_typ = "nl"  # type: "lin" = linear optics fields or "nl" = nonlinear r-z field

q7t1p3_zc = q7t1p2_zc + inter_quad_distance
#q7t1p3_str = 10000 # [V]
q7t1p3_sign = 1    # +1 for x_quad, -1 for y_quad
q7t1p3_typ = "nl"  # type: "lin" = linear optics fields or "nl" = nonlinear r-z field  

## 2nd triplet
q7t2p1_zc = 72.161896 # (q7: Q7 device type; t2: 2nd triplet; p1: part 1)
if IdealCSS == True:
	q7t2p1_zc = d5p1_zc + 2.58
#q7t1p1_str = 10000 # [V]
q7t2p1_sign = 1    # +1 for x_quad, -1 for y_quad
q7t2p1_typ = "nl"  # type: "lin" = linear optics fields or "nl" = nonlinear r-z field

q7t2p2_zc = q7t2p1_zc + inter_quad_distance
#q7t1p2_str = 10000 # [V]
q7t2p2_sign = -1   # +1 for x_quad, -1 for y_quad
q7t2p2_typ = "nl"  # type: "lin" = linear optics fields or "nl" = nonlinear r-z field

q7t2p3_zc = q7t2p2_zc + inter_quad_distance
#q7t1p3_str = 10000 # [V]
q7t2p3_sign = 1    # +1 for x_quad, -1 for y_quad
q7t2p3_typ = "nl"  # type: "lin" = linear optics fields or "nl" = nonlinear r-z field  

## --- linear element data  
##fi = PRpickle.PR("lat_q7/lat_q7.lin.????.pkl")
#fi = PRpickle.PR("lat_q7/lat_q7.lin.????.pkl")
#s4_dz  = fi.s4_dz 
#s4_nz  = fi.s4_nz  
#s4_z_m = fi.s4_z_m 
#s4_bz0_m   = fi.s4_bz0_m
#s4_bz0p_m  = fi.s4_bz0p_m
#fi.close() 

#s4_zlen = s4_z_m.max() - s4_z_m.min() 
#s4_lin_id = addnewmmltdataset(zlen=s4_zlen,ms=s4_bz0_m,msp=s4_bz0p_m,nn=0,vv=0)

# --- nonlinear element field data 
fi = PRpickle.PR('lat_q7/lat_q7.3d.20160607.pkl') 
##
#s4_len_coil   = fi.s4_len_coil 
#s4_len_magnet = fi.s4_len_magnet 
#s4_r_coil_i   = fi.s4_r_coil_i 
#s4_r_coil_o   = fi.s4_r_coil_o
#
#if fi.s4_nz != s4_nz: raise Exception("S4: Nonlinear field model nz not equal to linear field model nz") 

q7_dx   = fi.q7_dx
q7_dy   = fi.q7_dy
q7_dz   = fi.q7_dz
q7_nx   = fi.q7_nx
q7_ny   = fi.q7_ny
q7_nz   = fi.q7_nz
q7_x_m  = fi.q7_x_m 
q7_y_m  = fi.q7_y_m 
q7_z_m  = fi.q7_z_m 
q7_ex_m  = fi.q7_ex_m 
q7_ey_m  = fi.q7_ey_m 
q7_ez_m  = fi.q7_ez_m 
fi.close()

q7_zlen = q7_z_m.max() - q7_z_m.min()
q7_x_m_min = q7_x_m.min()
q7_y_m_min = q7_y_m.min()

q7_nl_id = addnewegrddataset(dx=q7_dx,dy=q7_dy,zlength=q7_zlen,ex=q7_ex_m,ey=q7_ey_m,ez=q7_ez_m)  # pass arb dy to avoid error trap  



if q7t1p1_typ == "ideal":
    q7t1p1_str = abs(q7t1p1_k1)*ref_brho_post_gap*ref_vel_post_gap
    q7t1p2_str = abs(q7t1p2_k1)*ref_brho_post_gap*ref_vel_post_gap
    q7t1p3_str = abs(q7t1p3_k1)*ref_brho_post_gap*ref_vel_post_gap

if q7t1p1_typ != "ideal" and q7_str_mode == 0:
    q7t1p1_str = abs(q7t1p1_k1)*ref_brho_post_gap*ref_vel_post_gap*q7_aper_r**2/2
    q7t1p2_str = abs(q7t1p2_k1)*ref_brho_post_gap*ref_vel_post_gap*q7_aper_r**2/2
    q7t1p3_str = abs(q7t1p3_k1)*ref_brho_post_gap*ref_vel_post_gap*q7_aper_r**2/2

if q7t1p1_typ != "ideal" and q7_str_mode == 1:
	dEx_array = q7_ex_m[51][50] - q7_ex_m[49][50]  # difference of E at two sets of mesh points along z
	dExdx_array = dEx_array / (2*q7_dx)            # approximate dExdx at the centre
	normalized_sum = sum(dExdx_array)*q7_dz        
	design_Gl_p1 = q7t1p1_k1*ref_brho_post_gap*ref_vel_post_gap*q7_design_len
	design_Gl_p2 = q7t1p2_k1*ref_brho_post_gap*ref_vel_post_gap*q7_design_len
	design_Gl_p3 = q7t1p3_k1*ref_brho_post_gap*ref_vel_post_gap*q7_design_len
	q7t1p1_str = abs(design_Gl_p1/normalized_sum)
	q7t1p2_str = abs(design_Gl_p2/normalized_sum)
	q7t1p3_str = abs(design_Gl_p3/normalized_sum)

q7t2p1_str = q7t1p3_str
q7t2p2_str = q7t1p2_str
q7t2p3_str = q7t1p1_str



# --- define esq q7 1st triplet part 1 
if q7t1p1_typ == "lin":
  q7t1p1 = addnewmmlt(zs=q7t1p1_zc-q7_zlen/2.,ze=q7t1p1_zc+q7_zlen/2.,id=q7_lin_id,sc=q7t1p1_str*q7t1p1_sign) 
elif q7t1p1_typ == "nl":
  q7t1p1 = addnewegrd(xs=q7_x_m_min,ys=q7_y_m_min,zs=q7t1p1_zc-q7_zlen/2.,ze=q7t1p1_zc+q7_zlen/2.,id=q7_nl_id,sc=q7t1p1_str*q7t1p1_sign)
elif q7t1p1_typ == "ideal":   # additional -ve sign in de because dEx/dx should be -ve for focusing in x
  q7t1p1 = addnewquad(zs = q7t1p1_zc-q7_design_len/2., ze = q7t1p1_zc+q7_design_len/2., de = -q7t1p1_str*q7t1p1_sign)
else:
  print("Warning: No S4 1st Solenoid Applied Fields Defined") 
  q7t1p1 = None

# --- define esq q7 1st triplet part 2
if q7t1p2_typ == "lin":
  q7t1p2 = addnewmmlt(zs=q7t1p2_zc-q7_zlen/2.,ze=q7t1p2_zc+q7_zlen/2.,id=q7_lin_id,sc=q7t1p2_str*q7t1p2_sign) 
elif q7t1p2_typ == "nl":
  q7t1p2 = addnewegrd(xs=q7_x_m_min,ys=q7_y_m_min,zs=q7t1p2_zc-q7_zlen/2.,ze=q7t1p2_zc+q7_zlen/2.,id=q7_nl_id,sc=q7t1p2_str*q7t1p2_sign) 
elif q7t1p2_typ == "ideal":
  q7t1p2 = addnewquad(zs = q7t1p2_zc-q7_design_len/2., ze = q7t1p2_zc+q7_design_len/2., de = -q7t1p2_str*q7t1p2_sign)
else:
  print("Warning: No S4 1st Solenoid Applied Fields Defined") 
  q7t1p2 = None

# --- define esq q7 1st triplet part 3
if q7t1p3_typ == "lin":
  q7t1p3 = addnewmmlt(zs=q7t1p3_zc-q7_zlen/2.,ze=q7t1p3_zc+q7_zlen/2.,id=q7_lin_id,sc=q7t1p3_str*q7t1p3_sign) 
elif q7t1p3_typ == "nl":
  q7t1p3 = addnewegrd(xs=q7_x_m_min,ys=q7_y_m_min,zs=q7t1p3_zc-q7_zlen/2.,ze=q7t1p3_zc+q7_zlen/2.,id=q7_nl_id,sc=q7t1p3_str*q7t1p3_sign)
elif q7t1p3_typ == "ideal":
  q7t1p3 = addnewquad(zs = q7t1p3_zc-q7_design_len/2., ze = q7t1p3_zc+q7_design_len/2., de = -q7t1p3_str*q7t1p3_sign) 
else:
  print("Warning: No S4 1st Solenoid Applied Fields Defined") 
  q7t1p3 = None

# --- define esq q7 2nd triplet part 1 
if q7t2p1_typ == "lin":
  q7t2p1 = addnewmmlt(zs=q7t2p1_zc-q7_zlen/2.,ze=q7t2p1_zc+q7_zlen/2.,id=q7_lin_id,sc=q7t2p1_str*q7t2p1_sign) 
elif q7t2p1_typ == "nl":
  q7t2p1 = addnewegrd(xs=q7_x_m_min,ys=q7_y_m_min,zs=q7t2p1_zc-q7_zlen/2.,ze=q7t2p1_zc+q7_zlen/2.,id=q7_nl_id,sc=q7t2p1_str*q7t2p1_sign)
elif q7t2p1_typ == "ideal":
  q7t2p1 = addnewquad(zs = q7t2p1_zc-q7_design_len/2., ze = q7t2p1_zc+q7_design_len/2., de = -q7t2p1_str*q7t2p1_sign) 
else:
  print("Warning: No S4 1st Solenoid Applied Fields Defined") 
  q7t2p1 = None

# --- define esq q7 2nd triplet part 2
if q7t2p2_typ == "lin":
  q7t2p2 = addnewmmlt(zs=q7t2p2_zc-q7_zlen/2.,ze=q7t2p2_zc+q7_zlen/2.,id=q7_lin_id,sc=q7t2p2_str*q7t2p2_sign) 
elif q7t2p2_typ == "nl":
  q7t2p2 = addnewegrd(xs=q7_x_m_min,ys=q7_y_m_min,zs=q7t2p2_zc-q7_zlen/2.,ze=q7t2p2_zc+q7_zlen/2.,id=q7_nl_id,sc=q7t2p2_str*q7t2p2_sign)
elif q7t2p2_typ == "ideal":
  q7t2p2 = addnewquad(zs = q7t2p2_zc-q7_design_len/2., ze = q7t2p2_zc+q7_design_len/2., de = -q7t2p2_str*q7t2p2_sign)  
else:
  print("Warning: No S4 1st Solenoid Applied Fields Defined") 
  q7t2p2 = None

# --- define esq q7 2nd triplet part 3
if q7t2p3_typ == "lin":
  q7t2p3 = addnewmmlt(zs=q7t2p3_zc-q7_zlen/2.,ze=q7t2p3_zc+q7_zlen/2.,id=q7_lin_id,sc=q7t2p3_str*q7t2p3_sign) 
elif q7t2p3_typ == "nl":
  q7t2p3 = addnewegrd(xs=q7_x_m_min,ys=q7_y_m_min,zs=q7t2p3_zc-q7_zlen/2.,ze=q7t2p3_zc+q7_zlen/2.,id=q7_nl_id,sc=q7t2p3_str*q7t2p3_sign) 
elif q7t2p3_typ == "ideal":
  q7t2p3 = addnewquad(zs = q7t2p3_zc-q7_design_len/2., ze = q7t2p3_zc+q7_design_len/2., de = -q7t2p3_str*q7t2p3_sign)  
else:
  print("Warning: No S4 1st Solenoid Applied Fields Defined") 
  q7t2p3 = None


### Scrapers ###

## Particles are not scraped if scraper is too thin w.r.t. simulation step size ##
# From trial and error: use 4mm-thick slits, 10mm-thick walls for 2mm step size



## Beam pipe (where the aperture is circular)

r_p_up   = 8.00*cm  # aperture   upstream of grated gap [m]
r_p_down = 7.62*cm  # aperture downstream of grated gap [m]

post_d5p1_pipe_r = 7.5*cm
post_d5p1_pipe_zs = d5p1_ze + 8*cm
post_d5p1_pipe_ze = d5p1_ze + 8*cm +  229*mm

q7t1_pipe_r = 12.4*cm
q7t1_pipe_zs = post_d5p1_pipe_ze
q7t1_pipe_ze = q7t1_pipe_zs + 952*mm

r_ap   = array([r_p_up,               gag_rp,       r_p_down,   post_d5p1_pipe_r,     12.4*cm,       7.5*cm  ])
v_ap   = array([SourceBias+StandBias, StandBias/2., 0.,         0.,                   0.,            0.]) 
z_ap_l = array([ecr_z_extr,           gag_col_zs,   gag_col_ze, post_d5p1_pipe_zs,    q7t1_pipe_zs,  q7t1_pipe_ze])
z_ap_u = array([gag_col_zs,           gag_col_ze,   d5p1_zs,    post_d5p1_pipe_ze,    q7t1_pipe_ze,  d5p2_ze + 0.4])

beampipe = [] 
for i in range(len(r_ap)):
  rp = r_ap[i] 
  v  = v_ap[i] 
  zl = z_ap_l[i] 
  zu = z_ap_u[i]
  #
  beampipe.append( ZCylinderOut(radius=rp,zlower=zl,zupper=zu,voltage=v,condid="next") )

r_p = max(r_ap)   # Max aperture in simulations 


## Rectangular aperture in D5 Dipole:

# Include an additional 8cm of rectangular aperture after 1m long D5 Dipole

added_len = 8*cm

d5p1_aperture_out = Box(xsize = 180*mm, ysize = 150*mm, zsize = d5p1_ideal_len + added_len, xcent = 0*mm, ycent = 0., zcent = d5p1_zc + added_len/2)
d5p1_aperture_in = Box(xsize = 150*mm, ysize = 120*mm, zsize = d5p1_ideal_len + added_len, xcent = 0*mm, ycent = 0., zcent = d5p1_zc + added_len/2)

d5p1_aperture = [d5p1_aperture_out - d5p1_aperture_in]

#d5p1_aperture_xplus = Box(xsize = 10*mm, ysize = 120*mm, zsize = d5p1_ideal_len + added_len, xcent = 80*mm, ycent = 0, zcent = d5p1_zc + added_len/2)
#d5p1_aperture_xminus = Box(xsize = 10*mm, ysize = 120*mm, zsize = d5p1_ideal_len + added_len, xcent = -80*mm, ycent = 0, zcent = d5p1_zc + added_len/2)
#d5p1_aperture_yplus = Box(xsize = 150*mm, ysize = 10*mm, zsize = d5p1_ideal_len + added_len, xcent = 0, ycent = 65*mm, zcent = d5p1_zc + added_len/2)
#d5p1_aperture_yminus = Box(xsize = 150*mm, ysize = 10*mm, zsize = d5p1_ideal_len + added_len, xcent = 0, ycent = -65*mm, zcent = d5p1_zc + added_len/2)

#d5p1_aperture = [d5p1_aperture_xplus, d5p1_aperture_xminus, d5p1_aperture_yplus, d5p1_aperture_yminus]

## Gate valve

valve_x_opening = 7*cm
valve_y_opening = 10*cm
valve_len = 4*mm
valve_zc = q7t1p1_zc - 140*mm

# These sizes are arbitrary, 
valve_xsize = 7*cm
valve_ysize = 10*cm


valve_x_plus = Box(xsize = valve_xsize, ysize = valve_ysize, zsize = valve_len, xcent = (valve_x_opening + valve_xsize)/2, ycent = 0, zcent = valve_zc)
valve_x_minus = Box(xsize = valve_xsize, ysize = valve_ysize, zsize = valve_len, xcent = -(valve_x_opening + valve_xsize)/2, ycent = 0, zcent = valve_zc)
valve_y_plus = Box(xsize = valve_xsize, ysize = valve_ysize, zsize = valve_len, xcent = 0, ycent = (valve_y_opening + valve_ysize)/2, zcent = valve_zc)
valve_y_minus = Box(xsize = valve_xsize, ysize = valve_ysize, zsize = valve_len, xcent = 0, ycent = -(valve_y_opening + valve_ysize)/2, zcent = valve_zc)

gate_valve = [valve_x_plus, valve_x_minus, valve_y_plus, valve_y_minus]

## End plates of the ESQs in the 1st triplet

endplate_len = 4*mm
endplate_aperture = 65*mm

q7t1_endplate_1 = ZCylinderOut(endplate_aperture, endplate_len, zcent= q7t1p1_zc - 19.5*mm)
q7t1_endplate_2 = ZCylinderOut(endplate_aperture, endplate_len, zcent= q7t1p1_zc + 19.5*mm)
q7t1_endplate_3 = ZCylinderOut(endplate_aperture, endplate_len, zcent= q7t1p2_zc - 19.5*mm)
q7t1_endplate_4 = ZCylinderOut(endplate_aperture, endplate_len, zcent= q7t1p2_zc + 19.5*mm)
q7t1_endplate_5 = ZCylinderOut(endplate_aperture, endplate_len, zcent= q7t1p3_zc - 19.5*mm)
q7t1_endplate_6 = ZCylinderOut(endplate_aperture, endplate_len, zcent= q7t1p3_zc + 19.5*mm)

q7t1_endplates = [q7t1_endplate_1,q7t1_endplate_2,q7t1_endplate_3,q7t1_endplate_4,q7t1_endplate_5,q7t1_endplate_6]

## Slits between ESQ in the 1st triplet

# Slit sizes are arbitrarily chosen at present, only making sure they extend into the beam tube in transverse directions
# The opening size in x is based on drawings

q7t1_mid_12 = (q7t1p1_zc+q7t1p2_zc)/2  # mid-point between 1st and 2nd ESQ 
q7t1_mid_23 = (q7t1p2_zc+q7t1p3_zc)/2  # mid-point between 2nd and 3rd ESQ 

slits1_x_opening = 7*cm
slits2_x_opening = 9*cm
slit_xsize = 8*cm
slit_ysize = 15*cm
slit_len = 4*mm

q7t1_slits1_xplus = Box(xsize = slit_xsize, ysize = slit_ysize, zsize = slit_len, xcent = (slits1_x_opening + slit_xsize)/2, ycent = 0, zcent = q7t1_mid_12)
q7t1_slits1_xminus = Box(xsize = slit_xsize, ysize = slit_ysize, zsize = slit_len, xcent = -(slits1_x_opening + slit_xsize)/2, ycent = 0, zcent = q7t1_mid_12)
q7t1_slits2_xplus = Box(xsize = slit_xsize, ysize = slit_ysize, zsize = slit_len, xcent = (slits2_x_opening + slit_xsize)/2, ycent = 0, zcent = q7t1_mid_23)
q7t1_slits2_xminus = Box(xsize = slit_xsize, ysize = slit_ysize, zsize = slit_len, xcent = -(slits2_x_opening + slit_xsize)/2, ycent = 0, zcent = q7t1_mid_23)

q7t1_slits = [q7t1_slits1_xplus,q7t1_slits1_xminus,q7t1_slits2_xplus,q7t1_slits2_xminus]

# Q7 Electrodes

# function returns a list of 4 cylinders that approximate the q7 electrodes
# required input: center position of the electrode

def q7_electrodes_conductor(zcenter):
	q7_electrode_len = 200*mm
	
	# approximate each electrode by a cylinder
	cylinder_offset = 170*mm
	cylinder_radius = 95*mm
	
	electrode_xplus = ZCylinder(cylinder_radius, q7_electrode_len, xcent = cylinder_offset, zcent= zcenter)
	electrode_xminus = ZCylinder(cylinder_radius, q7_electrode_len, xcent = -cylinder_offset, zcent= zcenter)
	electrode_yplus = ZCylinder(cylinder_radius, q7_electrode_len, ycent = cylinder_offset, zcent= zcenter)
	electrode_yminus = ZCylinder(cylinder_radius, q7_electrode_len, ycent = -cylinder_offset, zcent= zcenter)
	
	return [electrode_xplus, electrode_xminus, electrode_yplus, electrode_yminus]






#scraper = ParticleScraper([q7t1_endplate_1,q7t1_endplate_2,q7t1_endplate_3,q7t1_endplate_4,q7t1_endplate_5,q7t1_endplate_6])

scraperlist = beampipe + d5p1_aperture + gate_valve + q7t1_endplates + q7t1_slits + q7_electrodes_conductor(q7t1p1_zc) + q7_electrodes_conductor(q7t1p2_zc) + q7_electrodes_conductor(q7t1p3_zc)

scraper = ParticleScraper(scraperlist)

#scraper.registerconductor(scaperlist)
















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
 gag_zc + 22.28*cm,   # z of neut stop after  grated gap, set where 1% of gap E_z field reached
 d5p1_ze              # z of where the 1st D5 dipole would ideally end   
]    )

sp_neut_frac = \
array(
[0.75, 
 0., 
 0.75,
 0.
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

#r_p_up   = 8.00*cm  # aperture   upstream of grated gap [m]
#r_p_down = 7.62*cm  # aperture downstream of grated gap [m]

#r_ap   = array([r_p_up,               gag_rp,       r_p_down  ])
#v_ap   = array([SourceBias+StandBias, StandBias/2., 0.        ]) 
#z_ap_l = array([ecr_z_extr,           gag_col_zs,   gag_col_ze])
#z_ap_u = array([gag_col_zs,           gag_col_ze,   d5p1_zc   ])

#r_p = max(r_ap)   # Max aperture in simulations 

#def aperture_r(z):
  #index = sum(where(z/z_ap_l >= 1.,1,0))-1 
  #if index < 0: index = 0
  ## 
  #return(r_ap[index])  

#aperture = [] 
#for i in range(len(r_ap)):
  #rp = r_ap[i] 
  #v  = v_ap[i] 
  #zl = z_ap_l[i] 
  #zu = z_ap_u[i]
  ##
  #aperture.append( ZCylinderOut(radius=rp,zlower=zl,zupper=zu,condid="next") )

## --- Add a circular aperture particle scraper at pipe radius aperture. 
##       Could also use ParticleScraper(conductors=aperture) but 
##       setting prwall is faster for a simple cylinder. 
#top.prwall = r_p    # reset this later consistent with actual aperture in range simulated in advances 

