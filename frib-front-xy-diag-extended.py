# Diagnostics for FRIB front-end simulations using the x-y slice model 
# Notes:
#  * Slice model is not intrinsically well adapted to multi-species 
#    simulations so some diagnostics repeat (for clarity) what can be 
#    generated within Warp with other methods.
#  * Model allows easy generalization to include diagnostic quantities not 
#    in the usual Warp suite.  

##############################################################################
# Begin Inputs 
##############################################################################

# Diagnostic Parameters 
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
# --- Plot limits for particle phase space plots. If lframe = true (default
#     false) diagnostics such as ppxxp for x-x' particle phase space will
#     use these ranges.  
#      max/min x,y   plot coordinates (m) 
#      max/min x',y' plot coordinates (rad)
#l_diag = r_p
l_diag = 75*mm
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


# Particle phase space diagnostics.
# * The list diag_step_part contains all steps where diagnostics in
#   diag_part() are made.  
# * The list can contain repeated elements and need not be ordered.

diag_part_z = array([
  z_launch,
  s4p1_zc,
  s4p1_zc-20.*cm,
  s4p1_zc+20.*cm, 
  s4p2_zc,
  s4p2_zc-20.*cm,
  s4p2_zc+20.*cm, 
  gag_zc, 
  #gag_col_zs,
  #gag_col_zs-5.*cm,
  #gag_col_ze,
  gag_col_ze+5.*cm,
  d5p1_zs,
  (d5p1_zs+d5p1_zc)/2,
  d5p1_zc,
  (d5p1_zc+d5p1_ze)/2,
  d5p1_ze,
  valve_zc + 2*mm,
  q7t1p1_zc,
  q7t1_mid_12 + 2*mm,
  q7t1p2_zc,
  q7t1_mid_23 + 2*mm,
  q7t1p3_zc,
  z_adv]) 

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
  "D5 Dipole #1: z-start",
  "D5 Dipole #1: 1/4 of dipole length",
  "D5 Dipole #1: z-Center",
  "D5 Dipole #1: 3/4 of dipole length",
  "D5 Dipole #1: z-end",
  "after gate valve",
  "1st Q7 ESQ Triplet #1: z-Center",
  "after slits between Q7 #1 and #2",
  "1st Q7 ESQ Triplet #2: z-Center",
  "after slits between Q7 #2 and #3",
  "1st Q7 ESQ Triplet #3: z-Center",
  "Final position" 
                   ]      

diag_part_step = nint((diag_part_z-z_launch)/wxy.ds)

diag_part_z_names = {diag_part_step[i]:diag_part_z_name[i] for i in range(len(diag_part_step))}


# Field diagnostics.  
# * The list diag_step_field containins all steps where
#   diagnostics in diag_field() are made. 
# * The list can contain repeated elements and need not be ordered.   

diag_field_z = array([
  z_launch,
  s4p1_zc,gag_zc,
  s4p2_zc,
  d5p1_zc,
  z_adv 
                    ]) 

diag_field_z_name = [ 
  "Initial Launch", 
  "S4 Solenoid #1: z-Center", 
  "Grated Gap: z-Center",
  "S4 Solenoid #1: z-Center",
  "D5 Dipole #1: z-Center", 
  "Final position"
                     ]

diag_field_step = nint((diag_field_z-z_launch)/wxy.ds)

diag_field_z_names = {diag_field_step[i]:diag_field_z_name[i] for i in range(len(diag_field_step))}

# History diagnostics.  
# * Can be made at intermediate stages of the
#   run as well as at the end.  
# * The list diag_step_hist contains all
#   steps where diagnostics in diag_hsit() are made. 
# * The list can contain repeated elements and need not be ordered.

diag_hist_z    = array([z_adv]) #array([gag_col_zs,z_adv])
diag_hist_step = nint((diag_hist_z-z_launch)/wxy.ds)


######################################################################################################
# End Inputs 
######################################################################################################


# Diagnostic plot function of [B rho] vs Q/A for species.  
#  * Should work correctly at any point in the simulation while the beam 
#    accelerates.  

def plt_diag_bro(label=None):
  if label == None: label = " "
  brho_min =  largepos 
  brho_max = -largepos  
  for ii in sp.keys():
    s = sp[ii]
    js = s.js 
    #
    weight = sum(s.sw*s.w)   # total weight 
    #
    vbeam = sum( (s.sw*s.w)*s.getvz() )/weight  # avg axial velocity 
    gammabeam = 1./sqrt(1.-(vbeam/clight)**2)   # gamma from avg axial velocity 
    brho  = s.mass*gammabeam*vbeam/s.charge     # rigidity 
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




# Potential profile plot diagnostic for potential along x-y axes 
# * Primarily for initial beam but should work at any point in simulation. 

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


#  Augmented History Diagnostics for xy Slice Model 
# * Some by species, some all species 
# * Flag variables with prefix hl_ for "history local" 

# --- History variable accumulation arrays

hl_lenhist_max = 10000 # max accumulation points
#
hl_zbeam    = fzeros(hl_lenhist_max)           # z of beam at hl_ diagnostic accumulations (redundant with top.hzbeam)
#
hl_vbeam    = fzeros([hl_lenhist_max,top.ns])  # axial beam velocity [m/s]
hl_ekin     = fzeros([hl_lenhist_max,top.ns])  # axial beam NR kinetic energy [eV]
hl_brho     = fzeros([hl_lenhist_max,top.ns])  # rigidity [B rho]_js [Tesla-m] 
#
hl_xrms     = fzeros([hl_lenhist_max,top.ns])  # rms radius sqrt( <x*x>_js )  
hl_yrms     = fzeros([hl_lenhist_max,top.ns])  # rms radius sqrt( <y*y>_js )  
hl_rrms     = fzeros([hl_lenhist_max,top.ns])  # rms radius sqrt( <r*r>_js )  
#
hl_xrmst    = fzeros(hl_lenhist_max)           # Total species measures of above  
hl_yrmst    = fzeros(hl_lenhist_max)           #  
hl_rrmst    = fzeros(hl_lenhist_max)           #
#
hl_spnum    = fzeros([hl_lenhist_max,top.ns])  # number active simulation particles
hl_spnumt   = fzeros(hl_lenhist_max)           # number active simulation particles (all species) 
#
hl_ibeam_p  = fzeros([hl_lenhist_max,top.ns])  # beam current (particle)   
hl_ibeam_e  = fzeros([hl_lenhist_max,top.ns])  # beam current (electrical) 
hl_ibeam_pt = fzeros([hl_lenhist_max])         # total beam current (particle)   
hl_ibeam_et = fzeros([hl_lenhist_max])         # total beam current (electrical)   
#
hl_lambda_p = fzeros([hl_lenhist_max,top.ns])  # line charge (particle) 
hl_lambda_e = fzeros([hl_lenhist_max,top.ns])  # line charge (electrical) 
#
hl_ptheta   = fzeros([hl_lenhist_max,top.ns])  # canonical angular momentum <P_theta>_j (nonlinear appl field version)   
hl_pth      = fzeros([hl_lenhist_max,top.ns])  # <P_theta>_j in emittance units <P_theta>_j/(gamma_j*beta_j*m_j*c)  
hl_pthn     = fzeros([hl_lenhist_max,top.ns])  # <P_theta>_j in norm emittance units <P_theta>_j/(m_j*c)
#
hl_ptheta_l = fzeros([hl_lenhist_max,top.ns])  # Same canonical angular momentum measures with
hl_pth_l    = fzeros([hl_lenhist_max,top.ns])  #   linear applied magnetic field approximation.  
hl_pthn_l   = fzeros([hl_lenhist_max,top.ns])  #   (redundant with above for linear lattice) 
#
hl_lz       = fzeros([hl_lenhist_max,top.ns])  # mechanical angular momentum
hl_krot     = fzeros([hl_lenhist_max,top.ns])  # rotation wavenumber  
hl_lang     = fzeros([hl_lenhist_max,top.ns])  # Larmor rotation angle (from initial zero value)   
#
hl_epsx     = fzeros([hl_lenhist_max,top.ns])  # rms x-emittance (usual version)
hl_epsy     = fzeros([hl_lenhist_max,top.ns])  # rms y-emittance (usual version) 
#
hl_epsxn    = fzeros([hl_lenhist_max,top.ns])  # rms normalized x-emittance (usual version) 
hl_epsyn    = fzeros([hl_lenhist_max,top.ns])  # rms normalized y-emittance (usual version)  
#
hl_epsr     = fzeros([hl_lenhist_max,top.ns])  # rms radial emittance               (envelope model version) 
hl_epsrn    = fzeros([hl_lenhist_max,top.ns])  # rms normalized radial emittance    (envelope model version)
#
hl_epspv    = fzeros([hl_lenhist_max,top.ns])  # rms total phase volume emittance             (envelope model sense)
hl_epspvn   = fzeros([hl_lenhist_max,top.ns])  # rms normalized total phase volume emittance  (envelope model sense)
#
hl_temp     = fzeros([hl_lenhist_max,top.ns])  # Effective transverse ion temperature measure [eV]
#
hl_Qperv    = fzeros([hl_lenhist_max,top.ns])  # Generalized perveance Q_js for species: note matrix perv 
                                               #   Q_js,s calculable from this and line-charge densities  [1]   
hl_neutf    = fzeros([hl_lenhist_max,top.ns])  # Neutralization factor [1] 

hl_dz = top.nhist*wxy.ds  # Axial step size between diagnostic accumulations 

# ---- Function to Fill Auxillary History Arrays 
#      * Install after step in particle advance cycle 

@callfromafterstep
def diag_hist_hl():
  # check step in history accumulation cycle 
  if top.it%top.nhist != 0: return
  hl_zbeam[top.jhist] = top.zbeam  # z location of diagnostic accumulations 
  # accumulate history diagnostics by species
  weightt_work = 0.  
  xrmst_work = 0. 
  yrmst_work = 0. 
  rrmst_work = 0. 
  for ii in sp.keys():
    # --- species info and index js 
    s = sp[ii]
    js = s.js 
    # --- species weight: (real particle per macroparticle)/meter 
    weight = sum(s.sw*s.w) 
    # --- <v_z>_js, gamma_js and [B rho]_js calculated from result 
    vbeam     = sum( (s.sw*s.w)*s.getvz() )/weight
    gammabeam = 1./sqrt(1.-(vbeam/clight)**2)      
    brho      = s.mass*gammabeam*vbeam/s.charge
    hl_vbeam[top.jhist,js] = vbeam 
    hl_brho[top.jhist,js]  = brho  
    #
    # --- species quantities for later use 
    # --- avg_rsq = <r*r>_js
    r   = s.getr() 
    rsq = r*r  
    rsq_wsum = sum( (s.sw*s.w)*rsq )
    avg_rsq = rsq_wsum/weight
    # --- avg_xyp = <x*y'>_js and avg_yxp = <y*x'>_js 
    avg_xyp = sum( (s.sw*s.w)*s.getx()*s.getyp() )/weight
    avg_yxp = sum( (s.sw*s.w)*s.gety()*s.getxp() )/weight
    # --- avg_xpy = <x*p_y>_js and avg_ypx = <y*p_x>_js
    #       * Relativistically correct here  
    avg_xpy = s.mass*sum( (s.sw*s.w)*s.getx()*s.getuy() )/weight   
    avg_ypx = s.mass*sum( (s.sw*s.w)*s.gety()*s.getux() )/weight   
    # --- applied field B_z(r=0,z) at z location of beam 
    bz0  = getappliedfields(x=0.,y=0.,z=top.zbeam)[5]
    # --- Axial kinetic energy [eV], ekin_js, NR calcuation  
    hl_ekin[top.jhist,js] = (0.5*s.mass*sum( (s.sw*s.w)*s.getvz()**2 )/weight)/jperev         
                            # s.mass*clight**2*(gammabeam - 1.)/jperev 
    # --- rms x = <x*x>_js
    xsq_wsum = sum( (s.sw*s.w)*s.getx()**2 )
    hl_xrms[top.jhist,js] = sqrt( xsq_wsum/weight )
    # --- rms y = <y*y>_js
    ysq_wsum = sum( (s.sw*s.w)*s.gety()**2 )
    hl_yrms[top.jhist,js] = sqrt( ysq_wsum/weight )
    # --- rms r = <r*r>_js 
    hl_rrms[top.jhist,js] = sqrt( avg_rsq )     
    # --- Simulation Particle Number 
    hl_spnum[top.jhist,js] = s.getn()  
    # --- Current, electrical, Ie_js  [A]
    hl_ibeam_e[top.jhist,js] = s.charge*sum( (s.sw*s.w)*s.getvz() )     # slice code weight is particles/meter 
    # --- Current, particle, Ip_js [A]
    #       * Use way to calculate to remove neutralization factor 
    #       * Formula as given approx (paraxial) using appropriate weights 
    hl_ibeam_p[top.jhist,js] = s.charge*s.sw*(s.vbeam0/vbeam)*sum( s.getvz() ) 
    # --- line charge Lambda_js 
    hl_lambda_p[top.jhist,js] = hl_ibeam_p[top.jhist,js]/vbeam 
    hl_lambda_e[top.jhist,js] = hl_ibeam_e[top.jhist,js]/vbeam 
    # --- Mechanical angular momentum: <x*y'>_js - <y*x'>_js  
    hl_lz[top.jhist,js] = avg_xyp - avg_yxp
    # --- Canonical angular momentum <P_theta>_js 
    #       Notes: * Uses A_theta via getatheata() consistently with linear/nonlinear elements.   
    hl_ptheta[top.jhist,js] = avg_xpy - avg_ypx + sum( (s.sw*s.w)*s.charge*r*getatheta(r) )/weight
    # --- Normalized canonical angular momentum in emittance units. <P_theta>_js/(m_js*c) 
    #       * <P_theta>_j/(m_j*c) in envelope model scales as a normalized emittance 
    #         and should not vary with acceleration with linear forces.    
    #       * This employs the nonlinear definition of P_theta if the lattice is nonlinear !    
    hl_pthn[top.jhist,js] = hl_ptheta[top.jhist,js]/(s.mass*clight)
    # --- Canonical angular momentum of species in emittance units 
    hl_pth[top.jhist,js] = hl_pthn[top.jhist,js]/(gammabeam*(vbeam/clight))
    # --- Canonical angular momentum in linear applied field approx (all 3 versions above) 
    #       * These are redundant in linear field lattice 
    #       * Use _l for "linear" flag 
    hl_ptheta_l[top.jhist,js] = avg_xpy - avg_ypx + sum( (s.sw*s.w)*(s.charge*bz0/2.)*avg_rsq )/weight
    hl_pthn_l[top.jhist,js]   = hl_ptheta_l[top.jhist,js]/(s.mass*clight)
    hl_pth_l[top.jhist,js]    = hl_pthn_l[top.jhist,js]/(gammabeam*(vbeam/clight))
    # --- rms x- and y-emittances: account for factor of 4 diff between Warp rms edge and rms measures 
    hl_epsx[top.jhist,js] = top.hepsx[0,top.jhist,js]/4.
    hl_epsy[top.jhist,js] = top.hepsy[0,top.jhist,js]/4.
    # --- normalized rms x- and y-emittances: paraxial equivalent version 
    hl_epsxn[top.jhist,js] = (gammabeam*(vbeam/clight))*hl_epsx[top.jhist,js]
    hl_epsyn[top.jhist,js] = (gammabeam*(vbeam/clight))*hl_epsy[top.jhist,js]
    # --- rms radial thermal emittance eps_r_js as derived in envelope model: 
    #       * Warp accumulation used to extract has a factor of 2 diference from rms envelope model 
    #         due to use of an "edge" measure.  Note: this is different than the factor of 4 in epsx etc.  
    hl_epsr[top.jhist,js] = top.hepsr[0,top.jhist,js]/2.  
    # --- rms normalized radial thermal emittance epsn_r_js as derived in envelope model 
    hl_epsrn[top.jhist,js] = (gammabeam*(vbeam/clight))*hl_epsr[top.jhist,js]
    # --- rms total phase volume emittance including radial thermal and canonical angular momentum 
    #       contributions based on envelope model intrpretation of total phase-space area. 
    hl_epspv[top.jhist,js] = sqrt( (hl_epsr[top.jhist,js])**2 + (hl_pth[top.jhist,js])**2 ) 
    # --- rms normalized total phase volume emittance
    hl_epspvn[top.jhist,js] = sqrt( (hl_epsrn[top.jhist,js])**2 + (hl_pthn[top.jhist,js])**2 )  
    # --- ion temperature calculated from emittance [eV]
    hl_temp[top.jhist,js] = hl_ekin[top.jhist,js]*hl_epsr[top.jhist,js]**2/dvnz(hl_rrms[top.jhist,js]**2) 
    # --- Perveance, NR formula for species
    #     Note: * Define bare ... not accounting for neutralization fractions.
    #             Factor (s.charge/echarge) = Q accounts for charge state with particle line-charge to 
    #             get bare (unneutralized) electrical line charge.  
    #           * This is Q_js NOT the matrix perveance Q_j,s in the envelope model notes. 
    #           * Envelope model Q_js can be obtained from Q_j and line charges lambda_j: no need to save 
    hl_Qperv[top.jhist,js] = s.charge*(s.charge/echarge)*hl_lambda_p[top.jhist,js]/(2.*pi*eps0*s.mass*vbeam**2)
    # --- Ion rho electron neutralization factor [1] = No space-charge, [0] full space-charge 
    hl_neutf[top.jhist,js] = rho_neut_f(top.zbeam,ii) 
    # --- Rotation wavenumber 
    hl_krot[top.jhist,js] = hl_lz[top.jhist,js]/dvnz(avg_rsq)
    # --- Larmor Rotation angle: integrate from previous step  
    if top.jhist == 0:
      hl_lang[0,js] = 0.  # initial condition of zero angle 
    else:
      hl_lang[top.jhist,js] = hl_lang[top.jhist-1,js] + 0.5*hl_dz*(hl_krot[top.jhist-1,js]+hl_krot[top.jhist,js])
    # --- total (all species) accumulations 
    weightt_work = weightt_work + weight 
    xrmst_work = xrmst_work + xsq_wsum 
    yrmst_work = yrmst_work + ysq_wsum
    rrmst_work = rrmst_work + rsq_wsum
  # --- total number of simulation particles 
  hl_spnumt[top.jhist] = float(sum(hl_spnum[top.jhist,:]))
  # --- total currents 
  hl_ibeam_pt[top.jhist] = sum(hl_ibeam_p[top.jhist,:]) 
  hl_ibeam_et[top.jhist] = sum(hl_ibeam_e[top.jhist,:]) 
  # --- total species rms measures 
  hl_xrmst[top.jhist] = sqrt( xrmst_work/weightt_work ) 
  hl_yrmst[top.jhist] = sqrt( yrmst_work/weightt_work ) 
  hl_rrmst[top.jhist] = sqrt( rrmst_work/weightt_work ) 
 

# Particle Phase-Space Diagnostic Functions 
# * Make specified plots at location of simulation where diag_part() is called.  

### Make phase space projections of individual species using the syntax "ppxxp(js=js)" instead of "s.ppxxp"
### The latter had trouble with the argument ' slope="auto" '

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
      js = s.js
      ppxy(js=js,lframe=true,chopped=chop_fraction,titles=false,yscale=1./mm,xscale=1./mm,color=co)
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
      js = s.js
      ppxxp(js=js,lframe=true,chopped=chop_fraction,slope='auto',titles=false,yscale=1./mr,xscale=1./mm,color=co)
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
      js = s.js
      ppyyp(js=js,lframe=true,chopped=chop_fraction,slope='auto',titles=false,yscale=1./mr,xscale=1./mm,color=co)
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
      js = s.js
      ppxpyp(js=js,lframe=true,chopped=chop_fraction,slope='auto',titles=false,yscale=1./mr,xscale=1./mm,color=co)
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
    plg(den/cm**3, rmesh/mm)  # pos axis 
    plg(den/cm**3,-rmesh/mm)  # neg axis 
    ptitles("Radial Number Density: All Species, z = %5.2f m"%(top.zbeam),"radius r [mm]","rho [particles/cm**3]",z_label)
    ir = min(nr,sum(where(den>0,1,0)))      # index farthest radial extent of rho in radial mesh assuming no halo  
    rmmax = max(1.2*rmesh[ir],0.01) # set curoff to contain radial density  
    rmmax = cm*nint(rmmax/cm + 0.5) # round up to nearest cm to contain plot 
    denmax = 1.2*maxnd(den) 
    limits(-rmmax/mm,rmmax/mm,0.,denmax/cm**3)
    fma() 
    # --- for all species (common log scale)  
    for ii in sp.keys():
       s  = sp[ii]
       co = s.color 
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
       den[0]      = den[1]   # set origin by next grid up to remove distraction (origin location high noise) 
       # 
       plg(den/cm**3, rmesh/mm,color=co)
       plg(den/cm**3,-rmesh/mm,color=co) 
    #
    ptitles("Radial Number Density: All species, z = %5.2f m"%(top.zbeam),"radius r [mm]","rho [particles/cm**3]",z_label)
    limits(-rmmax/mm,rmmax/mm,1.e-4*denmax/cm**3,denmax/cm**3)
    logxy(0,1)  # specify log scale on y-axis 
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

# Field Diagnostic Functions 
# * Make specified plots at location of simulation where diag_field() is called.    

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


# History diagnostics.  
# * Makes specified history plots from begining of simulation at point called.  
# * Many additional history diagnostics can be added by looking for
#   relevant moments accumulated in the Warp (see the variable group
#   "Hist" in top.v for an extensive list of variables that can be
#   used) and using gist commands to make relevant plots

def diag_hist(
 plt_ekin    = False, 
 plt_spnum   = False, 
 plt_curr_p  = False,
 plt_curr_e  = False,
 plt_lam_p   = False,
 plt_lam_e   = False,
 plt_lz      = False, 
 plt_pth     = False, 
 plt_pthn    = False,
 plt_krot    = False,
 plt_lang    = False,
 plt_cen     = False, 
 plt_envrms  = False, 
 plt_envmax  = False,
 plt_envrmsp = False, 
 plt_emit    = False, 
 plt_emitn   = False, 
 plt_emitg   = False,
 plt_emitng  = False,
 plt_emitr   = False,
 plt_emitnr  = False,
 plt_emitpv  = False, 
 plt_emitpvn = False, 
 plt_temp    = False,
 plt_Qperv   = False,
 plt_neutf   = False):
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
      plg(hl_ekin[0:top.jhist+1,js]/(A*kV),hl_zbeam[0:top.jhist+1],color=co)        
      #hpekin(js=js,color=co,titles=false,yscale=1./A,lhzbeam=true)    
    ptitles("History: Kinetic Energy","z [m]","KeV/u", )
    fma()
    # --- Operating species, in keV/u
    for ii in sort(sp_Operate.keys()):
      s = sp[ii]
      js = s.js
      co = s.color
      A  = s.mass/amu
      #hpekin(js=js,color=co,titles=false,yscale=1./A,lhzbeam=true)
      plg(hl_ekin[0:top.jhist+1,js]/(A*kV),hl_zbeam[0:top.jhist+1],color=co)
    ptitles("History: Operating Species Kinetic Energy","z [m]","KeV/u", )
    fma()
    # --- Support species, in keV/u 
    for ii in sort(sp_Support.keys()):
      s = sp[ii]
      js = s.js
      co = s.color
      A  = s.mass/amu
      plg(hl_ekin[0:top.jhist+1,js]/(A*kV),hl_zbeam[0:top.jhist+1],color=co)        
      #hpekin(js=js,color=co,titles=false,yscale=1./A,lhzbeam=true) # Was getting wrong answer !!
    ptitles("History: Support Species Kinetic Energy","z [m]","KeV/u", )
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
      plg(hl_ekin[0:top.jhist+1,js]/(Q*kV),hl_zbeam[0:top.jhist+1],color=co)        
      #hpekin(js=js,color=co,titles=false,yscale=1./Q,lhzbeam=true)
    plg(array([ekin_t,ekin_t]),array([zi,zf]),type="dash") 
    ptitles("History: "+lab+"Kinetic Energy","z [m]","KeV/Q", )
    limits(zi,zf,0.,1.2*ekin_t) 
    fma()
  # --- simulation particle number (to check for lost particles)
  #     Comment: tried using hppnum() but was unclear what was being plotted 
  if plt_spnum:
    # --- All Species Combined  
    plg(hl_spnumt[0:top.jhist+1],hl_zbeam[0:top.jhist+1])    
    ptitles("History: Live Sim Particle Number (all species)", "z [m]","Particle Number (simulation)", )
    fma()
    # --- All Species Individually 
    for ii in sort(sp.keys()):
      s = sp[ii]
      js = s.js
      co = s.color
      plg(hl_spnum[0:top.jhist+1,js],hl_zbeam[0:top.jhist+1],color=co)        
    ptitles("History: Live Sim Particle Number (by species)","z [m]","Particle Number (simulation)", )
    fma() 
    # --- Target Species
    lab = "" 
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg(hl_spnum[0:top.jhist+1,js],hl_zbeam[0:top.jhist+1],color=co)        
    ptitles("History: "+lab+" Live Sim Particle Number","z [m]","Particle Number (simulation)", )
    fma()             
  # --- current (particle)  
  if plt_curr_p:
    # --- All Species Combined  
    for ii in sort(sp.keys()):
      s = sp[ii]        
      js = s.js
      co = s.color
      plg(hl_ibeam_p[0:top.jhist+1,js]*1.e6,hl_zbeam[0:top.jhist+1],color=co)    
    ptitles("History: Species Particle Current (approx)", "z [m]","Current (microA)", )
    fma() 
    # --- Target Species 
    lab = "" 
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg(hl_ibeam_p[0:top.jhist+1,js]*1.e6,hl_zbeam[0:top.jhist+1],color=co)    
    ptitles("History: "+lab+" Particle Current (approx)","z [m]","Current (microA)", )
    fma() 
    # --- Total
    plg(hl_ibeam_pt[0:top.jhist+1]*1.e3,hl_zbeam[0:top.jhist+1])    
    ptitles("History: Total Particle Current (approx)","z [m]","Current (mA)", )
    fma()             
  # --- current (electrical)  
  if plt_curr_e:
    # --- All Species Combined  
    for ii in sort(sp.keys()):
      s = sp[ii]        
      js = s.js
      co = s.color
      plg(hl_ibeam_e[0:top.jhist+1,js]*1.e6,hl_zbeam[0:top.jhist+1],color=co)    
    ptitles("History: Species Electrical Current", "z [m]","Current (microA)", )
    fma()
    # --- Target Species
    lab = "" 
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg(hl_ibeam_e[0:top.jhist+1,js]*1.e6,hl_zbeam[0:top.jhist+1],color=co)    
    ptitles("History: "+lab+" Electrical Current","z [m]","Current (microA)", )
    fma()
    # --- Total
    plg(hl_ibeam_et[0:top.jhist+1]*1.e3,hl_zbeam[0:top.jhist+1])    
    ptitles("History: Total Electrical Current","z [m]","Current (mA)", )
    fma()                          
  # --- line charge (particle)  
  if plt_lam_p:
    # --- All Species Combined  
    for ii in sort(sp.keys()):
      s = sp[ii]        
      js = s.js
      co = s.color
      plg(hl_lambda_p[0:top.jhist+1,js]*10**9,hl_zbeam[0:top.jhist+1],color=co)    
    ptitles("History: Species Particle Line Charge", "z [m]","Line Charge (nC/m)", )
    fma()
    # --- Target Species
    lab = "" 
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg(hl_lambda_p[0:top.jhist+1,js]*10**9,hl_zbeam[0:top.jhist+1],color=co)    
    ptitles("History: "+lab+" Particle Line Charge","z [m]","Line Charge (nC/m)", )
    fma()             
  # --- line charge (electrical)  
  if plt_lam_e:
    # --- All Species Combined  
    for ii in sort(sp.keys()):
      s = sp[ii]        
      js = s.js
      co = s.color
      plg(hl_lambda_e[0:top.jhist+1,js]*10**9,hl_zbeam[0:top.jhist+1],color=co)    
    ptitles("History: Species Electrical Line Charge", "z [m]","Line Charge (nC/m)", )
    fma()
    # --- Target Species
    lab = "" 
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg(hl_lambda_e[0:top.jhist+1,js]*10**9,hl_zbeam[0:top.jhist+1],color=co)    
    ptitles("History: "+lab+" Electrical Line Charge","z [m]","Line Charge (nC/m)", )
    fma()             
  # --- lz mechanical angular momentum  
  if plt_lz:
    # --- All Species Combined  
    for ii in sort(sp.keys()):
      s = sp[ii]        
      js = s.js
      co = s.color
      plg(hl_lz[0:top.jhist+1,js]*10**6,hl_zbeam[0:top.jhist+1],color=co)    
    ptitles("History: Species Mechanical Angular Mom", "z [m]","<xy'>-<yx'>  [mm-mrad]", )
    fma()
    # --- Target Species
    lab = "" 
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg(hl_lz[0:top.jhist+1,js]*10**6,hl_zbeam[0:top.jhist+1],color=co)    
    ptitles("History: "+lab+" Mechanical Angular Mom","z [m]","<xy'>-<yx'>  [mm-mrad]", )
    fma()             
  # --- canonical angular momentum <P_theta>_j/(gamma_j*beta_j*m_j*c) in mm-mrad units   
  if plt_pth:
    # --- All Species Combined  
    for ii in sort(sp.keys()):
      s = sp[ii]        
      js = s.js
      co = s.color
      plg(hl_pth[0:top.jhist+1,js]*10**6,hl_zbeam[0:top.jhist+1],color=co)    
    ptitles("History: Species Canonical Angular Mom <Ptheta>/(gamma*beta*m*c)", "z [m]",
            "Canonical Ang Mom [mm-mrad]", )
    fma()
    # --- Target Species
    lab = "" 
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg(hl_pth[0:top.jhist+1,js]*10**6,hl_zbeam[0:top.jhist+1],color=co)    
    ptitles("History: "+lab+" Canonical Angular Mom <Ptheta>/(gamma*beta*m*c)","z [m]",
            "Canonical Ang Mom [mm-mrad]", )
    fma()             
  # --- canonical angular momentum (normalized) <P_theta>_j/(m_j*c) in mm-mrad units  
  if plt_pthn:
    # --- All Species Combined  
    for ii in sort(sp.keys()):
      s = sp[ii]        
      js = s.js
      co = s.color
      plg(hl_pthn[0:top.jhist+1,js]*10**6,hl_zbeam[0:top.jhist+1],color=co)    
    ptitles("History: Species Norm Canonical Angular Mom <Ptheta>/(m*c)", "z [m]",
            "Canonical Ang Mom [mm-mrad]", )
    fma()
    # --- Target Species
    lab = "" 
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg(hl_pthn[0:top.jhist+1,js]*10**6,hl_zbeam[0:top.jhist+1],color=co)    
    ptitles("History: "+lab+" Norm Canonical Angular Mom <Ptheta>/(m*c)","z [m]",
            "Canonical Ang Mom [mm-mrad]", )
    fma()             
  # --- effective rotation wavenumber 
  if plt_krot:
    # --- All Species Combined  
    for ii in sort(sp.keys()):
      s = sp[ii]        
      js = s.js
      co = s.color
      plg(hl_krot[0:top.jhist+1,js],hl_zbeam[0:top.jhist+1],color=co)    
    ptitles("History: Species Effective Rot Wavenumber", "z [m]","krot  [rad/m]", )
    fma()
    # --- Target Species
    lab = "" 
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg(hl_krot[0:top.jhist+1,js],hl_zbeam[0:top.jhist+1],color=co)    
    ptitles("History: "+lab+" Effective Rot Wavenumber","z [m]","krot  [rad/m]", )
    fma()             
  # --- Larmor rotation angle   
  if plt_lang:
    # --- All Species Combined  
    for ii in sort(sp.keys()):
      s = sp[ii]        
      js = s.js
      co = s.color
      plg((180./pi)*hl_lang[0:top.jhist+1,js],hl_zbeam[0:top.jhist+1],color=co)    
    ptitles("History: Species Larmor Rot Angle", "z [m]","Rotation [deg]", )
    fma()
    # --- Target Species
    lab = "" 
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg((180./pi)*hl_lang[0:top.jhist+1,js],hl_zbeam[0:top.jhist+1],color=co)    
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
    # --- All Species Combined, x- and y-plane: Factor 4 in scale to account for Warp edge measure 
    hpepsx(titles=false,yscale=1./(4.*mm*mr),lhzbeam=true)
    hpepsy(titles=false,yscale=1./(4.*mm*mr),lhzbeam=true,color="red")
    ptitles("History: All Species RMS x-, y-Emittance: x[b],y[r]","z [m]","Emittance [mm-mr]", )
    fma()
    # --- Target Species, x-plane: Factor 4 in scale to account for Warp edge measure  
    hpepsx(titles=false,yscale=1./(4.*mm*mr),lhzbeam=true)
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsx(js=js,color=co,titles=false,yscale=1./(4.*mm*mr),lhzbeam=true)    
    ptitles("History: "+lab+"RMS x-Emittance","z [m]","Emittance [mm-mr]", )
    fma()
    #
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsx(js=js,color=co,titles=false,yscale=1./(mm*mr),lhzbeam=true)    
    ptitles("History: "+lab+"RMS x-Emittance","z [m]","Emittance [mm-mr]", )
    fma()
    # --- Target Species, y-plane 
    hpepsy(titles=false,yscale=1./(4.*mm*mr),lhzbeam=true)
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsy(js=js,color=co,titles=false,yscale=1./(4.*mm*mr),lhzbeam=true)    
    ptitles("History: "+lab+"RMS y-Emittance","z [m]","Emittance [mm-mr]", )
    fma()
    #
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsy(js=js,color=co,titles=false,yscale=1./(mm*mr),lhzbeam=true)    
    ptitles("History: "+lab+"RMS y-Emittance","z [m]","Emittance [mm-mr]", )
    fma()
  # --- emittance, normalized 
  if plt_emitn:
    # --- All Species Combined, x- and y-plane 
    #     ** warning norm emittance scaled mm-mrad by default in Warp **
    hpepsnx(titles=false,yscale=1./4.,lhzbeam=true)
    hpepsny(titles=false,yscale=1./4.,lhzbeam=true,color="red")
    ptitles("History: All Species Norm RMS x-, y-Emittance: x[b],y[r]","z [m]","Norm Emittance [mm-mr]", )
    fma()
    # --- By Target Species, x-plane 
    hpepsnx(titles=false,yscale=1./4.,lhzbeam=true)
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsnx(js=js,color=co,titles=false,yscale=1./4.,lhzbeam=true)    
    ptitles("History: "+lab+"Norm RMS x-Emittance","z [m]","Norm Emittance [mm-mr]", )
    fma()
    #
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsnx(js=js,color=co,titles=false,yscale=1./4.,lhzbeam=true)    
    ptitles("History: "+lab+"Norm RMS x-Emittance","z [m]","Norm Emittance [mm-mr]", )
    fma()
    # --- By Target Species, y-plane 
    hpepsny(titles=false,yscale=1./4.,lhzbeam=true)
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsny(js=js,color=co,titles=false,yscale=1./4.,lhzbeam=true)    
    ptitles("History: "+lab+"Norm RMS y-Emittance","z [m]","Norm Emittance [mm-mr]", )
    fma()
    #
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsny(js=js,color=co,titles=false,yscale=1./4.,lhzbeam=true)    
    ptitles("History: "+lab+"Norm RMS y-Emittance","z [m]","Emittance [mm-mr]", )
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
    ptitles("History: "+lab+"RMS g-Emittance","z [m]","Emittance [mm-mr]", )
    fma()
    #
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsg(js=js,color=co,titles=false,yscale=1./(mm*mr),lhzbeam=true)    
    ptitles("History: "+lab+"RMS g-Emittance","z [m]","Emittance [mm-mr]", )
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
  #     ** scaled mm-mrad by defualt in Warp ** 
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
    ptitles("History: "+lab+"RMS Norm g-Emittance","z [m]","Norm Emittance [mm-mr]", )
    fma()
    #
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsng(js=js,color=co,titles=false,yscale=1.,lhzbeam=true)    
    ptitles("History: "+lab+"RMS Norm g-Emittance","z [m]","Norm Emittance [mm-mr]", )
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
    hpepsr(titles=false,yscale=1./(2.*mm*mr),lhzbeam=true)
    ptitles("History: All Species RMS r-Emittance","z [m]","Emittance [mm-mr]", )
    fma()
    # --- By Target Species  
    hpepsr(titles=false,yscale=1./(2.*mm*mr),lhzbeam=true)
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsr(js=js,color=co,titles=false,yscale=1./(2.*mm*mr),lhzbeam=true)    
    ptitles("History: "+lab+"RMS r-Emittance","z [m]","Emittance [mm-mr]", )
    fma()
    #
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsr(js=js,color=co,titles=false,yscale=1./(2.*mm*mr),lhzbeam=true)    
    ptitles("History: "+lab+"RMS r-Emittance","z [m]","Emittance [mm-mr]", )
    fma()
  # --- emittance, generalized radial normalized ** warning norm emittance scaled mm-mrad by default **
  if plt_emitnr:
    # --- All Species Combined
    hpepsnr(titles=false,yscale=1./2.,lhzbeam=true)
    ptitles("History: All Species Norm RMS r-Emittance","z [m]","Norm Emittance [mm-mr]", )
    fma()
    # --- By Target Species  
    hpepsnr(titles=false,yscale=1./2.,lhzbeam=true)
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsnr(js=js,color=co,titles=false,yscale=1./2.,lhzbeam=true)    
    ptitles("History: "+lab+"RMS Norm r-Emittance","z [m]","Norm Emittance [mm-mr]", )
    fma()
    #
    lab = ""    
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpepsnr(js=js,color=co,titles=false,yscale=1./2.,lhzbeam=true)    
    ptitles("History: "+lab+"RMS Norm r-Emittance","z [m]","Norm Emittance [mm-mr]", )
    fma()
  # --- emittance, total phase volume, unnormalized 
  if plt_emitpv:
    # --- All Species Combined  
    for ii in sort(sp.keys()):
      s = sp[ii]        
      js = s.js
      co = s.color
      plg(hl_epspv[0:top.jhist+1,js]/(mm*mr),hl_zbeam[0:top.jhist+1],color=co)    
    ptitles("History: Species Total Phase Volume Emittance", "z [m]","Emittance [mm-mrad]", )
    fma()
    # --- Target Species
    lab = "" 
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg(hl_epspv[0:top.jhist+1,js]/(mm*mr),hl_zbeam[0:top.jhist+1],color=co)    
    ptitles("History: "+lab+" Total Phase Volume Emittance","z [m]","Emittance [mm-mrad]", )
    fma()             
  # --- emittance, total phase volume, normalized  
  if plt_emitpvn:
    # --- All Species Combined  
    for ii in sort(sp.keys()):
      s = sp[ii]        
      js = s.js
      co = s.color
      plg(hl_epspvn[0:top.jhist+1,js]/(mm*mr),hl_zbeam[0:top.jhist+1],color=co)    
    ptitles("History: Species Total Phase Volume Norm Emittance", "z [m]","Norm Emittance [mm-mrad]", )
    fma()
    # --- Target Species
    lab = "" 
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg(hl_epspvn[0:top.jhist+1,js]/(mm*mr),hl_zbeam[0:top.jhist+1],color=co)    
    ptitles("History: "+lab+" Total Phase Volume Norm Emittance","z [m]","Norm Emittance [mm-mrad]", )
    fma()             
  # --- Effective ion temperature calculated from radial thermal emittance  
  if plt_temp:
    # --- All Species Combined  
    for ii in sort(sp.keys()):
      s = sp[ii]        
      js = s.js
      co = s.color
      plg(hl_temp[0:top.jhist+1,js],hl_zbeam[0:top.jhist+1],color=co)    
    ptitles("History: Species Transverse Thermal Temperature", "z [m]","Temp [eV]", )
    fma()
    # --- Target Species
    lab = "" 
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg(hl_temp[0:top.jhist+1,js],hl_zbeam[0:top.jhist+1],color=co)    
    ptitles("History: "+lab+" Species Transverse Thermal Temperature","z [m]","Temp [eV]", )
    fma()             
  # --- Perveance   
  if plt_Qperv:
    # --- All Species Combined  
    for ii in sort(sp.keys()):
      s = sp[ii]        
      js = s.js
      co = s.color
      plg(hl_Qperv[0:top.jhist+1,js],hl_zbeam[0:top.jhist+1],color=co)    
    ptitles("History: Species Bare Perveance Q", "z [m]","Perveance [1]", )
    fma()
    # --- Target Species
    lab = "" 
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg(hl_Qperv[0:top.jhist+1,js],hl_zbeam[0:top.jhist+1],color=co)    
    ptitles("History: "+lab+" Species Bare Perveance Q","z [m]","Perveance [1]", )
    fma()             
  # --- Neutralization Factor 
  if plt_neutf:
    # --- All Species Combined  
    for ii in sort(sp.keys()):
      s = sp[ii]        
      js = s.js
      co = s.color
      plg(hl_neutf[0:top.jhist+1,js],hl_zbeam[0:top.jhist+1],color=co)    
    ptitles("History: Species Electron Neutralization Fractions", "z [m]","Fraction [1]", )
    fma()
    # --- Target Species
    lab = "" 
    for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      plg(hl_neutf[0:top.jhist+1,js],hl_zbeam[0:top.jhist+1],color=co)    
    ptitles("History: "+lab+" Electron Neutralization Factors","z [m]","Fraction [1]", )
    fma()               


#  -- Install diagnostics at appropriate intervals after steps
#       Add options to generate plots desired 
#  -- Install diagnostics at appropriate intervals after steps
#       Add options to generate plots desired 

# Function to call diagnostics at a timestep in step control lists 
def diag_calls():
  if top.it in diag_part_step:
    diag_part(plt_xy=true,plt_xxp=true,plt_yyp=false,plt_xpyp=true,
              plt_trace=false,plt_denxy=true,plt_denr=true)
  if top.it in diag_field_step: 
    diag_field(plt_pc=true,plt_pc_xy=true,plt_pa=true)
  if top.it in diag_hist_step:
    diag_hist(plt_ekin=true,plt_spnum=true,plt_curr_e=true,plt_curr_p=true,plt_lam_p=true,plt_lam_e=true,
              plt_lz=true,plt_pth=true,plt_pthn=true,plt_krot=true,plt_lang=true, 
              plt_cen=true,plt_envrms=true,plt_envmax=true,plt_envrmsp=true,  
              plt_emit=true,plt_emitn=true,plt_emitg=true,plt_emitng=true,plt_emitr=true,plt_emitnr=true, 
              plt_emitpv=true,plt_emitpvn=true,plt_temp=true,plt_Qperv=true,plt_neutf=true)



