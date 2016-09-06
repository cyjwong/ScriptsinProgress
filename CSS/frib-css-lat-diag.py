# Diagnostics for applied field lattice  
# FRIB Front End 

# Input Parameters 
# Setup mesh line for extraction of applied field components  
#  dl_ .... suffix for "diagnostic lattice"

dl_nz     = 1000    # number z-steps for x-y field plots        
dl_z_mmin = 65.5    # min z to make plots 
dl_z_mmax = 69.     # max z to make plots 
dl_r      = 4.*mm   # radius to use for some field plots 

dl_sym_check = True # make symmetry check plots of fields to check lattice descriptions  

# Script calculated parameters 
 
# --- mesh line 
dl_dz = (dl_z_mmax-dl_z_mmin)/dl_nz
dl_z_mesh = dl_z_mmin + arange(dl_nz+1)*dl_dz  

# --- field variables 
dl_ex_mesh = zeros(dl_nz+1)
dl_ey_mesh = zeros(dl_nz+1)
dl_ez_mesh = zeros(dl_nz+1)

dl_bx_mesh = zeros(dl_nz+1)
dl_by_mesh = zeros(dl_nz+1)
dl_bz_mesh = zeros(dl_nz+1)


# Function to extract applied E and B fields
#   * Returns fields on x=r,y=0,z=dl_zmesh line 
#   * Fields will be linear or nonlinear depending how setup in lattice 
def dl_getappliedfields(r): 
  x_mesh = r*ones(dl_nz+1) 
  y_mesh = zeros(dl_nz+1) 
  (ex,ey,ez,bx,by,bz) = getappliedfields(x=x_mesh,y=y_mesh,z=dl_z_mesh)
  return(ex,ey,ez,bx,by,bz)

## Nonlinear S4 Solenoid using Warp built-in script to plot gridded field elements 
##  * Plots over extent of mesh definition by default so should be ranged appropriately to see structure  
##  * ONLY USED if nonlinear description of elements set 

#if s4p1_typ == "nl":

	#plotbgrd(ib=s4_nl_id,component='z',ix=0,iy=0,titles=false)
	#ptitles("NL S4 Solenoid Data: Scaled on-axis Bz vs z","z [m]","Bz [1]",)
	#fma()
	
	#plotbgrd(ib=s4_nl_id,component='x',ix=s4_nr/2,iy=0,titles=false)
	#ptitles("NL S4 Solenoid Data: Scaled Br vs z at r = %s mm"%(s4_r_m[s4_nr/2]/mm),"z [m]","Br [1]",)
	#fma()

## Nonlinear Grated Accel Gap using Warp built-in script to plot gridded field elements 
##  * Plots over extent of mesh definition by default so should be ranged appropriately to see structure  
##  * ONLY USED if nonlinear description of elements set 

#if gag_typ == "nl":

	#plotegrd(ie=gag[0],component='z',ix=0,iy=0,titles=false,color=green)
	#ptitles("NL Grated Accel Gap Data: on-axis Ez vs z","z [m]","Ez [V/m]",)
	#fma()
	
	#plotegrd(ie=gag[0],component='x',ix=gag_nr/2,iy=0,titles=false,color=green)
	#ptitles("NL Grated Accel Gap Data: Er vs z at r = %s mm"%(gag_r_m[gag_nr/2]/mm),"z [m]","Er [V/m]",)
	#fma()

# Calculate applied fields on r=0 axis 
(dl_ex_mesh,dl_ey_mesh,dl_ez_mesh,dl_bx_mesh,dl_bymesh,dl_bz_mesh) = dl_getappliedfields(0.) 

# Axial Magnetic Field on r=0 axis 
plg(dl_bz_mesh,dl_z_mesh)
ptitles("Axial Magnetic Field at r=0","z [m]","Bz [Tesla]",)
fma() 

# Axial Electric Field on r=0 axis 
plg(dl_ez_mesh,dl_z_mesh,color=green) 
ptitles("Axial Electric Field at r=0","z [m]","Ez [V/m]") 
fma() 

# Calculate applied fields on r=dl_r axis 
(dl_ex_mesh,dl_ey_mesh,dl_ez_mesh,dl_bx_mesh,dl_bymesh,dl_bz_mesh) = dl_getappliedfields(dl_r) 

# Radial Magnetic Field on r=dl_r axis 
plg(dl_bx_mesh,dl_z_mesh)
ptitles("Radial Magnetic Field at r=%s mm"%(dl_r/mm),"z [m]","Er [V/m]",)
fma() 

# Radial Electric Field on r=dl_r axis 
plg(dl_ex_mesh,dl_z_mesh,color=green) 
ptitles("Radial Electric Field at r=%s mm"%(dl_r/mm),"z [m]","Er [V/m]",)
fma() 

# Symmetry check of Br and Er that signs are correct along principal axes  
#  * plot field at x = +- dl_rf to make sure that field structure has mirror image 
#  * program inefficient: not used much so just lay out easy to read 

if dl_sym_check:
  dl_y_mesh = zeros(dl_nz+1) 
  dl_x_mesh = dl_r*ones(dl_nz+1) 
  (ex_meshp,ey_meshp,ez_meshp,bx_meshp,by_meshp,bz_meshp) =\
    getappliedfields(x= dl_x_mesh,y=dl_y_mesh,z=dl_z_mesh)  
  dl_x_mesh = -dl_r*ones(dl_nz+1) 
  (ex_meshm,ey_meshm,ez_meshm,bx_meshm,by_meshm,bz_meshm) =\
    getappliedfields(x= dl_x_mesh,y=dl_y_mesh,z=dl_z_mesh)  
  # Solenoids: x-plane    
  plg(bx_meshp,dl_z_mesh,color="black")
  plg(bx_meshm,dl_z_mesh,color="red"  )
  ptitles("Bx Magnetic Field at y=0, x = +-%s mm B/R"%(dl_r/mm),"z [m]","Bx [Tesla]",)
  fma()
  # Grated Gap: x-plane  
  plg(ex_meshp,dl_z_mesh,color="black")
  plg(ex_meshm,dl_z_mesh,color="red"  )
  ptitles("Ex Electric Field at y=0, x = +-%s mm B/R"%(dl_r/mm),"z [m]","Ex [V/m]",)
  fma()
  # 
  dl_x_mesh = zeros(dl_nz+1)
  dl_y_mesh = dl_r*ones(dl_nz+1) 
  (ex_meshp,ey_meshp,ez_meshp,bx_meshp,by_meshp,bz_meshp) =\
    getappliedfields(x= dl_x_mesh,y=dl_y_mesh,z=dl_z_mesh)  
  dl_y_mesh = -dl_r*ones(dl_nz+1) 
  (ex_meshm,ey_meshm,ez_meshm,bx_meshm,by_meshm,bz_meshm) =\
    getappliedfields(x= dl_x_mesh,y=dl_y_mesh,z=dl_z_mesh)  
  # Solenoids: y-plane 
  plg(by_meshp,dl_z_mesh,color="black")
  plg(by_meshm,dl_z_mesh,color="red"  )
  ptitles("By Magnetic Field at x=0, y = +-%s mm B/R"%(dl_r/mm),"z [m]","By [Tesla]",)
  fma()
  # Grated Gap: y-plane 
  plg(ey_meshp,dl_z_mesh,color="black")
  plg(ey_meshm,dl_z_mesh,color="red"  )
  ptitles("Ey Electric Field at x=0, y = +-%s mm B/R"%(dl_r/mm),"z [m]","Ey [V/m]",)
  fma()

## On axis Bz element by element 
##  * Plots appropriate for linear or nonlinear elements 
##  * Scale should range on the coordinates defined to better visualize field structure 

#def diag_field_ele_Bz(type):
  #dict_id  = {'s4p1':s4p1[0],  's4p2':s4p2[0], 'ecr':ecr[0]}
  #dict_typ = {'s4p1':s4p1_typ, 's4p2':s4p2_typ,'ecr':ecr_typ}
  #id  = dict_id[type]
  #typ = dict_typ[type]
  #if typ == "lin":
    #plotmmlt(im=id,titles=false)
  #elif typ == "nl":
    #plotbgrd(ib=id,component='z',ix=0,iy=0,titles=false)
  #else:
    #return

#diag_field_ele_Bz('ecr')
#ptitles("On Axis Bz of ECR","z [m]","Bz [Tesla]",)
#fma() 

#diag_field_ele_Bz('s4p1')
#ptitles("On Axis Bz of 1st S4","z [m]","Bz [Tesla]",)
#fma() 

#diag_field_ele_Bz('s4p2')
#ptitles("On Axis Bz of 2nd S4","z [m]","Bz [Tesla]",)
#fma()

#diag_field_ele_Bz('ecr')
#diag_field_ele_Bz('s4p1')
#diag_field_ele_Bz('s4p2')
#ptitles("On Axis Bz of Elements","z [m]","Bz [Tesla]",)
#fma() 


# Plot of Bz and Ez on Axis on a similar scale to help understand field structure
#  * Take over range of simulation 

dlr_nz = 1000  
dlr_z_mesh = linspace(ecr_z_extr,z_adv,dlr_nz+1)

dlr_x_mesh = zeros(dlr_nz+1)
dlr_y_mesh = zeros(dlr_nz+1) 

(dlr_ex_mesh,dlr_ey_mesh,dlr_ez_mesh,dlr_bx_mesh,dlr_by_mesh,dlr_bz_mesh) =\
  getappliedfields(x=dlr_x_mesh,y= dlr_y_mesh,z=dlr_z_mesh)

dlr_ez_max = maxnd( abs(dlr_ez_mesh) ) 
dlr_bz_max = maxnd( abs(dlr_bz_mesh) ) 

plg([1.,0.],[z_launch,z_launch], type = "dash")
plg(dlr_bz_mesh/dlr_bz_max, dlr_z_mesh)
plg(dlr_ez_mesh/dlr_ez_max, dlr_z_mesh, color = "green")  
ptitles("Scaled Axial Field of Elements (B: Bz) and (G: Ez)","z [m]","Bz and Ez",) 
fma() 

# Repeat Plot of Bz and Ez on Axis on a similar scale to help understand field structure
# with added range over ECR and launch location indicated.  
#  * Take over range of simulation 

dlr_z_mesh = linspace(ecr_z_extr - 1.1,z_adv,dlr_nz+1)

(dlr_ex_mesh,dlr_ey_mesh,dlr_ez_mesh,dlr_bx_mesh,dlr_by_mesh,dlr_bz_mesh) =\
  getappliedfields(x=dlr_x_mesh,y= dlr_y_mesh,z=dlr_z_mesh)

dlr_ez_max = maxnd( abs(dlr_ez_mesh) ) 
dlr_bz_max = maxnd( abs(dlr_bz_mesh) ) 

plg([1.,0.],[z_launch,z_launch], type = "dash")
plg(dlr_bz_mesh/dlr_bz_max, dlr_z_mesh)
plg(dlr_ez_mesh/dlr_ez_max, dlr_z_mesh, color = "green")  
ptitles("Scaled Axial Field of Elements (B: Bz) and (G: Ez)","z [m]","Bz and Ez",) 
fma() 

# Combined plot of Bz and Ez with launch/ECR extractor point indicated 
diag_field_ele_Bz('ecr')
diag_field_ele_Bz('s4p1')
diag_field_ele_Bz('s4p2')
plg(2.*dlr_ez_mesh/dlr_ez_max, dlr_z_mesh, color = "green")  
plg([3.5,0.],[z_launch,z_launch], type = "dash")
limits('e','e',0.,3.5)
ptitles("On Axis Bz (B) of Elements and Scaled Ez (G)","z [m]","Bz [Tesla]",)
fma() 
  
