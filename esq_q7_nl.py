"""
Script to generate electric field data for Q7 type ESQ in WARP by reading in data from an ascii text file provided by CST Studio  

Nonlinear 3D E field on Cartesian mesh

Steve Lund 
lund@frib.msu.edu  
"""

from warp import *

setup() 


# Data from mechanical drawing 
len_coil   = 29.59*cm 
len_magnet = 39.51*cm 
r_coil_i   = 7.75*cm
r_coil_o   = 22.5*cm  

aspect = r_coil_i/len_coil 





electrode_voltage = 10000 # [V]

# Read in data file from CST Studio
#   Data has E_x(x,y,z),E_y(x,y,z) and E_z(x,y,z) from on a uniform Cartesian mesh.  
#   Use getdatafromtextfile() Warp function.
#    d_ prefix for data,  _m suffix for mesh 

date      = "20160607"
data_file = "Q7_mesh2mm_eighth.txt"

[d_x,d_y,d_z,d_ex,d_ey,d_ez] = getdatafromtextfile(data_file,nskip=2,dims=[6,None],)


# Plot of raw data 
#plg(d_bz0_m,d_z_m)
#limits('e','e','e',0.)
#ptitles('S4 Solenoid Raw Data: On-Axis B_z_','z [cm]','B_z_(r=0,z) [Gauss]',) 
#fma() 

# Scale data for z_m in meters and bz in fraction of peak on-axis field at center  

d_x  = d_x*mm 
d_y  = d_y*mm
d_z  = d_z*mm 




#def ext_field(x):
  #x1 = x.max()
  #a1 = abs(x1)  
  #x2 = x.min()
  #a2 = abs(x2) 
  #if a1 >= a2:
    #return(x1) 
  #else: 
    #return(x2) 

# d_bz_peak = ext_field(d_bz)   Danger: peak field was in outer radius!  Nonlinear strength increases, pick off on-axis value 
# d_br = d_br/d_bz_peak 
# d_bz = d_bz/d_bz_peak



d_ex /= electrode_voltage   # Normalize against electrode voltage
d_ey /= electrode_voltage
d_ez /= electrode_voltage


# Calculate dx,dy,dz for uniform mesh data   

dx = average(diff(unique(d_x)))
dy = average(diff(unique(d_y)))
dz = average(diff(unique(d_z)))

## Generate data mesh and field arrays with +z and -z (reflected) coordinates 
## excluding boundary points.  df_ prefix for "data full".
##    Note:  z coordinates flipped in sign since d_z contains data from z = 0 to z = max (downstream) 
##           br data flipped in sign consistently 
     
#df_z   = concatenate(( -d_z[-1:0:-1], d_z[:])) 
#df_r   = concatenate((  d_r[-1:0:-1], d_r[:])) 
#df_br  = concatenate((-d_br[-1:0:-1],d_br[:]))
#df_bz  = concatenate(( d_bz[-1:0:-1],d_bz[:]))




# Generate data mesh and field arrays from reduced data set that exploits symmetry
# Data set only covers one quadrant of the transverse plane and one half of the length, i.e. 1/8 of the total volume
# Duplication of mesh points in the array apparently does not matter

## Extend data set with x -> -x, E_x -> -E_x

#d_x   = concatenate(( -d_x[-1:0:-1], d_x[:])) 
#d_y   = concatenate(( d_y[-1:0:-1], d_y[:])) 
#d_z   = concatenate(( d_z[-1:0:-1], d_z[:])) 
#d_ex  = concatenate(( -d_ex[-1:0:-1],d_ex[:]))
#d_ey  = concatenate(( d_ey[-1:0:-1],d_ey[:]))
#d_ez  = concatenate(( d_ez[-1:0:-1],d_ez[:]))


## Extend data set with y -> -y, E_y -> -E_y

#d_x   = concatenate(( d_x[-1:0:-1], d_x[:])) 
#d_y   = concatenate(( -d_y[-1:0:-1], d_y[:])) 
#d_z   = concatenate(( d_z[-1:0:-1], d_z[:])) 
#d_ex  = concatenate(( d_ex[-1:0:-1],d_ex[:]))
#d_ey  = concatenate(( -d_ey[-1:0:-1],d_ey[:]))
#d_ez  = concatenate(( d_ez[-1:0:-1],d_ez[:]))


## Extend data set with z -> -z, E_z -> -E_z

#d_x   = concatenate(( d_x[-1:0:-1], d_x[:])) 
#d_y   = concatenate(( d_y[-1:0:-1], d_y[:])) 
#d_z   = concatenate(( -d_z[-1:0:-1], d_z[:])) 
#d_ex  = concatenate(( d_ex[-1:0:-1],d_ex[:]))
#d_ey  = concatenate(( d_ey[-1:0:-1],d_ey[:]))
#d_ez  = concatenate(( -d_ez[-1:0:-1],d_ez[:]))





# Avoiding duplication:

#num_reduced_x_grid = nint((d_x.[-1]-d_x.[0])/dx) + 1
#num_reduced_y_grid = nint((d_y.[-1]-d_y.[0])/dy) + 1
#num_reduced_z_grid = nint((d_z.[-1]-d_z.[0])/dz) + 1

#assert num_reduced_x_grid * num_reduced_y_grid * num_reduced_z_grid == len(d_x)

#k = num_reduced_y_grid*num_reduced_z_grid

df_x = d_x.tolist()[:]
df_y = d_y.tolist()[:]
df_z = d_z.tolist()[:]
df_ex = d_ex.tolist()[:]
df_ey = d_ey.tolist()[:]
df_ez = d_ez.tolist()[:]


for i in range(len(df_x)):
	if df_x[i] != 0:
		df_x.append(-df_x[i])
        df_y.append(df_y[i])
        df_z.append(df_z[i])
        df_ex.append(-df_ex[i])
        df_ey.append(df_ey[i])
        df_ez.append(df_ez[i])

for i in range(len(df_y)):
	if df_y[i] != 0:
		df_x.append(df_x[i])
        df_y.append(-df_y[i])
        df_z.append(df_z[i])
        df_ex.append(df_ex[i])
        df_ey.append(-df_ey[i])
        df_ez.append(df_ez[i])

for i in range(len(df_z)):
	if df_z[i] != 0:
		df_x.append(df_x[i])
        df_y.append(df_y[i])
        df_z.append(-df_z[i])
        df_ex.append(df_ex[i])
        df_ey.append(df_ey[i])
        df_ez.append(-df_ez[i])



df_x = asarray(df_x)
df_y = asarray(df_y)
df_z = asarray(df_z)
df_ex = asarray(df_ex)
df_ey = asarray(df_ey)
df_ez = asarray(df_ez)


# Calculate mesh bounds and increments and radial/axial mesh vectors 
#   r_m = radial mesh vector 
#   z_m = axial  mesh vector 

nx = nint((df_x.max()-df_x.min())/dx)
x_mmin = df_x.min()
x_mmax = df_x.max()

ny = nint((df_y.max()-df_y.min())/dy)
y_mmin = df_y.min()
y_mmax = df_y.max()

nz = nint((df_z.max()-df_z.min())/dz)
z_mmin = df_z.min()
z_mmax = df_z.max()


x_m = x_mmin + dx*arange(0,nx+1,dtype=float) 
y_m = y_mmin + dy*arange(0,ny+1,dtype=float) 
z_m = z_mmin + dz*arange(0,nz+1,dtype=float) 

# Load field data in arrays 
ex_m = fzeros((nx+1,ny+1,nz+1)) 
ex_m = fzeros((nx+1,ny+1,nz+1)) 
ez_m = fzeros((nx+1,ny+1,nz+1)) 

ix_m = nint((df_x - x_mmin)/dx)
iy_m = nint((df_y - y_mmin)/dy)
iz_m = nint((df_z - z_mmin)/dz)

ii = ix_m + iy_m*(nx+1) + iz_m*(nx+1)*(ny+1)

ex_m.ravel(order='F').put(ii,df_ex)
ey_m.ravel(order='F').put(ii,df_ey)
ez_m.ravel(order='F').put(ii,df_ez)


#ir_half = sum(where(r_m < 0.5*r_coil_i, 1, 0)) 
#ir_90   = sum(where(r_m < 0.9*r_coil_i, 1, 0)) 

## Plots of processed data 
## --- B_z(r,z) vs z at a few values of r  
#plg(bz_m[0,      :],z_m/mm,color='black')
#plg(bz_m[ir_half,:],z_m/mm,color='blue' ) 
#plg(bz_m[ir_90,  :],z_m/mm,color='red'  ) 
#plg([1.1,-0.1],[-(len_coil/2.)/mm,-(len_coil/2.)/mm],type='dash',color='green')
#plg([1.1,-0.1],[ (len_coil/2.)/mm, (len_coil/2.)/mm],type='dash',color='green')
#plg([1.1,-0.1],[-(len_magnet/2.)/mm,-(len_magnet/2.)/mm],type='dash',color='magenta')
#plg([1.1,-0.1],[ (len_magnet/2.)/mm, (len_magnet/2.)/mm],type='dash',color='magenta')
#limits('e','e',-0.1,1.1)
#ptitles('S4 Bz: r=0,rc/2,09*rc Black,Blue,Red','z [mm]','Bz(r,z)/Bz(0,0)',) 
#fma() 

## --- Linear optic B_r vs z 
#ymax = 0.6
#plg(br_m[0,      :],z_m/mm,color='black')
#plg(br_m[ir_half,:],z_m/mm,color='blue' )
#plg(br_m[ir_90,  :],z_m/mm,color='red'  )
#plg([ymax,-ymax],[-(len_coil/2.)/mm,-(len_coil/2.)/mm],type='dash',color='green')
#plg([ymax,-ymax],[ (len_coil/2.)/mm, (len_coil/2.)/mm],type='dash',color='green')
#plg([ymax,-ymax],[-(len_magnet/2.)/mm,-(len_magnet/2.)/mm],type='dash',color='magenta')
#plg([ymax,-ymax],[ (len_magnet/2.)/mm, (len_magnet/2.)/mm],type='dash',color='magenta')
#limits('e','e',-ymax,ymax)
#ptitles('S4 Br: r=0,rc/2,0.9*rc Black,Blue,Red','z [mm]','Br(r,z)/Bz(0,0)',) 
#fma() 

# Save linear field data to an external binary array 
fo = PWpickle.PW("q7.3d."+date+".pkl")

#fo.s4_len_coil   = len_coil
#fo.s4_len_magnet = len_magnet
#fo.s4_r_coil_i   = r_coil_i
#fo.s4_r_coil_o   = r_coil_o

fo.q7_dx = dx
fo.q7_dy = dy
fo.q7_dz = dz
fo.q7_nx = nx   
fo.q7_ny = ny   
fo.q7_nz = nz 
fo.q7_x_m     = x_m 
fo.q7_y_m     = y_m 
fo.q7_z_m     = z_m
fo.q7_z_m_cen = 0.      # center position on z_m   
fo.q7_ex_m    = ex_m
fo.q7_ey_m    = ey_m
fo.q7_ez_m    = ez_m


fo.close() 

