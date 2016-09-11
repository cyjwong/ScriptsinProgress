# hl_pthn, hl_pth, hl_esprn are redefined in frib-scan-xy-params-U.py


# Axisymmetric envelope model solution.
#   * Designed to use in a variety of simulations with inputs set accordingly. 
#   * Setup below for use in front end simulation script frib-front-xy.py
#   * Separate diagnostic scripts in frib-front-env-diag.py used to make 
#     diagnostic plots of the solution.   

#########################
# Begin Inputs 
#########################

# Range to advance envelope and increment of advance 
#  env_zs = start envelope advance [m]
#  env_ze = end   envelope advance [m]
#  env_ds = envelope solution step size [m]

env_zs = z_launch
env_ze = z_adv 		# (halt the envelope solver before the first bend)
env_ds = 5.*mm  

# Longitudinal velocity correction terms in envelope model 
# 
#  CorrectionMode = 0 - On-axis E_z only  
#                   1 - On-axis E_z + dBdz term from magnetic field  
#                   2 - On-axis E_z + dBdz term from magnetic field + d2Edz2 term from electric field
#

CorrectionMode = 1 


# Solve Envelope Model using Warp simulation data 
#   integratewarp = True  : Use pic data for all z in integration range for energy, ps areas 
#                 = False : Use pic data at z_begin only  
#   * Case 1 also generates case 0 solution 

integratewarp = False


# Numerical Integration Control
#
#  mxstep = number of steps to take in integration between env_zs and env_ze 
#  rtol   = relative numerical tolerance 
#  atol   = absolute numerical tolerance 

mxstep = 5000
rtol = 1.49012e-8
atol = 1.49012e-8


#########################
# End Inputs 
#########################

# Import ode solver 
from scipy.integrate import odeint

# Make time array for solution based on advance range and step size

stepnum = int(round((env_ze - env_zs)/(1.*cm)))

sss = linspace(env_zs, env_ze, stepnum)

stepsize = (env_ze - env_zs)/(stepnum - 1)


# Data needed in Env. Model

speciesq    = zeros(top.ns)
speciesI    = zeros(top.ns)
specieslist = [0]*top.ns 

for ii in sp.keys():
  s = sp[ii]
  js = s.js
  speciesq[js] = s.charge
  speciesI[js] = ibeam[ii]
  specieslist[js] = sp[ii]


# Form state vector to advance
#   Length 3*ns where ns is the total of species
# 
#   1st ns block: KE of species
#   2nd ns block: sigma_r  of species 
#   3rd ns block: d(sigma_r)/dz of species
# 
# * Same ordering used in rest of code 
 
#state_vector = [0]*3*top.ns
state_vector = zeros(3*top.ns) 

deltaz = stepsize/2.

termdata = []


## Set up the derivative of the state vector

def f(state_vector, rrr):
	
	efieldz = getappliedfields(0, 0, rrr)[2][0]
	bfieldz = getappliedfields(0, 0, rrr)[5][0]
	
	dEdz = (getappliedfields(0, 0, rrr + deltaz/2)[2][0] - getappliedfields(0, 0, rrr - deltaz/2)[2][0])/deltaz
	dBdz = (getappliedfields(0, 0, rrr + deltaz/2)[5][0] - getappliedfields(0, 0, rrr - deltaz/2)[5][0])/deltaz
	
	d2Edz2 = (getappliedfields(0, 0, rrr + deltaz/2)[2][0] - getappliedfields(0, 0, rrr)[2][0]*2 + getappliedfields(0, 0, rrr - deltaz/2)[2][0])/(deltaz**2)*4
	#dVdz = 0
	
	derivs = []	
	
	speciesbeta = []
	
	for i in range(top.ns):
		speciesbeta.append(sqrt((2*state_vector[i]*jperev)/(specieslist[i].mass*clight**2)))
		
	
## build first lot in deriv output (i.e. dKEdz)
	
	for j in range(top.ns):
		
		if CorrectionMode == 0:
			derivs.append(speciesq[j]*(efieldz)/jperev )
			
		if CorrectionMode == 1:
			derivs.append(speciesq[j]/jperev*(efieldz + (speciesbeta[j]*clight*hl_pth[j]/2 - speciesq[j]*state_vector[j+top.ns]**2*bfieldz/4/specieslist[j].mass)*dBdz) )
		
		if CorrectionMode == 2:
			derivs.append(speciesq[j]/jperev*(efieldz - state_vector[j+top.ns]**2/4*d2Edz2 + (speciesbeta[j]*clight*hl_pth[j]/2 - speciesq[j]*state_vector[j+top.ns]**2*bfieldz/4/specieslist[j].mass)*dBdz) )

## build second lot in deriv output (i.e. dsigmadz)

	for i in range(top.ns):
		derivs.append(state_vector[i+2*top.ns]) 
	
## build third lot in deriv output (i.e. sigma'' )

	# generate list of neutralization factors	
	#    although the current structure loops over dictionaries, looping over arrays (defined for this sole purpose)
	#    have been tried and showed no appreciable improvement in speed
	
	species_neut_f   = zeros(top.ns)
	
	for ii in sp.keys():
	  js = sp[ii].js
	  species_neut_f[js] = rho_neut_f(rrr,ii)

	for j in range(top.ns):
		
		scterm = 0
		
		for s in range(top.ns):
			QQQ = (speciesq[j]*speciesI[s])/(2*pi*eps0*specieslist[j].mass*speciesbeta[j]**2*speciesbeta[s]*clight**3)
			scterm += QQQ*(1-species_neut_f[s])*state_vector[j+top.ns]/(state_vector[j+top.ns]**2 + state_vector[s+top.ns]**2)
		
		term1 = (speciesq[j]*-efieldz)/(2*state_vector[j]*jperev) * state_vector[j+2*top.ns]
		
		term2 = (speciesq[j]*-dEdz)/(4*state_vector[j]*jperev) * state_vector[j+top.ns]
		
		term3 = ((speciesq[j]*bfieldz)/(2*specieslist[j].mass*speciesbeta[j]*clight))**2*state_vector[j+top.ns]

		emitterm = ((hl_epsrn[j]/speciesbeta[j])**2 + (hl_pthn[j] /speciesbeta[j])**2) / state_vector[j+top.ns]**3
		
		# top.hepsr equals two times the rms-r-emittance
		
		d2sigmadz2 = term1 + term2 - term3 + scterm + emitterm
		
		derivs.append(d2sigmadz2)
		
		termdata.append([j, rrr, term1, term2, -term3, scterm, emitterm, state_vector[j]])
		
	return derivs




# Set up initial states 
#  * Loop over species and fill the elements corresponding to 
#    their respective "js" values
#  * Must fill the list in the correct order to match other input data

#initialstates = [0]*3*top.ns
initialstates = zeros(3*top.ns) 

for ii in sp.keys():
	  s = sp[ii]
	  js = s.js
## Kinetic energy
	  initialstates[js] = ekin[ii]
## Initial rms-radius
	  initialstates[js + top.ns] = sqrt((s.a0/2)**2 + (s.b0/2)**2)
## Initial envelope angle
	  initialstates[js + 2*top.ns] = (s.a0*s.ap0 + s.b0*s.bp0) / initialstates[js+top.ns]


## Numerical Solution of the Env. Model

psoln = odeint (f, initialstates, sss, hmax = stepsize , mxstep= mxstep, atol = atol, rtol = rtol)





### Set up function used to integrate the Env. Model using real-time WARP data

#from scipy import interpolate

##state_vector_2 = [0]*3*top.ns
#state_vector_2 = zeros(3*top.ns)

#deltaz = stepsize/2.

#last_step_number = int((env_ze-env_zs)/ds_diag) + 1 # employed in out-of-range interpolation

#zlist = []
		
#zlist = array([top.hzbeam[kkk] for kkk in range(0, last_step_number + 1)])

#def fwarp(state_vector_2, rrr):
	
	#efieldz = getappliedfields(0, 0, rrr)[2][0]
	#bfieldz = getappliedfields(0, 0, rrr)[5][0]
	
	#dEdz = (getappliedfields(0, 0, rrr + deltaz/2)[2][0] - getappliedfields(0, 0, rrr - deltaz/2)[2][0])/deltaz
	##dVdz = 0
	
	#derivs = [0]*top.ns
	
	#for ii in sp.keys():
	  #s = sp[ii]
	  #js = s.js
	  #derivs[js] = s.charge/jperev*efieldz
	
	#for i in range(top.ns):
		#derivs.append(state_vector_2[i+2*top.ns]) #build second lot in deriv output
	
	#speciesbeta = []
	
	#for i in range(top.ns):
		#speciesbeta.append(sqrt((2*state_vector_2[i]*jperev)/(specieslist[i].mass*clight**2)))

	## generate list of neutralization factors
	
	#species_neut_f   = zeros(top.ns)
	
	#for ii in sp.keys():
	  #js = sp[ii].js
	  #species_neut_f[js] = rho_neut_f(rrr,ii)

	#for j in range(top.ns):
		
		#scterm = 0
		
		#emittancelist = []
		
		##for kkk in range(len(top.hepsr)):
			##emittancelist = emittancelist.append(top.hepsr[0,kkk,j])
		
		#emittancelist = array([hl_epsrn[kkk,j] for kkk in range(0, last_step_number+1)])
		
		#pthetaLIST = []
		
		#pthetaLIST = array([hl_pthn[kkk,j] for kkk in range(0, last_step_number+1)])
		
		#kineticenergylist = []
		
		#kineticenergylist = array([hl_ekin[kkk,j] for kkk in range(0, last_step_number+1)])
		
		#emitinter = interpolate.interp1d(zlist, emittancelist, kind='slinear')
		
		#pthetainter = interpolate.interp1d(zlist, pthetaLIST, kind='slinear')
		
		#keinter = interpolate.interp1d(zlist, kineticenergylist, kind='slinear')	
			
		#if rrr <= env_zs:
			#emittance_j = hl_epsrn[j]
			#ptheta_j = hl_pthn[j]
			#ke_j = hl_ekin[j]
		#elif rrr >= env_ze:
			#emittance_j = hl_epsrn[last_step_number,j]
			#ptheta_j = hl_pthn[last_step_number,j]
			#ke_j = hl_ekin[last_step_number,j]		
		#else:
			#emittance_j = emitinter(rrr)
			#ptheta_j = pthetainter(rrr)
			#ke_j = keinter(rrr)				
	
		#for s in range(top.ns):
			#QQQ = (speciesq[j]*speciesI[s])/(2*pi*eps0*specieslist[j].mass*speciesbeta[j]**2*speciesbeta[s]*clight**3)
			#scterm += QQQ*(1-species_neut_f[s])*state_vector_2[j+top.ns]/(state_vector_2[j+top.ns]**2 + state_vector_2[s+top.ns]**2)
		
		#term1 = (speciesq[j]*-efieldz)/(2*state_vector_2[j]*jperev) * state_vector_2[j+2*top.ns]
		
		#term2 = (speciesq[j]*-dEdz)/(4*state_vector_2[j]*jperev) * state_vector_2[j+top.ns]
		
		#term3 = ((speciesq[j]*bfieldz)/(2*specieslist[j].mass*speciesbeta[j]*clight))**2*state_vector_2[j+top.ns]

		#emitterm = ((emittance_j/speciesbeta[j])**2 + (ptheta_j /speciesbeta[j])**2) / state_vector_2[j+top.ns]**3
		
		#d2sigmadz2 = term1 + term2 - term3 + scterm + emitterm
		
		#derivs.append(d2sigmadz2)
		
	#return derivs
		
		
#if integratewarp == 1:
	#psolnwarp = odeint (fwarp, initialstates, sss, hmax = stepsize, mxstep= mxstep, atol = atol, rtol = rtol)
