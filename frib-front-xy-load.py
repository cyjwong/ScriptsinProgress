# Load adjustment to account for beam canonical angular momentum.  
# Used in front end simulation script frib-front-xy.py 


# Beam loading for the 0th birth mode

#if birth_mode == 0, nothing needs to be done
		

# Beam loading for the 1st birth mode

if birth_mode == 1:

	bz0_launch = getappliedfields(x=0.,y=0.,z=z_launch)[5]      # B_z on-axis at simulation launch location 
	
	sp_krot_launch = {}
	sp_krot_v      = {} 
	for ii in sp.keys():
	  s = sp[ii]
	  # --- rigidity 
	  gamma = 1./sqrt(1.-(s.vbeam/clight)**2)
	  brho  = gamma*s.mass*s.vbeam/s.charge
	  # --- rms calculation
	  rms_launch = sqrt(average( (s.xp)**2 + (s.yp)**2 ))
	  # --- rot wavenumbers at launch and in vacuum v
		  
	  # bz0_birth*rms_birth**2 of the species under question:
	  b_r2_j = ptheta[ii] / s.charge*2*s.mass*clight
	  
	  krot_launch = (b_r2_j/rms_launch**2 - bz0_launch)/(2.*brho)
	  krot_v      = b_r2_j/rms_launch**2/(2.*brho)
	  # 
	  sp_krot_launch.update({ii:krot_launch})
	  sp_krot_v.update({ii:krot_v}) 
	  #
	  s.uxp -= krot_launch*s.yp*s.uzp
	  s.uyp += krot_launch*s.xp*s.uzp
		 

# Beam loading for the 2nd birth mode

if birth_mode == 2:

	bz0_launch = getappliedfields(x=0.,y=0.,z=z_launch)[5]      # B_z on-axis at simulation launch location 
	
	sp_krot_launch = {}
	sp_krot_v      = {} 
	for ii in sp.keys():
	  s = sp[ii]
	  # --- rigidity 
	  gamma = 1./sqrt(1.-(s.vbeam/clight)**2)
	  brho  = gamma*s.mass*s.vbeam/s.charge
	  # --- rms calculation
	  rms_launch = sqrt(average( (s.xp)**2 + (s.yp)**2 ))
	  # --- rot wavenumbers at launch and in vacuum v
	  
	  krot_launch = (bz0_birth[ii]*rms_birth[ii]**2/rms_launch**2 - bz0_launch)/(2.*brho)
	  krot_v      = bz0_birth[ii]/(2.*brho)
	  # 
	  sp_krot_launch.update({ii:krot_launch})
	  sp_krot_v.update({ii:krot_v}) 
	  #
	  s.uxp -= krot_launch*s.yp*s.uzp
	  s.uyp += krot_launch*s.xp*s.uzp


# Function that locates the two peaks in the ECR B-field and extracts information in-between
# input: starting position of the ECR B-field data
# return: [average field, distance bewteen two peaks, z_peak1, z_peak2, z_trough]

def field_two_peak(ecr_zmin):
	
	step_size = 0.001
	
	cumulativefield = 0.
	stepcounter = 0
	peak1 = 0.
	peak2 = 0.
	trough = 0.
	zzz = ecr_zmin

	region = 0
	
	while True:
		
		aaa = getappliedfields(0, 0, zzz)[5][0]
		bbb = getappliedfields(0, 0, zzz + step_size)[5][0]
		
		if region == 0:
			if aaa > bbb:
				peak1 = zzz
				region = 1
		
		if region == 1:
			if aaa < bbb:
				trough = zzz
				region = 2		
		
		if region == 1:
			cumulativefield += getappliedfields(0, 0, zzz)[5][0]
			stepcounter += 1
			
		if region == 2:
			if aaa > bbb:
				peak2 = zzz
				break
		
		if region == 2:
			cumulativefield += getappliedfields(0, 0, zzz)[5][0]
			stepcounter += 1	
						
		zzz += step_size
		
	return [cumulativefield / stepcounter, stepcounter*step_size , peak1, peak2, trough]

bfieldinfo = field_two_peak(ecr_zmmin)


# Beam loading for the 3rd birth mode
# Define a number of slots corresponding to even spaced positions between and including the two peaks
# For each species, divide number of particles by slot number and place the quotient into each slot
# Distribute the remainder evenly between the two peaks

if birth_mode == 3:

	peak1 = bfieldinfo[2]
	peak2 = bfieldinfo[3]
	
	slotnum = 100
	
	slot_alpha = []
	slot_B = []
	
	bz0_launch = getappliedfields(x=0.,y=0.,z=z_launch)[5]      # B_z on-axis at simulation launch location 
	
	inj_ang_mom = true
	
	sp_krot_launch = {}
	sp_krot_v      = {} 
	for ii in sp.keys():
	  s = sp[ii]
	  # --- rigidity 
	  gamma = 1./sqrt(1.-(s.vbeam/clight)**2)
	  brho  = gamma*s.mass*s.vbeam/s.charge
	  
	  #Main Allocation (place bulk of particles into slots defined bewteen two peaks)
	  
	  numberinslot = len(s.uxp)/slotnum
	  slotz = linspace(peak1, peak2, slotnum)
	  
	  for jjj in range(0,slotnum):
		  zzz = slotz[jjj]
		  bz0_birth = getappliedfields(x=0.,y=0.,z=zzz)[5]
		  rms_birth = sqrt(average( (s.xp)**2 + (s.yp)**2 ))
		  rms_launch = sqrt(average( (s.xp[jjj*slotnum:(jjj+1)*slotnum])**2 + (s.yp[jjj*slotnum:(jjj+1)*slotnum])**2 ))
		  
		  if ii == 'U33':
			  slot_alpha.append((rms_launch/rms_birth)**2)
			  slot_B.append(bz0_birth)
		  
		  krot_launch = (bz0_birth*rms_birth**2/rms_launch**2 - bz0_launch)/(2.*brho)
		  krot_v      = bz0_birth/(2.*brho)
		  
		  for kkk in range(0, numberinslot):
			  pname = kkk + jjj*numberinslot
			  s.uxp[pname] -= krot_launch*s.yp[pname]*s.uzp[pname]
			  s.uyp[pname] += krot_launch*s.xp[pname]*s.uzp[pname]
	  
	  #Remainder Allocation (Evenly distribute the remainder between two peaks, one in each new slot)
	
	  slotnum2 = len(s.uxp)%slotnum
	  numberinmain = len(s.uxp) - slotnum2
	  
	  if slotnum2 != 0:
	  
		  slotz2 = linspace(peak1, peak2, slotnum2)
		  
		  for jjj in range(0,slotnum2):
			  zzz = slotz2[jjj]
			  bz0_birth = getappliedfields(x=0.,y=0.,z=zzz)[5]
			  rms_birth = sqrt(average( (s.xp)**2 + (s.yp)**2 ))
			  
			  pname = numberinmain + jjj
			  
			  rms_launch = sqrt(average( (s.xp[pname])**2 + (s.yp[pname])**2 ))
			  
			  krot_launch = (bz0_birth*rms_birth**2/rms_launch**2 - bz0_launch)/(2.*brho)
			  krot_v      = bz0_birth/(2.*brho)
			  
			  s.uxp[pname] -= krot_launch*s.yp[pname]*s.uzp[pname]
			  s.uyp[pname] += krot_launch*s.xp[pname]*s.uzp[pname]


# Plots of initial rotation by species at launch point

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
