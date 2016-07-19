# Solenoid Field Strength (in Telsa)

SolenoidStrength = 0.4

# Lattice Length

LatticeLength = 1.2 # (times 5)

# Number of Periods

NumberofPeriods = 1

# Initial Radius

InitialRadius = 0.008

# set linear or nonlinear fields

abc = "lin"

# set neutralization factor

NeutralizationFactor = 0.75

#set velocity correction method: 0 - no correction, 1 - dBdz only, 2 - dBdz + d2Edz2

CorrectionMode = 1

# Number of particles per grid (default = 400, adjust for quick checks)

ParticleNumber = 400  #(i.e. n_grid)

# Envelope Plot Limits:
plotxmin = "e"
plotxmax = "e"
plotymin = 0.
plotymax = "e"

# Center position of 1st solenoid

s4p2_zc  = 0.3

#Number of steps per period

nstep = 500
