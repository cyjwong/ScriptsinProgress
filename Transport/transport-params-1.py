# Solenoid Field Strength (in Telsa)

SolenoidStrength = 0.4

# Lattice Length

PeriodLength = 1.2 # (times 5)

# Number of Periods

NumberofPeriods = 5

# Initial Radius

InitialRadius = 0.008*sqrt(2.)

# set linear or nonlinear fields

fieldtype = "lin"

# set neutralization factor

NeutralizationFactor = 0.75

#set velocity correction method: 0 - no correction, 1 - dBdz only, 2 - dBdz + d2Edz2

CorrectionMode = 1

# Number of particles per grid (default = 400, adjust for quick checks)

GridNumber = 400

# Envelope Plot Limits:
plotxmin = "e"
plotxmax = "e"
plotymin = 0.
plotymax = "e"

# Center of 1st Solenoid

s4p2_zc  = PeriodLength / 2.

# Number of steps per period

nstep_period = 200
