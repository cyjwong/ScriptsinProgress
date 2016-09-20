import numpy as np

noSClist = np.loadtxt('NoSC.out')
withSClist = np.loadtxt('withSC.out')

len1 = len(noSClist)
len2 = len(withSClist)

min_diff = 10.
min_pos_1 = 0
min_pos_2 = 0

for ii in range(len1):
	for jj in range(len2):
		penalty = abs(noSClist[ii,2] - withSClist[jj,2]) + 5*abs(noSClist[ii,3] - withSClist[jj,3])
		
		if penalty < min_diff:
			min_diff = penalty
			min_pos_1 = ii
			min_pos_2 = jj

print noSClist[min_pos_1, :]
print withSClist[min_pos_2, :]
print min_diff
		
