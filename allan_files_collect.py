import hp8780a
import numpy as np
import os
import time
from spec_base import FPGA

Nout = 4000
sampleint = 0.5

syn = hp8780a.HP8780A()
syn.output_off

noise = 29

fpga = FPGA('172.30.51.97')
fpga.configure()

spec_matrix = np.zeros((Nout,2048))
times = np.zeros(Nout)


for i in range(0, Nout):
	times[i] = time.time()
	spec_matrix[i] = fpga.integrate_for_allan(sampleint)

	#time = '#%s' %(times[i])
	#dex = i + 1
	#np.savetxt('allan_%ddbnoise_q_%d.txt' % (noise, dex), spec_matrix[i],
	#	   comments = time)
	

np.savez('allan_%ddbnoise_q' %noise, spec_matrix = spec_matrix, times = times)
 


