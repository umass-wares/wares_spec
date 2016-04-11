import spec_base
import hp8780a
import numpy as np
import os
import time
from spec_base import FPGA

power = 4.0

syn = hp8780a.HP8780A()
freq = 100
freqset = freq*1000000.0
syn.set_freq(freqset)

noise = 29
frequency = (800.0/2048)*np.arange(2048)
fpga = FPGA('172.30.51.97')
fpga.configure()

syn.output_off()
#syn.set_power_level(power)

for i in range(0, 4000):
	#time.sleep(0.5)
	data = fpga.integrate_for_allan(0.5)
	dex = i + 1
	np.savetxt('unquant_rfoff_%d_dbnoise_%d.txt' % (noise, dex), data)
	
