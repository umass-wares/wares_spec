import hp8780a
import numpy as np
import os
import time
from spec_base import FPGA

Nout = 4000
sampleint = 1.0

syn = hp8780a.HP8780A()
syn.output_off

#noise = 20

fpga = FPGA('172.30.51.97', 
	    bitstream='rspec_1600mhz_r112_asiaa_adc_6_sb1_2016_Jul_14_1242.bof.gz',
	    cal=False, gain=0xffff, numchannels=2048, bandwidth=800)

spec_matrix = np.zeros((Nout,2048))
start_t = np.zeros(Nout)
read_t = np.zeros(Nout)

for i in range(0, Nout):

	spec_matrix[i], start_t[i], read_t[i] =\
	         fpga.integrate_for_allan(0, sampleint)

	print 'Done with integ # %i' %(i+1)

np.savez('allan_800mhz_nocal_really', spec_matrix = spec_matrix, start_t = start_t, read_t = read_t)
 


