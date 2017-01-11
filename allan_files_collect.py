#import hp8780a
import numpy as np
#import os
import time
from spec_base import Spec, bof_q_800
from corr import katcp_wrapper

Nout = 50
sampleint = 2.0

#roach = katcp_wrapper.FpgaClient('172.30.51.97')
#roach.progdev(bof_q_800)

spec = Spec('172.30.51.97', bitstream=bof_q_800, numchannels=2048, bandwidth=800.0)

spec.fpga.wait_connected()

spec_matrix = np.zeros((Nout,2048))
start_t = np.zeros(Nout)
read_t = np.zeros(Nout)


for i in range(0, Nout):

	spec_matrix[i], start_t[i], read_t[i]= spec.integrate_for_allan(chan=0, integtime=sampleint, read_acc_n=False)

        #np.savetxt('allan_spectra_800mhz_4000samps_1int_aug24_1_samp%i.txt' %(i+1), spec_matrix[i])  

	print 'Done with integ # %i' %(i+1)


np.savez('allan_spectra_800mhz_4000samps_1int_aug24_1', 
         spec_matrix = spec_matrix, start_t = start_t, read_t = read_t)
 


