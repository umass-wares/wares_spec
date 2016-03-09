from adc5g_cal import ADC5g_Calibration_Tools
import numpy as np
from corr import katcp_wrapper

roach = katcp_wrapper.FpgaClient('172.30.51.97')
adc_cal = ADC5g_Calibration_Tools(roach, clk=1600)

df = np.load('ogp_noise_default.npz')

zdok0_ogp = df['zdok0_ogp']

zdok = 0
cores = [1,2,3,4]

print
print "Setting ogp for zdok0..."
print zdok0_ogp
print

#Set the ogp
adc_cal.set_ogp(zdok0_ogp, zdok, cores)
