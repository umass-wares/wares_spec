from adc5g_cal import ADC5g_Calibration_Tools
import numpy as np
from corr import katcp_wrapper

roach = katcp_wrapper.FpgaClient('172.30.51.97')
adc_cal = ADC5g_Calibration_Tools(roach, clk=1600)

#df = np.load('ogp_default.npz')
df = np.load('ogp_update_zdok0.npz')

zdok0_ogp = df['zdok0_ogp']

zdok = 0
cores = [1,2,3,4]

print
print "Setting ogp for zdok0..."
print zdok0_ogp
print

adc_cal.set_ogp(zdok, cores, zdok0_ogp)
