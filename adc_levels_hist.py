from adc5g_cal import ADC5g_Calibration_Tools
import numpy as np
from corr import katcp_wrapper
import adc5g

roach = katcp_wrapper.FpgaClient('172.30.51.97')
roach.wait_connected()
roach.progdev('adc5g_tim_aleks_test_2015_Oct_14_1208.bof.gz')

adc_cal = ADC5g_Calibration_Tools(roach)

raw = adc5g.get_snapshot(roach, 'scope_raw_a0_snap')

adc_cal.levels_hist(raw)

