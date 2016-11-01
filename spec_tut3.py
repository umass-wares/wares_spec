import time
import struct
from corr import katcp_wrapper
import numpy as np
from adc5g_cal import ADC5g_Calibration_Tools

roach = katcp_wrapper.FpgaClient('172.30.51.97')
roach.wait_connected()
roach.progdev('tut3_adc5g_2016_Oct_13_1431.bof.gz')

adc_cal_tools = ADC5g_Calibration_Tools(roach, program=False)

nbram = 2
numchannels = 2048 
gain = 0xffff
sync_time = 167772160/(roach.est_brd_clk()*1e6)
acc_len= (sync_time)*(200e6)/(1024)

def configure():
        
        print 'Configuring accumulation period to %d...' %acc_len,
        roach.write_int('acc_len', acc_len)
        print 'done'

        #time.sleep(0.5)
        print 'Setting digital gain of all channels to %i...' %gain,
        roach.write_int('gain', gain)
	
	sync_time = 167772160/(roach.est_brd_clk()*1e6)
	
	reset()
        time.sleep(0.05)
        
def adc_cal():

        print '------------------------'
        print 'Loading default OGP/INL corrections to ADCs'
        adc_cal_tools.update_ogp()
        adc_cal_tools.update_inl()
        
def reset():

        print 'Resetting counters...',
        roach.write_int('cnt_rst',1) 
        roach.write_int('cnt_rst',0) 
        print 'done'


def get_acc_n():

        acc_n = roach.read_uint('acc_cnt')
        return acc_n

def get_sync_cnt():

        sync_n = roach.read_uint('sync_cnt')
        return sync_n    

def read_bram(chan, bramNo):
        
	name = None

	if bramNo == 0:
		name = 'even'

	if bramNo == 1:
		name = 'odd'

        bram_array = np.array(struct.unpack('>%dl' %(numchannels/nbram), 
                                     roach.read('%s' %(name),
                                     4.*(numchannels/nbram), 0)),
                                     dtype='float')

        return bram_array

    
def integrate(chan, noisebox_atten, plt=True):

        #reset()
	#time.sleep(integtime)

        t1 = time.time()
        
	acc_n = get_acc_n()
	sync_n = get_sync_cnt()

        interleave_a = np.empty((numchannels,), dtype=float)

        for bramNo in range(nbram):

            interleave_a[bramNo::nbram] = read_bram(chan, bramNo)
	
        read_time = time.time() - t1

	print 'Done with integration'
	print 'acc_n = %i, sync_n = %i' %(acc_n, sync_n)
	
	if plt:
		plot(10.*np.log10(interleave_a[10:]), label='%s' %(noisebox_atten))
		xlabel('FFT Channel')
		ylabel('dB')
		                 
        return acc_n, sync_n, interleave_a, read_time
   
def snap(hist=True):

	raw = adc_cal_tools.get_snap(0, None)
	
	if hist:
		adc_cal_tools.levels_hist(raw)

	return raw

	
# MAIN
configure()
adc_cal()
