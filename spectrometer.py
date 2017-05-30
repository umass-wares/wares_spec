import numpy as np
import struct, datetime, time
import adc5g
from pylab import plot,xlabel,ylabel
from adc5g_cal import ADC5g_Calibration_Tools
from corr import katcp_wrapper
from spectrometer_modes import mode_800, mode_400, mode_200

class Spectrometer(object):

    def __init__(self, roach_id='172.30.51.101', katcp_port=7147, mode=800, scale=1024):

        self.roach_id = roach_id
        self.katcp_port = katcp_port 

        roach = katcp_wrapper.FpgaClient(roach_id)
        roach.wait_connected()

        self.roach = roach

        self.sync_scale = scale
        self.sync_period = None
        self.sync_time = None
        self.acc_len = None

        if (mode==800):
            self.mode = mode_800()

        if (mode==400):
            self.mode = mode_400()

        if (mode==200):
            self.mode = mode_200()

        self.adc_cal_tools = ADC5g_Calibration_Tools(self.roach, program=False)

        self.set_sync_period()
        self.set_acc_len()

        self.program_device()
        self.configure()
        self.calADC()

    def calc_sync_period(self, scale):

        period = ((scale*self.mode.sync_LCM*self.mode.pfb_taps*self.mode.FFTsize)/self.mode.FFTinputs)-2
        return period

    def set_sync_period(self):
    
        self.sync_period = self.calc_sync_period(self.sync_scale)
        self.sync_time = self.sync_period/(self.mode.clk/self.mode.ADCstreams)
    
    def set_acc_len(self):

        self.acc_len = self.sync_period/(self.mode.numchannels/self.mode.nbram)


    def program_device(self):

        bitcode = self.mode.bitcode
        self.roach.progdev(bitcode)

        print 'Programming bitcode %s' %(bitcode)

    def configure(self):

        print 'Configuring accumulation period to %d...' %self.acc_len,
        self.roach.write_int('acc_len', self.acc_len)
        print 'done'
                                                                                       
        print 'Setting digital gain of all channels to %i...' %self.mode.gain,
        self.roach.write_int('gain', self.mode.gain)
        print 'done'

        print 'Setting fft shift schedule to %i...' %self.mode.shift,
        self.roach.write_int('fftshift', self.mode.shift)
        print 'done'

        print 'Setting sync period to %i...' %self.sync_period,
        self.roach.write_int('sync_constant', self.sync_period)
        print 'done'

        self.reset()
        time.sleep(0.1)
    
    def calADC(self):

        print '------------------------'
        print 'Loading default OGP/INL corrections to ADCs'
        self.adc_cal_tools.update_ogp()
        self.adc_cal_tools.update_inl()

    def reset(self):

        print 'Resetting counters...',
        self.roach.write_int('cnt_rst',1)
        self.roach.write_int('cnt_rst',0)
        print 'done'

    def get_acc_n(self):

        acc_n = self.roach.read_uint('acc_cnt')
        return acc_n

    def get_sync_cnt(self):

        sync_n = self.roach.read_uint('sync_cnt')
        return sync_n

    def read_bram(self, inp, bramNo):                                                                                                                                                    

        bram_array = np.array(struct.unpack('>%dl' %(self.mode.numchannels/self.mode.nbram),
                                     self.roach.read('bram%i%i' %(inp,bramNo),
                                     4.*(self.mode.numchannels/self.mode.nbram), 0)),
                                     dtype='float')

        return bram_array
    

    def integrate(self, inp, plt=True, write_nc=True):

        t1 = time.time()

        acc_n = self.get_acc_n()
        sync_n = self.get_sync_cnt()

        interleave = np.empty((self.mode.numchannels,), dtype=float)

        for bramNo in range(self.mode.nbram):

            interleave[bramNo::self.mode.nbram] = self.read_bram(inp, bramNo)

        read_time = time.time() - t1

        print 'Done with integration'
        print 'acc_n = %i, sync_n = %i' %(acc_n, sync_n)

        if plt:
                plot(10.*np.log10(interleave[10:]))
                xlabel('FFT Channel')
                ylabel('dB')

        return acc_n, sync_n, interleave, read_time


    def snap(self, inp, hist=True):
                                                                                                                                     
        raw = adc5g.get_snapshot(self.roach, 'snap%i' %(inp))

        if hist:
                self.adc_cal_tools.levels_hist(raw)

        return raw

#    def snap_bad(self, hist=True):

#        raw = adc5g.get_snapshot(self.roach, 'snap')

#        if hist:
#            self.adc_cal_tools.levels_hist(raw)
        
#        return raw
