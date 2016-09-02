#!/usr/bin/env python

import corr, logging, struct
import numpy as np
import time
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.animation as animation
import sys
from scheduler import Task, Scheduler
import math
from adc5g_cal import ADC5g_Calibration_Tools

test_bof = 'adc5g_tim_aleks_test_2015_Oct_14_1208.bof.gz'

# 800 MHz
bof_q_800 = 'rspec_1600mhz_r112_asiaa_adc_6_sb1_2016_Jul_14_1242.bof.gz'
bof_uq_800 = 'rspec_1600mhz_r112_asiaa_adc_5_2016_Feb_05_1528.bof.gz'

# 400 MHz
bof_q_400 = 'rspec_800mhz_r112_asiaa_adc_6_2016_May_16_1549.bof.gz'

# 200 MHz 
# bof_q_200

katcp_port = 7147
intnum = 1

class Spec(object):

    def __init__(self, roach_id, katcp_port=7147,
                 bitstream = bof_q_800,
                 cal_bitstream = test_bof,
                 skip = False,
                 quant = True,
                 acc_len = int(.2*(2**28)/2048),
                 gain = 0xffff,
                 shift = 0xffff,
                 numchannels = 2048,
                 bandwidth = 800.0 
                 ):

        self.roach_id = roach_id
        self.fpga = None
        self.katcp_port = katcp_port
        self.bitstream = bitstream
        self.cal_bistream = cal_bitstream
        self.skip = skip

        self.quant = quant
        self.acc_len = acc_len
        self.gain = gain
        self.shift = shift
        self.acc_n = 0

        self.numchannels = numchannels
        self.bandwidth = bandwidth
        self.freq = np.arange(self.numchannels)\
            *float(self.bandwidth)/self.numchannels
        self.set_logging()

        if self.quant:
            self.bram_size = 4 #bytes
            #self.format = '>%dl' %(self.bram_vec_len)
        else:
            self.bram_size = 8 #bytes
            #self.format = '>%dq' %(self.bram_vec_len)

        if self.bandwidth == 200.0:
            self.nbram = 2
        else:
            self.nbram = 4

        self.connect()
        self.program()

        self.adc_cal_tools = ADC5g_Calibration_Tools(self.fpga, program=False)
        self.adc_cal()

        self.configure()

        
    def set_logging(self):

        self.lh = corr.log_handlers.DebugLogHandler()
        self.logger = logging.getLogger(self.roach_id)
        self.logger.addHandler(self.lh)
        self.logger.setLevel(10)

    def connect(self):

        print('Connecting to server %s on port %i... ' % (self.roach_id, 
                                                          self.katcp_port)),

        self.fpga = corr.katcp_wrapper.FpgaClient(self.roach_id, self.katcp_port,
                                                  timeout=10,
                                                  logger=self.logger)
        time.sleep(1)
        if self.fpga.is_connected():
            print 'ok\n'
        else: 
            print 'ERROR connecting to server %s on port %i.\n' % (self.roach_id, self.katcp_port)
            self.exit_fail()

    def program(self):

        print '------------------------'
        print 'Programming FPGA with %s...' % self.bitstream,
        if not self.skip:
            self.fpga.progdev(self.bitstream)
            print 'done'
        else:
            print 'Skipped.'

    def adc_cal(self):

        print '------------------------'
        print 'Loading default OGP/INL corrections to ADCs'
        self.adc_cal_tools.update_ogp()
        self.adc_cal_tools.update_inl()

    def configure(self):
        
        print 'Configuring accumulation period to %d...' % self.acc_len,
        self.fpga.write_int('acc_len', self.acc_len)
        print 'done'

        #time.sleep(0.5)
        #print 'Setting digital gain of all channels to %i...' % self.gain,
        if self.quant and not self.skip:
            self.fpga.write_int('gain', self.gain) 
            print 'done'
        else:   
            print 'Skipped.'

        #time.sleep(0.5)
        print 'Setting Shift schedule for FFT to %i...' % self.shift,
        if not self.skip:
            self.fpga.write_int('fft_shift', self.shift)
            print 'done'
        else:
            print 'Skipped.'
        #time.sleep(0.5)

        self.reset()
        time.sleep(0.05)
        
    def reset(self):

        print 'Resetting counters...',
        self.fpga.write_int('cnt_rst',1) 
        self.fpga.write_int('cnt_rst',0) 
        print 'done'

    def stop(self):

        print 'Stopping accumulation...',
        self.fpga.write_int('cnt_rst',1)

    def go(self):

        print 'Starting accumulation...',
        self.fpga.write_int('cnt_rst',0)
        
    def exit_fail(self):

        print 'FAILURE DETECTED. Log entries:\n', self.lh.printMessages()
        try:
            self.fpga.stop()
        except: pass
        raise
        sys.exit()

    def exit_clean(self):

        try:
            self.fpga.stop()
        except: pass
        sys.exit()

    def get_acc_n(self):

        self.acc_n = self.fpga.read_uint('acc_cnt')
        return self.acc_n

    def get_sync_cnt(self):

        self.acc_n = self.fpga.read_uint('sync_cnt')
        return self.acc_n    


    def read_bram(self, chan, bramNo):
        
        bram_array = np.array(struct.unpack('>%dl' %(self.numchannels/self.nbram), 
                                     self.fpga.read('bram%i%i' %(chan, bramNo),
                                     self.bram_size*(self.numchannels/self.nbram), 0)),
                                     dtype='float')

        return bram_array

    

    def get_data(self, chan, read_acc_n=False):

        """
           For quantized 32-bit vacc bram (unquantized 64-bit vacc needs '>%dq')                                         
           Retrieves data from channel 'chan'                                                                            
           Interleaves bram blocks into full spectrum 
        """

        t1 = time.time()

        if read_acc_n:
            acc_n = self.get_acc_n()
        
        interleave_a = np.empty((self.numchannels,), dtype=float)

        for bramNo in range(self.nbram):

            interleave_a[bramNo::self.nbram] = self.read_bram(chan, bramNo)

        read_time = time.time() - t1
        
        if read_acc_n:
            print "acc_n = %d; Got data in %s secs" % (acc_n, read_time)
            return acc_n, interleave_a, read_time

        else:
            print "read Channel %i in %s secs" %(chan, read_time)
            return interleave_a, read_time


    def get_data_all(self):
        
        """
           Gets data from all 4 channels
        """
        
        t1 = time.time()
        data = {}

        for i in range(4):

            acc_n, dt = self.get_data(i)
            data[i] = dt

        print "Got all data in %s secs" % (time.time() - t1)            
        return data



    def get_data_task(self):

        acc_n, data = self.get_data()
        print acc_n, self.prev_acc_n
        if self.accum_count > 0:
            if acc_n - self.prev_acc_n > 1:
                print "Missed a frame!"
        self.prev_acc_n = acc_n
        self.accum_count += 1
        self.integtime += self.dumptime
        print "Accum count: %s; Integ Time: %s" % (self.accum_count, self.integtime)
        self.data += data
        
        #Save Data        

        if self.integtime >= self.required_integ_time:
            fpga.stop()
            print "Stopping scheduler"
            self.sched.StopAllTasks()
            np.savetxt('data.txt', self.data/float(self.accum_count))
        
    def integrate(self, init_delay, dumptime, integtime):
        """
        Start an integration. Sleep for init_delay, dump every dumptime
        seconds, and integrate for integtime. Uses a scheduler and checks
        for acc_n to increase linearly
        """
        self.dumptime = dumptime
        self.required_integ_time = integtime
        self.integtime = 0.0
       # self.acc_len = int(1.5*dumptime*(2**28)/1024) #make sure that dumptime is longer than read interval
        self.data = np.zeros(self.numchannels, dtype='float')
        self.accum_count = 0
        self.prev_acc_n = 0
        self.sched = Scheduler()
        self.sched.AddTask(self.get_data_task, self.dumptime, init_delay)
        #self.skip = False #ensure gain and acc_len are written
        self.configure()
        self.sched.StartAllTasks()


    def integrate_blocking(self, init_delay, dumptime, integtime):
        """
        Start an integration. Sleep for init_delay, dump every dumptime
        seconds, and integrate for integtime. 
        """
        self.dumptime = dumptime
        self.required_integ_time = integtime
        self.integtime = 0.0
       #self.acc_len = int(1.5*dumptime*(2**28)/1024) #make sure that dumptime is longer than read interval
        self.data = np.zeros(self.numchannels, dtype='float')
        self.accum_count = 0
        self.prev_acc_n = 0
        self.configure()
        time.sleep(init_delay)

        while self.integtime < self.required_integ_time:
            acc_n, data = self.get_data_fast()
            #fp = open('data_%d.txt' % self.accum_count, 'wb')
            #fp.write(data)
            #fp.close()
            print acc_n, self.prev_acc_n
            if self.accum_count > 0:
                if acc_n - self.prev_acc_n > 1:
                    print "Missed a frame!"
            self.prev_acc_n = acc_n
            self.accum_count += 1
            self.integtime += self.dumptime
            print "Accum count: %s; Integ Time: %s" % (self.accum_count, self.integtime)
            self.data += data
            time.sleep(self.dumptime)

        # Done with integ
        fpga.stop()
        print "Stopping"
        np.savetxt('data.txt', self.data/float(self.accum_count))
    
    def integrate_single_accum(self, init_delay, dumptime, integtime): 

        self.dumptime = dumptime
        self.required_integ_time = integtime
        self.integtime = 0.0
        self.data = np.zeros(self.numchannels, dtype='float')
        self.accum_count = 0
        self.prev_acc_n = 0

        self.configure()

        while self.integtime < self.required_integ_time:
            t1 = time.time()

            self.reset()
            time.sleep(self.dumptime)
            acc_n, data = self.get_data()
            #fp = open('data_%d.txt' % self.accum_count, 'wb')
            #fp.write(data)
            #fp.close()
            print acc_n, self.prev_acc_n
            if self.accum_count > 0:
                if acc_n - self.prev_acc_n > 1:
                    print "Missed a frame!"
            self.prev_acc_n = acc_n
            self.accum_count += 1
            self.integtime += self.dumptime
            print "Accum count: %s; Integ Time: %s" % (self.accum_count, self.integtime)
            self.data += data
            
            t2 = time.time()

            print "Loop time: %s \n" % (t2-t1)

        fpga.stop()
        print "Stopping"
        np.savetxt('data.txt', self.data/float(self.accum_count))


    def integrate_for_allan(self, chan, integtime, read_acc_n=False):

        self.reset()
        start_time = time.time()
        time.sleep(integtime)

        if read_acc_n:
            acc_n, data, read_time = self.get_data(chan, read_acc_n=read_acc_n)
            return acc_n, data, start_time, read_time 
        
        else:
            data, read_time = self.get_data(chan, read_acc_n=read_acc_n)
            return data, start_time, read_time 

    def integrate_all(self, integtime):

        t1 = time.time()
        self.reset()
        time.sleep(integtime)
        data = {}
        for i in range(4):
            data[i], read_time = self.get_data(i, read_acc_n=False)
        print "Took %s seconds for all 4 pixels" % (time.time() - t1)
        return data
    
    def map_dump(self, dumptime, maptime, chans=[0, 1, 2, 3]):

        t0 = time.time()
        ndmp = int(maptime/dumptime)
        mapdata = np.zeros((len(chans), self.numchannels, ndmp),
                              dtype='float')
        for i in range(ndmp):
            self.reset()
            time.sleep(dumptime)
            t1 = time.time()
            for chan in chans:
                mapdata[chan, :, i], read_time = self.get_data(chan, read_acc_n=False)
            print "Took %s seconds to read out %d pixels" % (time.time() - t1, len(chans))
        print "Finished Map in %s seconds" % (time.time() - t0)
        return mapdata
        
        
    def do_sfdr_sinad(self, power, tone):
	    """
	    Calculates the sinad and sfdr of the raw data. It first gets the PSD.
	    Then, finds the fundamental peak and the maximum spurious peak. 
	    (Harmonics are also spurs). The DC level is excluded.
	    """

	    freqs = self.freq
	    db = 10*np.log10(power)
	    tot_pwr = 0.0
	    in_peak = False
	    spur_pwr = 0.0
	    peak_freq = 0.0
	    pwr_in_peak = 0.0
	    peak_db = 0.0
	    sig_pwr=0.0
     
	    for i in range(4,len(freqs)):
		    if abs(freqs[i] - tone) < 4:
			    test = -70
		    else:
			    test = -90
		    pwr = 10**(float(db[i])/10.)
		    tot_pwr += pwr

		    if in_peak:
			    if db[i] < test:
				    in_peak = False
				    if abs(peak_freq - tone) < 1:
					    sig_pwr = pwr_in_peak
				    else:
					    if pwr_in_peak > spur_pwr:
						    spur_pwr = pwr_in_peak	
			    else:
				    pwr_in_peak += pwr
				    if db[i] > peak_db:
					    peak_db = db[i]
					    peak_freq = freqs[i]
		    elif db[i] > test:
			    pwr_in_peak = pwr
			    peak_freq = freqs[i]
			    peak_db = db[i]
			    in_peak = True

	    try:
    	        sfdr = 10.0*math.log10(sig_pwr / spur_pwr)
    	        sinad = 10.0*math.log10(sig_pwr/(tot_pwr - sig_pwr))
	    except:
    	        sfdr=0
    	        sinad=0         
    	        print ("Math error at freq %f Mhz." %tone)
    	        print ("sig_pwr=%f, spur_pwr= %f, tot_pwr=%f" %(sig_pwr, spur_pwr, tot_pwr))
	    return sfdr, sinad
    
    
    
class Integration(object):

    def __init__(self, intnum):
        self.intnum = intnum
        self.acc_n = 0
        self.total_acc_n = 0

    def update_line(frames, fpga, integ, line, ax):
        acc_n, data = fpga.get_data()
        integ.intnum += 1
        integ.acc_n = acc_n
        integ.total_acc_n += integ.acc_n
        print "acc_n = %d, intnum=%d, total = %d" % (integ.acc_n, integ.intnum, integ.total_acc_n)
        line.set_ydata(np.array(data))
        ax.set_title("Integration number= %d, acc_n = %d" % (integ.intnum, integ.acc_n))
        plt.draw()
        return line,

#
#START OF MAIN:
#
if __name__ == '__main__':
    from optparse import OptionParser


    p = OptionParser()
    p.set_usage('Usage: %prog <ROACH_HOSTNAME_or_IP> [options]')
    p.set_description(__doc__)
    p.add_option('-q', '--quant', dest='quant', action='store_true',
                 help='Quantized (True) or Unquantized (False) spectrometer bof')

    p.add_option('-g', '--gain', dest='gain', type='int',
                 default=0xffff,
                 help='Set the digital gain (6bit quantisation scalar).Default is 0xffffffff (max), good for wideband noise. Set lower for CW tones.')

    p.add_option('-s', '--skip', dest='skip', action='store_true',
                 help='Skip reprogramming the FPGA and configuring EQ.')

    p.add_option('-b', '--bof', dest='boffile',
                 type='str', default='',
                 help='Specify the bof file to load')

    p.add_option('-n', '--chans', dest='chans',
                 type='int', default='2048',
                 help='Number of frequency channels in bof')

    p.add_option('-l', '--acc_len', dest='acc_len',
                 type='float', default='0.3',
                 help='Set the number of vectors to accumulate between dumps. Format is  t*(2^28)/numchannels, where user input t corresponds to acc_len ~ t seconds')

    p.add_option('-y', '--delay', dest='init_delay',
                 type='float', default='0.2',
                 help='Initial delay after beginning accumulation (sec)')

    p.add_option('-d', '--dumptime', dest='dumptime',
                 type='float', default='0.2',
                 help='Time between bram data dumps (sec)')

    p.add_option('-i', '--integtime', dest='integtime',
                 type='float', default='1.0',
                 help='Total integration period before reset/stop (sec)')
                 
    opts, args = p.parse_args(sys.argv[1:])

    if args==[]:
        print 'Please specify a ROACH board. Run with the -h flag to see all options.\nExiting.'
        exit()
    else:
        roach = args[0]
        
    if opts.boffile != '':
        bitstream = opts.boffile

    #fpga = FPGA(roach, katcp_port, bitstream = bitstream, quant=opts.quant,
    #            skip=opts.skip, acc_len=opts.acc_len, gain=opts.gain)

    chans = opts.chans
    acc_len=int(opts.acc_len*(2.**28)/2048.)
    print 'acc_len = %g' %acc_len

    fpga = FPGA('172.30.51.97', acc_len=acc_len)

    init_delay = opts.init_delay
    dumptime = opts.dumptime
    integtime = opts.integtime

    fpga.integrate_single_accum(init_delay, dumptime, integtime)
 
    #set up the figure with a subplot to be plotted
    #fig = matplotlib.pyplot.figure()
    #ax = fig.add_subplot(1,1,1)

    #start the process
    #time.sleep(2)
    #fig.canvas.manager.window.after(1000, plot_spectrum, fpga)
    
    # scope = SpecScope(ax, fpga)
    
    #ani = animation.FuncAnimation(fig, scope.update, interval=10,
    #                              blit=True)
    

    #integ = Integration(intnum)
    #freq = numpy.arange(chans)*BW/chans        
    #acc_n, data = fpga.get_data()
    #print acc_n
    #line, = ax.plot(freq, numpy.array(data), linestyle='steps-mid')
    #ax.set_title("Integration number= %d, acc_n = %d" % (intnum, acc_n))
    #ine_ani = animation.FuncAnimation(fig, update_line, 200, fargs=(fpga, integ, line, ax),
    #                                  interval=5000, blit=True)

    #plt.show()
    #plot_spectrum(fpga)
    #matplotlib.pyplot.show()
    #print 'Plot started.'
