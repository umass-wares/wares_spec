#!/usr/bin/env python
'''
This script demonstrates programming an FPGA, configuring a wideband
spectrometer and plotting the received data using the Python KATCP
library along with the katcp_wrapper distributed in the corr
package. Designed for use with TUT3 at the 2009 CASPER workshop.

You need to have KATCP and CORR installed. Get them from
http://pypi.python.org/pypi/katcp and
http://casper.berkeley.edu/svn/trunk/projects/packetized_correlator/corr-0.4.0/

Author: Jason Manley, November 2009.
'''

#TODO: add support for ADC histogram plotting.
#TODO: add support for determining ADC input level 

import corr,time,numpy,struct,sys,logging,pylab,matplotlib
import time
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.animation as animation
import sys
from scheduler import Task, Scheduler

bitstream = 'r_spec_1600mhz_r112_asiaa_adc_2015_Nov_24_1220.bof.gz'
katcp_port = 7147

intnum = 1

class FPGA(object):
    def __init__(self, roach_id, katcp_port=7147,
                 bitstream = bitstream,
                 skip = True,
                 acc_len = int(0.1*(2**28)/1024),
                 gain = 0xff,
                 shift = 0xaaa,
                 numchannels = 2048,
                 bandwidth = 800.0
                 ):
        self.roach = roach_id
        self.katcp_port = katcp_port
        self.bitstream = bitstream
        self.skip = skip
        self.acc_len = acc_len
        self.gain = gain
        self.shift = shift
        self.fpga = None
        self.acc_n = 0
        self.numchannels = numchannels
        self.bandwidth = bandwidth
        self.freq = numpy.arange(self.numchannels)*self.bandwidth/self.numchannels
        self.set_logging()
        self.connect()
        self.program()
        self.configure()
        
    def set_logging(self):
        self.lh = corr.log_handlers.DebugLogHandler()
        self.logger = logging.getLogger(self.roach)
        self.logger.addHandler(self.lh)
        self.logger.setLevel(10)

    def connect(self):
        print('Connecting to server %s on port %i... ' % (self.roach, self.katcp_port)),
        self.fpga = corr.katcp_wrapper.FpgaClient(self.roach, self.katcp_port,
                                                  timeout=10,
                                                  logger=self.logger)
        time.sleep(1)
        if self.fpga.is_connected():
            print 'ok\n'
        else: 
            print 'ERROR connecting to server %s on port %i.\n' % (self.roach, self.katcp_port)
            self.exit_fail()

    def program(self):
        print '------------------------'
        print 'Programming FPGA with %s...' % self.bitstream,
        if not self.skip:
            self.fpga.progdev(self.bitstream)
            print 'done'
        else:
            print 'Skipped.'

    def configure(self):
        print 'Configuring accumulation period to %d...' % self.acc_len,
        self.fpga.write_int('acc_len', self.acc_len)
        print 'done'


        #time.sleep(0.5)
        print 'Setting digital gain of all channels to %i...' % self.gain,
        if not self.skip:
            self.fpga.write_int('gain', self.gain) #write the same gain for all inputs, all channels
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
    
    def get_data(self):
        t1 = time.time()
        acc_n = self.get_acc_n()
        
        #
        # Retrieves data from bram0-bram3 registers
        # Interleaves into full spectrum
        #
        a_0 = numpy.array(struct.unpack('>%dl' % (self.numchannels/4.), 
                                        self.fpga.read('bram0',
                                                       self.numchannels, 0)
                                        ),
                          dtype='float')

        a_1 = numpy.array(struct.unpack('>%dl' % (self.numchannels/4.), 
                                        self.fpga.read('bram1', 
                                                       self.numchannels, 0)
                                        ),
                          dtype='float')

        a_2 = numpy.array(struct.unpack('>%dl' % (self.numchannels/4.),
                                        self.fpga.read('bram2',
                                                       self.numchannels, 0)
                                        ),
                          dtype='float')

        a_3 = numpy.array(struct.unpack('>%dl' % (self.numchannels/4.),
                                        self.fpga.read('bram3',
                                                       self.numchannels, 0)
                                        ),
                          dtype='float')
      
        interleave_a = numpy.empty((self.numchannels,), dtype=float)
        
        interleave_a[0::4] = a_0
        interleave_a[1::4] = a_1
        interleave_a[2::4] = a_2
        interleave_a[3::4] = a_3

        print "acc_n = %d; Got data in %s secs" % (acc_n, time.time() - t1)
        #return acc_n, numpy.array(interleave_a, dtype=float)
        return acc_n, interleave_a

    
    def get_data_old(self):
        t1 = time.time()        
        acc_n = self.get_acc_n()
        
        a_0 = numpy.array(struct.unpack('>1024l', self.fpga.read('even',1024*4,0)),
                          dtype=float)
        a_1 = numpy.array(struct.unpack('>1024l', self.fpga.read('odd',1024*4,0)),
                          dtype=float)
        
        interleave_a = numpy.empty((2048,), dtype=float)

        #for i in range(1024):
        #    interleave_a.append(a_0[i])
        #    interleave_a.append(a_1[i])
        interleave_a[0::2] = a_0
        interleave_a[1::2] = a_1

        print "Got data in %s secs" % (time.time() - t1)        
        return acc_n, interleave_a

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
        if self.integtime >= self.required_integ_time:
            print "Stopping scheduler"
            self.sched.StopAllTasks()
        
        
    def integrate(self, init_delay, dumptime, integtime):
        """
        Start an integration. Sleep for init_delay, dump every dumptime
        seconds, and integrate for integtime. Uses a scheduler and checks
        for acc_n to increase linearly
        """
        self.dumptime = dumptime
        self.required_integ_time = integtime
        self.integtime = 0.0
        self.acc_len = int(1.5*dumptime*(2**28)/2048) #make sure that dumptime is longer than read interval
        self.data = numpy.zeros(2048, dtype='float')
        self.accum_count = 0
        self.prev_acc_n = 0
        self.sched = Scheduler()
        self.sched.AddTask(self.get_data_task, self.dumptime, init_delay)
        self.skip = False #ensure gain and acc_len are written
        self.configure()
        self.sched.StartAllTasks()
        
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
    line.set_ydata(numpy.array(data))
    ax.set_title("Integration number= %d, acc_n = %d" % (integ.intnum, integ.acc_n))
    plt.draw()
    return line,

#START OF MAIN:

if __name__ == '__main__':
    from optparse import OptionParser


    p = OptionParser()
    p.set_usage('Usage: %prog <ROACH_HOSTNAME_or_IP> [options]')
    p.set_description(__doc__)
    p.add_option('-l', '--acc_len', dest='acc_len',
                 type='int', default=2*(2**28)/2048,
                 help='Set the number of vectors to accumulate between dumps. default is 2*(2^28)/2048, or just under 2 seconds.')
    p.add_option('-g', '--gain', dest='gain', type='int',
                 default=0xffffffff,
                 help='Set the digital gain (6bit quantisation scalar). Default is 0xffffffff (max), good for wideband noise. Set lower for CW tones.')
    p.add_option('-s', '--skip', dest='skip', action='store_true',
                 help='Skip reprogramming the FPGA and configuring EQ.')
    p.add_option('-b', '--bof', dest='boffile',
                 type='str', default='',
                 help='Specify the bof file to load')
    opts, args = p.parse_args(sys.argv[1:])

    if args==[]:
        print 'Please specify a ROACH board. Run with the -h flag to see all options.\nExiting.'
        exit()
    else:
        roach = args[0]
        
    if opts.boffile != '':
        bitstream = opts.boffile

    fpga = FPGA(roach, katcp_port, bitstream=bitstream,
                skip=opts.skip, acc_len=opts.acc_len,
                gain=opts.gain)

    #set up the figure with a subplot to be plotted
    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(1,1,1)

    # start the process
    time.sleep(2)
    #fig.canvas.manager.window.after(1000, plot_spectrum, fpga)
    
    #scope = SpecScope(ax, fpga)
    
    #ani = animation.FuncAnimation(fig, scope.update, interval=10,
    #                              blit=True)
    

    integ = Integration(intnum)
    freq = numpy.arange(2048)*400./2048        
    acc_n, data = fpga.get_data()
    print acc_n
    line, = ax.plot(freq, numpy.array(data), linestyle='steps-mid')
    ax.set_title("Integration number= %d, acc_n = %d" % (intnum, acc_n))
    ine_ani = animation.FuncAnimation(fig, update_line, 200, fargs=(fpga, integ, line, ax),
                                      interval=5000, blit=True)

    plt.show()
    #plot_spectrum(fpga)
    #matplotlib.pyplot.show()
    print 'Plot started.'

