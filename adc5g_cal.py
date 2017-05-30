
from corr import katcp_wrapper
import adc5g
from hp8780a import HP8780A
import time
import datetime
import numpy as np
from numpy import array, mean
from scipy.optimize import leastsq, curve_fit, fmin
import math
from matplotlib import pyplot as plt
from matplotlib import mlab
from matplotlib.mlab import psd, detrend_mean

#
# Functions used by the fitting routines
#

def syn_func( x, off, a1, a2):
        return off +  a1*np.sin(x) + a2*np.cos(x)

def fitsin(p, sin, cos):
	return p[0] + p[1]*sin + p[2]*cos

def sin_residuals(p, sin, cos, raw):
        res = raw - fitsin(p, sin, cos)
        for i in range(raw.size):
            if raw[i] == -128 or raw[i] == 127:
                res[i] = 0
        return res

def gaussian(x,a,mu,sig): return a*np.exp(-(x-mu)**2 / (2. * sig**2))

def chisq(par, x, y, yerr):
    (a, mu, sig) = par
    return np.sum((gaussian(x,a,mu,sig)-y)**2/yerr**2)



#
# ADC5g Calibration Tools
#
# This class contains functions to calibrate two ASIAA ADC5G samplers, 
# each configured in 2-channel mode, for a total of 4 signal channels
# 
# In this mode, each channel is sampled by two independent adc 'cores'
# which are interleaved to achieve the desired sample rate 
#
# ADCs are denoted zdok0, zdok1 
# Channels are denoted 0,1 (zdok0) 2,3 (zdok1)
# Each channel has associated OGPs and INLs for each core
# Cores are denoted A, B, C, D 
# (in 2-channel mode, A and B sample together, C and D sample together)
#

class ADC5g_Calibration_Tools ():
    
    def __init__(self, roach, program = False, roach_id = '172.30.51.97', 
		 bitstream = 'adc5g_tim_aleks_test_2015_Oct_14_1208.bof.gz',
		 clk=1600):

	    self.roach_id = roach_id
	    self.bitstream = bitstream
	    self.clk = clk
	    self.freqs_cores = np.linspace(100,clk/4,6,endpoint=False)
	    self.freqs_bw = np.arange(50,clk/2,25)
	    self.syn = None
	    self.roach = roach

	    if program:
                    self.program()
                   

    def program(self):
            
         print 'Programming ROACH with calibration bitstream:'
         print '%s' %(self.bitstream)
         self.roach.progdev(self.bitstream)
   
    def get_channel_snap_reg(self, chan):

         """
	    Gets the snap register string based on channel No. (0,1,2,3)
	 """
	    
	 reg = None

	 if chan == 0:
		 reg = 'snap0'

	 if chan == 1:
		 reg = 'snap1'

	 if chan == 2:
		 reg = 'snap2'
	    
	 if chan == 3:
		 reg = 'snap3'

	 return reg
     
    
    def get_channel_core_spi(self, chan):

        """
	   Gets the spi parameters based on channel No. (0,1,2,3)
	"""

        if chan == 0:
                zdok = 0
                cores = range(1,3)

        if chan == 1:
                zdok = 0
                cores = range(3,5)

        if chan == 2:
                zdok = 1
                cores =range(1,3)

        if chan == 3:
                zdok = 1
                cores =range(3,5)

        return zdok, cores

		    
    def get_snap(self, chan, freq): #MHz

        """
	   Takes a snap shot on channel 'chan' of a cw tone
	   of frequency 'freq' (MHz)
	"""

#	if freq != None:

#		self.syn.output_on()
#		self.syn.set_freq(freq*1e6)
#		time.sleep(1)

	reg = self.get_channel_snap_reg(chan)
	raw = array(adc5g.get_snapshot(self.roach, reg))

	return raw, freq


    def fit_snap(self, raw, freq): 

        """
	   Takes a snap 'raw' of cw frequency 'freq', minimizes least squares
	   to find best fit parameters for linear sine model
	"""

	del_phi = 2*math.pi*freq/self.clk

	x  = [del_phi*i for i in range(len(raw))]

	core1 = raw[0::2]
	core2 = raw[1::2]

	x1 = x[0::2]
	x2 = x[1::2]
        
	p0 = [128.,90.,90.]
                
	params1 = curve_fit(syn_func, x1, core1, p0)

	params2 = curve_fit(syn_func, x2, core2, p0)
    
	return params1, params2
    

    def calc_ogp(self, params1, params2, freq):

        """
	   Calculates OGP values for both cores given best fit parameters
	   from a snap fit
	"""

	z_fact = -500.0/256.0
	true_zero = 0.0 * z_fact
	d_fact = 1e12/(2*np.pi*freq*1e6)

	z1 = z_fact * params1[0][0]
	sin1a = params1[0][1]
	cos1a = params1[0][2]
	amp1 = math.sqrt(sin1a**2 + cos1a**2)
	dly1 = d_fact*math.atan2(sin1a, cos1a)

	z2 = z_fact * params2[0][0]
	sin2a = params2[0][1]
	cos2a = params2[0][2]
	amp2 = math.sqrt(sin2a**2 + cos2a**2)
	dly2 = d_fact*math.atan2(sin2a, cos2a)

	avz = (z1+z2)/2.0
	avamp = (amp1+amp2)/2.0
	a1p = 100*(avamp-amp1)/avamp 
	a2p = 100*(avamp-amp2)/avamp 
	avdly = (dly1+dly2)/2.0

	ogp1 = (z1-true_zero, a1p, dly1-avdly)
	ogp2 = (z2-true_zero, a2p, dly2-avdly)

	return ogp1, ogp2
        
    def calc_sinad(self, raw, freq):

        """
	   Calculates the SINAD from the residuals of a fit. 
	"""

	del_phi = 2*math.pi*freq/self.clk

	x  = [del_phi*i for i in range(len(raw))]
	p0 = [128.,90.,90.]
        
	params0 = curve_fit(syn_func, x , raw  , p0)

	off = params0[0][0]
	sin0a = params0[0][1]
	cos0a = params0[0][2]
	amp0 = math.sqrt(sin0a**2 + cos0a**2)
       
	fit0 = syn_func(x, off, sin0a, cos0a)
	ssq0 = 0.0
	data_cnt = len(raw)

	for i in range (data_cnt):
		ssq0 += (raw[i]-fit0[i])**2

	pwr_sinad = (amp0**2) / (2*ssq0/data_cnt)
	sinad = 10.0 * math.log10(pwr_sinad)
        
	return sinad
        
    def do_sfdr(self, raw, freq):

	    """
	    Calculates the sinad and sfdr of the raw data. It first gets the PSD.
	    Then, finds the fundamental peak and the maximum spurious peak. 
	    (Harmonics are also spurs). The DC level is excluded.
	    """

	    samp_freq=self.clk
	    nfft=len(raw)
	    power, freqs = psd(raw, nfft, Fs=samp_freq*1e6, 
			       detrend=detrend_mean, scale_by_freq=True)
	    freqs = freqs/1e6
	    db = 10*np.log10(power)
	    tot_pwr = 0.0
	    in_peak = False
	    spur_pwr = 0.0
	    peak_freq = 0.0
	    pwr_in_peak = 0.0
	    peak_db = 0.0
	    sig_pwr=0.0
     
        
	    for i in range(4,len(freqs)):

		    if abs(freqs[i] - freq) < 4:
			    test = -70
		    else:
			    test = -90

		    pwr = 10**(float(db[i])/10.)
		    tot_pwr += pwr

		    if in_peak:

			    if db[i] < test:
				    in_peak = False

				    if abs(peak_freq - freq) < 1:
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
    	        print ("Math error at freq %f Mhz." %freq)
    	        print ("sig_pwr=%f, spur_pwr= %f, tot_pwr=%f" %(sig_pwr, spur_pwr, tot_pwr))
	    return sfdr, sinad
        
    def do_sfdr_sinad_cw_sweep( self, chans=[0,1], save=False, 
                               fname='sfdr_sinad.npz',
                               freqarray=[]):

        """
        Calculates the SFDR and SINAD from a sweep in frequency from 50 Mhz to
        final_freq Mhz.
        """

	#final_freq = self.clk /2
	if not freqarray:
        	freqs = self.freqs_bw
	else: 
        	freqs=freqarray
      
	multi_sfdr = []
	multi_sinad_psd=[]

	for chan in chans:
        	sfdr_chan = {}
        	sinad_psd_chan= {}
        
        	for i in range(len(freqs)): # MHz
        		freq = freqs[i]
        		raw, f = self.get_snap(chan, freq)
        		sfdr_chan[freq], sinad_psd_chan[freq] = self.do_sfdr(raw,freq)
          
        	multi_sfdr.append(sfdr_chan)
        	multi_sinad_psd.append(sinad_psd_chan)
        
	if save:		
		np.savez(fname, sfdr=multi_sfdr, sinad_psd=multi_sinad_psd)

	return multi_sfdr, multi_sinad_psd


    def do_sinad_cw_sweep( self, chans=[0,1], save=False, fname='sinad.npz',
                          freqarray=[]):

        """
        Calculates the SINAD from a sweep in frequency from 50 Mhz to
        final_freq Mhz. 
        freqarray array of freqs
        """
 	if not freqarray: #if freqarray is empty then...
        	freqs = self.freqs_bw
	else: 
        	freqs = freqarray
 
	multi_sinad=[]
 
	for chan in chans:
        	sinad_chan= {}
        
        	for i in range(len(freqs)): # MHz
        		freq = freqs[i]
        		raw, f = self.get_snap(chan, freq)      
        		sinad_chan[freq] = self.calc_sinad(raw,freq)
      
    		multi_sinad.append(sinad_chan)
            
	if save:		
        	np.savez(fname, sinad=multi_sinad)
    
	return multi_sinad     


    def get_ogp_chan(self, chan):

        """
	   Returns ogp values for a channel 'chan'
	"""

        ogp_chan = []
    
        zdok, cores = self.get_channel_core_spi(chan)

        for core in cores:

            off = adc5g.get_spi_offset(self.roach, zdok, core) 
            gain = adc5g.get_spi_gain(self.roach, zdok, core)
            phase = adc5g.get_spi_phase(self.roach, zdok, core)
            ogp_core = [off, gain, phase]

            ogp_chan.append(ogp_core)

        return ogp_chan

    def get_ogp(self, chans=[0,1,2,3]):

        """
	   Get OGP for all channels in 'chans'
	"""

	for chan in chans:

		ogp_chan = self.get_ogp_chan(chan)

		print
		print 'OGP set on channel %i' %chan
		print
		
		print 'Core 1'
		print ogp_chan[0]
		print

		print 'Core 2'
		print ogp_chan[1]
		print
	

    def set_ogp(self, ogp_chan, chan):

        """
	   Sets ogp for two cores of channel 'chan'
	   multi_ogp is format (ogp1, ogp2)
        """

        zdok, cores = self.get_channel_core_spi(chan)

	i = 0
        for core in cores:
            
            off, gain, phase = ogp_chan[i]

            off_spi = math.floor(.5 + off*255/100.) + 0x80
            adc5g.set_spi_offset(self.roach, zdok, core, float(off))

            gain_spi = math.floor(.5 + gain*255/36.) + 0x80
            adc5g.set_spi_gain(self.roach, zdok, core, float(gain))

            phase_spi = math.floor(.5 + phase*255/28.) + 0x80
            adc5g.set_spi_phase(self.roach, zdok, core, float(phase)*0.65)

	    i += 1

        #self.roach.progdev(self.bitstream)


    def do_ogp_cw_sweep(self, chans = [0,1,2,3], set_ogp=True, 
			save=False, fname='ogp_default.npz'):

        """
	   Generates frequency averaged set of ogp for channels in 'chans'
	   Can save to .npz file and/or set generated OGP values to ADC
	   
	"""

	freqs = self.freqs_cores

	multi_ogp = []

	for chan in chans:

		ogp1 = {}
		ogp2 = {}

		for i in range(len(freqs)): # MHz

			freq = freqs[i]
 
			raw, f = self.get_snap(chan, freq) 
		
			params1, params2 = self.fit_snap(raw, freq)
  
			ogp1[freq], ogp2[freq] = self.calc_ogp(params1, params2, 
							       freq)
		
			ogp1_m = tuple(map(lambda y: sum(y)/float(len(y)), 
					   zip(*ogp1.values())))

			ogp2_m = tuple(map(lambda y: sum(y)/float(len(y)), 
					   zip(*ogp2.values())))
        
			ogp_chan = (ogp1_m, ogp2_m)

			if set_ogp:
				
				self.set_ogp(ogp_chan, chan)
				
			

		multi_ogp.append(ogp_chan)
			
	if save:		
		np.savez(fname, multi_ogp = multi_ogp, chans = chans)       


	return multi_ogp, chans



    def do_ogp_noise_source(self, rpt, chans = [0,1], set_ogp = False, 
			    save=False, fname='ogp_noise_chans01.npz'):

        """
	   Calculates OGP from noise source for channels in 'chans'
	   averaged over 'rpts', can set OGP and/or save to file
	"""

	multi_ogp = []

	for chan in chans:

		ogp1 = {}
		ogp2 = {}

		for i in range(rpt):

			raw, f = self.get_snap(chan, None)
			ogp1[i], ogp2[i] = self.calc_ogp_noise(raw)

		ogp1_m = tuple(map(lambda y: sum(y)/float(len(y)), 
				   zip(*ogp1.values())))

		ogp2_m = tuple(map(lambda y: sum(y)/float(len(y)), 
				   zip(*ogp2.values())))

		ogp_chan = (ogp1_m, ogp2_m) 

		if set_ogp:

			self.set_ogp(ogp_chan, chan)

		multi_ogp.append(ogp_chan)

        if save:
                np.savez(fname, multi_ogp = multi_ogp, chans = chans)

        return multi_ogp, chans



    def repeat_snaps(self, rpt=10):

	    counts = np.zeros((2,4,256))

	    for r in range(rpt):

		    for chan in [0,1,2,3]:

			    zdok,cores = self.get_channel_core_spi(chan)
			    raw,f = self.get_snap(chan, freq=None)
 
			    bc0 = np.bincount((raw[0::2]+128))
			    bc1 = np.bincount((raw[1::2]+128))
			    
			    counts[zdok, cores[0]-1, :len(bc0)] += bc0
			    counts[zdok, cores[1]-1, :len(bc1)] += bc1

	    return counts

    def og(self):
	    
	    counts = self.repeat_snaps(100)

	    #corder = [1, 3, 2, 4]
	    x = np.arange(-128, 128, 1) # integer values assigned to bins                                                            
	    xx = np.linspace(-128, 128, 200)

	    gr = 75
	    doplot = True
	    setog = True

	    for i in [0,1]:
		    means = np.zeros(4)
		    stds = np.zeros(4)
		    for j in [0,1,2,3]:
			    y = counts[i,j]
			    yerr = np.sqrt(1+y+.10*y**2) # 10% systematic error                                                              
			    p0=(np.max(y), 0., 30.)
			    # do fit and ignore first and last bins (saturation)                                                             
			    ret = fmin(chisq, (np.max(y), 0, 40), args=(x[gr:-gr], y[gr:-gr], yerr[gr:-gr]), disp=False)
        
			    if doplot:
				    plt.subplot(4,2,1+4*i+j)
				    iflabel = 'IF%d core %d' % (i,j)
				    statslabel = r'$\mu$:%.1f, $\sigma$:%.1f' % (ret[1], ret[2])
				    # h0 = plt.errorbar(x, y, yerr, fmt='.', label='IF%d core %d' % (i,j))                                       
				    h0 = plt.plot(x, y, '.', label='zdok%d core %d' % (i,j))
				    h1 = plt.plot(xx, gaussian(xx, *ret), label=r'$\mu$:%.1f, $\sigma$:%.1f' % (ret[1], ret[2]))
				    plt.text(0.05, 0.95, iflabel, ha='left', va='top', transform=plt.gca().transAxes)
				    plt.text(0.95, 0.95, statslabel, ha='right', va='top', transform=plt.gca().transAxes)
				    plt.xlim(-128, 128)
				    plt.ylim(0, 1.05 * np.max(counts))
				    plt.yticks([])
				    plt.xticks([])

			    means[j] = ret[1]
			    stds[j] = ret[2]
			    print "IF%d Core %d: mean %5.2f std %5.2f" % (i, j, ret[1], ret[2])
		    avg_std = np.mean(stds) # target std
                                                                                 
		    for j in [0,1,2,3]:
			    #orig_off = adc5g.get_spi_offset(self.roach, i, j+1)
			    #orig_gain = adc5g.get_spi_gain(self.roach, i, j+1)
			    #new_off = orig_off - means[j] * 500./256.
			    new_off = - means[j]*500./256.
			    #new_gain = (100. + orig_gain) * (avg_std / stds[j]) - 100.
			    new_gain = 100.*(avg_std/stds[j])-100.
			    print
			    print new_off
			    print new_gain
			    print

			    if setog:
				    adc5g.set_spi_offset(self.roach, i, j+1, new_off)
				    adc5g.set_spi_gain(self.roach, i, j+1, new_gain)
			    #sol[i,j,0] = new_off
			    #sol[i,j,1] = new_gain

	    #if doplot:
	    #	    plt.suptitle('%s ADC 8bit population\n%s' % (open('/etc/hostname').read().strip(), tag))
	    #	    plt.subplots_adjust(hspace=0, wspace=0)
	    #	    plt.setp(plt.gcf(), figwidth=8, figheight=12)
		    #figfile = 'ogplot-%s.png' % tag
		    #print "saving figure to: %s" % figfile
		    #plt.savefig(figfile)


    def calc_ogp_noise(self, raw):

        """
	   Determines OGP coefficients from snapshot of noise source
	"""
        N = len(raw)

        raw_off = sum(raw)/float(N)
	    
	raw_amp = sum(abs(raw-raw_off))/float(N)

	ogp_chan = []

	for i in range(0,2):

		core = raw[i::2]
		n = len(core)

		off = (sum(core)/n)*(-500.0/256.0)
		amp = sum(abs(core-off))/float(n)
		rel_amp = 100.0*(raw_amp-amp)/raw_amp
		
		ogp_core = (off, rel_amp, 0)
		ogp_chan.append(ogp_core)

	return ogp_chan


    def update_ogp(self, fname='ogp_chans01.npz'):
	    
        """
	   Updates ogp from .npz file
	"""

	df = np.load(fname)
	multi_ogp = df['multi_ogp']
	chans = df['chans']

	i = 0
	for chan in chans:

		ogp_chan = multi_ogp[i]

		
		print
		print "Setting ogp for chan %i..." %chan
		print ogp_chan
		print
	
		self.set_ogp(ogp_chan = ogp_chan, chan = chan)
		i += 1

	#self.roach.progdev(self.bitstream)


    def clear_ogp(self, chans = [0,1,2,3]):

        """
	   Sets OGP to 0 for channels in 'chans'
        """

	for chan in chans:

		zdok, cores = self.get_channel_core_spi(chan)

		for core in cores:

			adc5g.set_spi_offset(self.roach, zdok, core, 0)

			adc5g.set_spi_gain(self.roach, zdok, core, 0)

			adc5g.set_spi_phase(self.roach, zdok, core, 0)

        #self.roach.progdev(self.bitstream)

#
# INL Functions
#

    def get_code_errors(self, raw, freq):

         """
	    Determines residual errors for all 256 output codes
	 """

         s = []
	 c = []
	 x = []
	 del_phi = 2*math.pi*freq/self.clk

	 for i in range(len(raw)):

		 s += [math.sin(del_phi*i)]
		 c += [math.cos(del_phi*i)]
		 x += [del_phi*i]

	 s1 = s[0::2]
	 s2 = s[1::2]
	 c1 = c[0::2]
	 c2 = c[1::2]
	 x1 = x[0::2]
	 x2 = x[1::2]
         core1 = raw[0::2]
	 core2 = raw[1::2]

	 params1 = curve_fit(syn_func, x1, core1, p0=[120., 90., 90.]) 
	 params2 = curve_fit(syn_func, x2, core2, p0=[120., 90., 90.]) 

	 args1 = (np.array(s1), np.array(c1), np.array(core1)) 
	 args2 = (np.array(s2), np.array(c2), np.array(core2))

	 Fit1 = fitsin(params1[0], args1[0], args1[1])
	 Fit2 = fitsin(params2[0], args2[0], args2[1])

	 code_errors = np.zeros((256,2), dtype='float')
	 ce_counts = np.zeros((256,2), dtype='int32')

	 for i in range(len(core1)):

		 code1 = core1[i]
		 code2 = core2[i]

		 code_errors[code1][0] += code1 - Fit1[i]
		 code_errors[code2][1] += code2 - Fit2[i]

		 ce_counts[code1][0] += 1
		 ce_counts[code2][1] += 1

	 errors = np.zeros((256,2))

	 for code in range(256):

		 if ce_counts[code][0] != 0:

			 errors[code][0] = code_errors[code][0]/ce_counts[code][0]

		 else:

			 errors[code][0] = 0

		 if ce_counts[code][1] != 0:

			 errors[code][1] = code_errors[code][1]/ce_counts[code][1]
		 
		 else:

			 errors[code][1] = 0

	 return errors



    def fit_inl(self, raw, freq):

         """
	    Finds 17 INL code correction coefficients for each core in a channel
	    based on residual code errors
	 """

         errors1, errors2 = zip(*self.get_code_errors(raw, freq))
         
         corrections = np.zeros((17, 3), dtype = 'float')

         wts = array([1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,
                 15.,14.,13.,12.,11.,10.,9.,8.,7.,6.,5.,4.,3.,2.,1.])

	 start_data = 0
	 file_limit = len(errors1)
	 data_limit = start_data + file_limit

	 for corr_level in range(17):

		 a = corr_level*16 - 15 - start_data
		 b = a + 31

		 if a < 0:
			 a = 0

		 if a == 0 and start_data == 0:
			 a = 1

		 if b > file_limit:
			 b = file_limit

		 if b == file_limit and data_limit == 256:
			 b -= 1

		 if a > b:
			 continue

		 wt_a = a - corr_level*16 + 15 + start_data
		 wt_b = wt_a -a + b
		 wt = sum(wts[wt_a:wt_b])
		 av1 = sum(errors1[a:b]*wts[wt_a:wt_b])/wt
		 av2 = sum(errors2[a:b]*wts[wt_a:wt_b])/wt
		 corrections[corr_level][0] = 16*corr_level
		 corrections[corr_level][1] = av1
		 corrections[corr_level][2] = av2
        
	
         codes, inl1, inl2 = zip(*corrections)
	 inl_chan = (inl1, inl2)

	 return inl_chan


    def get_inl_chan(self, chan):

         """
	    Returns INL correction coefficients (17 coeff/core) 
	    for channel 'chan'
	 """

	 inl_chan = []

	 zdok, cores = self.get_channel_core_spi(chan)

	 for core in cores:

		 inl = adc5g.get_inl_registers(self.roach, zdok, core)
            
		 inl_chan.append(inl)

	 return inl_chan


    def get_inl(self, chans=[0,1]):
    
         """
	    Get INL for all channels in 'chans'
	 """
     
	 for chan in chans:

	        inl_chan = self.get_inl_chan(chan)

	        print
                print 'INL set on channel %i' %chan
                print

                print 'Core 1'
                print inl_chan[0]
                print

                print 'Core 2'
                print inl_chan[1]
                print   


    def set_inl(self, chan, inl_chan):

        """
	   Sets INL for two cores of channel 'chan'
	   inl_chan is format (inl1, inl2)
	"""

        zdok, cores = self.get_channel_core_spi(chan)

        i = 0

	for core in cores:

		inl1, inl2 = inl_chan

		adc5g.set_inl_registers(self.roach, zdok, core, inl1)
		adc5g.set_inl_registers(self.roach, zdok, core, inl2)
		
		i += 1
	
	#self.roach.progdev(self.bitstream)
        

    def do_inl(self, chans, freq, set_inl = False, save = False, 
	       fname = 'inl_chans01.npz'):

       """
          Calculates INL from snap shot for all channels in 'chans'
	  Can set INL and/or save to file
       """
       
       multi_inl = []

       for chan in chans:

                raw, f = self.get_snap(chan, freq)

		inl1, inl2 = self.fit_inl(raw, f)

                inl_chan = (inl1, inl2)

                if set_inl:

                        self.set_inl(chan, inl_chan)

                multi_inl.append(inl_chan)

       if save:
	       np.savez(fname, multi_inl = multi_inl, chans = chans)

       return multi_inl, chans


    def update_inl(self, fname='inl_chans01.npz'):

        """
	   Updates INL from .npz file
	"""

	df = np.load(fname)
	multi_inl = df['multi_inl']
	chans = df['chans']

	i = 0
        for chan in chans:

                inl_chan = multi_inl[i]

                print
                print "Setting inl for chan %i..." %chan
                print inl_chan
                print

                self.set_inl(inl_chan = inl_chan, chan = chan)
                i += 1

        #self.roach.progdev(self.bitstream)
		


    def clear_inl(self, chans = [0,1,2,3]):

        """
	   Sets INL to 0 for all channels in 'chans'
	"""

        inl0 = np.zeros(17)

        for chan in chans:

              zdok, cores = self.get_channel_core_spi(chan)

	      for core in cores:

		      adc5g.set_inl_registers(self.roach, zdok, core, inl0) 

        #self.roach.progdev(self.bitstream)

#
# Misc. 
#        
    
    #def save_raw_npz(self, raw, freq, chan):

	#t = datetime.datetime.now()
        #header = t.strftime("%Y-%m-%d %H:%M") + '_'
	#snap_file = '%s_%i_MHZ_raw_%s_snap.npz' %(header, freq, chan)
	#np.savez(snap_file, header=header, freq=freq, raw=raw)
	    

    def levels_hist(self, raw):

       	bins = np.arange(-128.5, 128.5, 1)
	
	th = 32
	lim1 = -th-0.5
	lim2 = -0.5
	lim3 = th-0.5

	ideal_gauss = mlab.normpdf(bins, 0, th)

	plt.subplot(111)
	plt.plot((lim1,lim1), (0,1), 'k--')
	plt.plot((lim2,lim2), (0,1), 'k--')
	plt.plot((lim3,lim3), (0,1), 'k--')
	plt.plot(bins, ideal_gauss, 'gray', linewidth=1)

	plt.hist(raw, bins, normed=1, facecolor='blue', 
		 alpha=0.9, histtype='stepfilled')

	plt.xlim(-129, 128)
	plt.ylim(0, 0.06)
	plt.xlabel('ADC Value')
	plt.ylabel('Normalized Count')
	plt.title('zdok0 levels')
	plt.grid()
	
	plt.show()

    def plot_compare_ogp(self, new_ogp, old_ogp, freq, pts=50):

	    """
	    This function plots a comparison between a snap and a fitted sine
	    for two different sets of ogp
  
	    Inputs: 
	    - 'new_ogp' is an array containing measured ogp values,
	    - 'old_ogp' is an array containing current ogp values loaded on ADC5g
	    (if old_ogp = None, will use 0 ogp for all 4 cores),
	    - 'freq' is frequency of input CW tone (in MHz)
	    - 'pts' is for how many points in snap to plot (max 8192)

	    Outputs:
	    -
	    -
	    -
	    """
	    tot_pts = 8192
	    del_phi = 2*math.pi*freq/self.clk #freq in MHz
	    x = [del_phi*i for i in range(tot_pts)]
	    r = np.arange(0, tot_pts*del_phi, del_phi/10.)
	    p0 = [128., 90., 90.]

	    self.syn.set_freq(freq*1e6)

	    self.clear_ogp()
	    time.sleep(0.5)
	    raw_u = adc5g.get_snapshot(self.roach, 'scope_raw_a0_snap')
	    params_u = curve_fit(syn_func, x, raw_u, p0)

	    self.set_ogp(new_ogp, 0, cores=[1,2])
	    time.sleep(0.5)
	    raw_c = adc5g.get_snapshot(self.roach, 'scope_raw_a0_snap')
	    params_c = curve_fit(syn_func, x, raw_c, p0)

	    f_u =  syn_func(r, params_u[0][0], params_u[0][1], params_u[0][2])
	    f_c =  syn_func(r, params_c[0][0], params_c[0][1], params_c[0][2])

	    fig = plt.figure()
	    ax1 = fig.add_subplot(211)
	    ax2 = fig.add_subplot(212)	    
	    
	    ax1.plot(r[:(pts*10)], f_u[:(pts*10)])
	    ax1.plot(x[:pts], raw_u[:pts], 'o', color='r')
	    ax1.set_title('OGP Uncorrected Fit, freq = %d MHz' %freq)

	    ax2.plot(r[:(pts*10)], f_c[:(pts*10)])
            ax2.plot(x[:pts], raw_c[:pts], 'o', color='r')
            ax2.set_title('OGP Corrected Fit, freq = %d MHz' %freq)

	    plt.tight_layout()
	    plt.show()

	    #return 
	    
    
#if __name__ == '__main__':

    # Program ROACH in ipython shell...
    #
    #roach = katcp_wrapper.FpgaClient('172,30.51.97')
    #roach.progdev('adc5g_tim_aleks_test_2015_Oct_14_1208.bof.gz')
