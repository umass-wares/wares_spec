from corr import katcp_wrapper
import adc5g
from hp8780a import HP8780A
import time
import datetime
import numpy as np
from numpy import array, mean
from scipy.optimize import leastsq
from scipy.optimize import curve_fit
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

class ADC5g_Calibration_Tools (object):
    
    def __init__(self, roach_id = '172.30.51.97', 
		 bitstream = 'adc5g_tim_aleks_test_2015_Oct_14_1208.bof.gz',
		 clk=1600):

	    self.roach_id = roach_id
	    self.bitstream = bitstream
	    self.clk = clk
	    self.freqs = np.arange(6)*50 + 100.0
	    self.syn = HP8780A()

	    roach = katcp_wrapper.FpgaClient(roach_id)
	    roach.wait_connected()
	    roach.progdev(self.bitstream)
	    self.roach = roach

    def get_channel_snap_reg(chan):
	    """
	       Gets the snap register string based on channel No. (0,1,2,3)
	    """
	    reg = None

	    if chan == 0:
		    reg = 'scope_raw_a0_snap_bram'

	    if chan == 1:
		    reg = 'scope_raw_c0_snap_bram'

	    if chan == 2:
		    reg = 'scope_raw_a1_snap_bram'
	    
	    if chan == 3:
		    reg = 'scope_raw_c1_snap_bram'

	    return reg
     
    
    def get_channel_core_spi(chan):

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
		    
    def get_snap(chan, freq): #MHz
	    """
	       Takes a snap shot on channel 'chan' of a cw tone
	       of frequency 'freq' (MHz)
	    """
	    
	    self.syn.output_on()
	    self.syn.set_freq(freq*1e6)
	    time.sleep(.5)

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
	    Calculates the SINAD from the residuals of the fit. 
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
        
    def do_sfdr(self, raw, freq, nfft=1024):

	    """
	    Calculates the sinad and sfdr of the raw data. It first gets the PSD.
	    Then, finds the fundamental peak and the maximum spurious peak. 
	    (Harmonics are also spurs). The DC level is excluded.
	    """

	    samp_freq=self.clk
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
        
	    for i in range(4,len(freqs)):

		    if abs(freqs[i] - sig_freq) < 4:
			    test = -70
		    else:
			    test = -90

		    pwr = 10**(float(db[i])/10.)
		    tot_pwr += pwr

		    if in_peak:

			    if db[i] < test:
				    in_peak = False

				    if abs(peak_freq - sig_freq) < 1:
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

	    sfdr = 10.0*math.log10(sig_pwr / spur_pwr)
	    sinad = 10.0*math.log10(sig_pwr/(tot_pwr - sig_pwr))
	    return sfdr, sinad
        
    def do_sfdr_sinad_cw_sweep( self, zdok=0, save=False, 
		       fname='sfdr_sinad.npz'):

        """
        Calculates the SFDR and SINAD from a sweep in frequency from 50 Mhz to
        final_freq Mhz. 
        """

	final_freq = self.clk /2
	freqs = np.linspace(50,final_freq,16)
 
     #When freq is 400Mhz, it has no spurrious frequency with clk=1600. SFDR 
     #can't be calculated in this case. 
	if 400 in freqs:
	     x=np.where(freqs==400)
	     freqs[x[0][0]]=401

	sfdrAB = {}
	sfdrCD = {}
	sinadAB_psd= {}
	sinadCD_psd= {}

	for i in range(len(freqs)): # MHz

		freq = freqs[i]
		self.syn.set_freq(freq*1e6)
		time.sleep(2.0)
	
		rawAB = adc5g.get_snapshot(self.roach, 'scope_raw_a%i_snap' %(zdok))
		rawCD = adc5g.get_snapshot(self.roach, 'scope_raw_c%i_snap' %(zdok))
  
		sfdrAB[freq], sinadAB_psd[freq] = self.do_sfdr(freq,rawAB)
		sfdrCD[freq], sinadCD_psd[freq] = self.do_sfdr(freq,rawCD)
  
  		multi_sfdr = (sfdrAB, sfdrCD)
		multi_sinad_psd = (sinadAB_psd, sinadCD_psd)
        
	if save:		
		np.savez(fname, zdok_sfdr=multi_sfdr, zdok_sinad_psd=multi_sinad_psd)

	return multi_sfdr, multi_sinad_psd


    def do_sinad_cw_sweep( self, zdok=0, save=False, fname='sinad.npz'):

        """
        Calculates the SINAD from a sweep in frequency from 50 Mhz to
        final_freq Mhz. 
        """

        final_freq = self.clk /2
        freqs = np.linspace(50,final_freq,16)     
        sinadAB= {}
        sinadCD= {}
    
        for i in range(len(freqs)): # MHz
    
    		freq = freqs[i]
    		self.syn.set_freq(freq*1e6)
    		time.sleep(2.0)
    	
    		rawAB = adc5g.get_snapshot(self.roach, 'scope_raw_a%i_snap' %(zdok))
    		rawCD = adc5g.get_snapshot(self.roach, 'scope_raw_c%i_snap' %(zdok))
      
    		sinadAB[freq] = self.calc_sinad(freq,rawAB)
    		sinadCD[freq] = self.calc_sinad(freq,rawCD)
      
    		multi_sinad = (sinadAB, sinadCD)
            
        if save:		
    		np.savez(fname, zdok_sinad=multi_sinad)
    
        return multi_sinad     


    def get_ogp(self, chan):

        """
	   Returns ogp values for a channel 'chan'
	"""

        multi_ogp = []
    
        zdok, cores = self.get_channel_cores_spi(chan)

        for core in cores:

            off = adc5g.get_spi_offset(self.roach, zdok, core) 
            gain = adc5g.get_spi_gain(self.roach, zdok, core)
            phase = adc5g.get_spi_phase(self.roach, zdok, core)
            ogp_core = [off, gain, phase]

            multi_ogp.append(ogp_core)

        return multi_ogp

    def set_ogp(self, multi_ogp, chan):

        """
	   Sets ogp for two cores of channel 'chan'
	   multi_ogp is format (ogp1, ogp2)
        """

        zdok, cores = self.get_channel_cores_spi(chan)

        for core in cores:
            
            off, gain, phase = multi_ogp[core-1]

            off_spi = math.floor(.5 + off*255/100.) + 0x80
            adc5g.set_spi_offset(self.roach, zdok, core, float(off))

            gain_spi = math.floor(.5 + gain*255/36.) + 0x80
            adc5g.set_spi_gain(self.roach, zdok, core, float(gain))

            phase_spi = math.floor(.5 + phase*255/28.) + 0x80
            adc5g.set_spi_phase(self.roach, zdok, core, float(phase)*0.65)

        self.roach.progdev(self.bitstream)


    def do_ogp_cw_sweep(self, chan, set=True, 
			save=False, fname='ogp_default.npz', 
			save_raws=False, raw_fname = 'raw_snaps.npz', 
			raw_len=8192, sinad=True, sfdr=False):

	freqs = self.freqs

	ogpA = {}
	ogpB = {}
	ogpC = {}
	ogpD = {}

	for i in range(len(freqs)): # MHz

		freq = freqs[i]
 
		raw, f = self.get_snap(chan, freq) 
		
		#rawsAB[i] = rawAB
		#rawsCD[i] = rawCD

		params1, params2 = self.fit_snap(freq, raw)
  
		ogp1[freq], ogp2[freq] = self.convert_fit_to_ogp(freq, 
								 paramsA, paramsB)
		
        ogpA_m = tuple(map(lambda y: sum(y)/float(len(y)), zip(*ogpA.values())))
        ogpB_m = tuple(map(lambda y: sum(y)/float(len(y)), zip(*ogpB.values())))
        ogpC_m = tuple(map(lambda y: sum(y)/float(len(y)), zip(*ogpC.values())))
        ogpD_m = tuple(map(lambda y: sum(y)/float(len(y)), zip(*ogpD.values())))
        
        multi_ogp = (ogpA_m,ogpB_m,ogpC_m,ogpD_m)        
        
        if save:		
		np.savez(fname, zdok0_ogp = multi_ogp)       

        if save_raws:		
		np.savez(raw_fname, freqs = freqs, rawsAB = rawsAB, rawsCD = rawsCD)

	return multi_ogp, multi_sinad, multi_sfdr


    #
    # Generates average ogp values from noise source
    #
    def do_ogp_noise_source(self, rpt, save=False, fname='ogp_noise_default.npz'):

        ogpA = {}
        ogpB = {}
        ogpC = {}
        ogpD = {}

        for i in range(rpt):

		rawAB = np.array(adc5g.get_snapshot(self.roach, 'scope_raw_a0_snap')) 
		rawCD = np.array(adc5g.get_snapshot(self.roach, 'scope_raw_c0_snap'))

		ogpA[i], ogpB[i] = self.convert_noise_to_ogp(rawAB)
		ogpC[i], ogpD[i] = self.convert_noise_to_ogp(rawCD)

        ogpA_m = tuple(map(lambda y: sum(y)/float(len(y)), zip(*ogpA.values())))
        ogpB_m = tuple(map(lambda y: sum(y)/float(len(y)), zip(*ogpB.values())))
        ogpC_m = tuple(map(lambda y: sum(y)/float(len(y)), zip(*ogpC.values())))
        ogpD_m = tuple(map(lambda y: sum(y)/float(len(y)), zip(*ogpD.values())))

        multi_ogp = (ogpA_m,ogpB_m,ogpC_m,ogpD_m)

        if save:
                np.savez(fname, zdok0_ogp = multi_ogp)

        return multi_ogp
	
    #
    # Takes noise source measurement and determines values for ogp corrections 
    # on ADC
    #
    def convert_noise_to_ogp(self, raw):

        N = len(raw)

        raw_off = sum(raw)/float(N)
	    
	raw_amp = sum(abs(raw-raw_off))/float(N)

	multi_ogp = []

	for i in range(0,2):

		core = raw[i::2]
		n = len(core)

		off = (sum(core)/n)*(-500.0/256.0)
		amp = sum(abs(core-off))/float(n)
		rel_amp = 100.0*(raw_amp-amp)/raw_amp
		
		ogp_core = (off, rel_amp, 0)
		multi_ogp.append(ogp_core)

	return multi_ogp

    #
    # Updates ogp from file
    #
    def update_ogp(self, fname='ogp_noise_default.npz'):

        zdok = 0

	df = np.load(fname)
	zdok0_ogp = df['zdok0_ogp']

	print
	print "Setting ogp for zdok0..."
	print zdok0_ogp
	print
	
	self.set_ogp(multi_ogp = zdok0_ogp, zdok = zdok)

    #
    # Sets ogp to 0 for a list of cores
    #
    def clear_ogp(self, cores=[1,2,3,4]):

	for core in cores:

              adc5g.set_spi_offset(self.roach, 0, core, 0)
              adc5g.set_spi_gain(self.roach, 0, core, 0)
	      adc5g.set_spi_phase(self.roach, 0, core, 0)

	self.roach.progdev(self.bitstream)
	
    #
    # Determines residuals from raw snapshot and fitted sine
    #
    def get_code_errors(self, freq, raw):

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

    #
    # Determines INL corrections based on code error data
    #
    def fit_inl(self, freq, raw):
         
         errors1, errors2 = zip(*self.get_code_errors(freq,raw))
         
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
        
	 return corrections


    def get_inl(self, zdok = 0, cores = [1,2,3,4]):
	    """
	        Get INL correction coefficients (17 coeff/core)
	    """

	    multi_inl = []

	    for core in cores:

		    inl = adc5g.get_inl_registers(self.roach, zdok, core)
            
		    multi_inl.append(inl)

	    return multi_inl


    def do_inl_sweep(self, freq, zdok = [0], save = False, fname = 'inl_default.npz'):
    	
    	inlA = {}
    	inlB = {}
    	inlC = {}
    	inlD = {}
    	
    	rawAB = adc5g.get_snapshot(self.roach, 'scope_raw_a%i_snap' %(zdok))
	rawCD = adc5g.get_snapshot(self.roach, 'scope_raw_c%i_snap' %(zdok))
	
	residA, residB = self.get_resid(freq, rawAB)
	residC, residD = self.get_resid(freq, rawCD)
	
	correctionsAB = self.convert_residuals_to_inl(freq, rawAB)
	correctionsCD = self.convert_residuals_to_inl(freq, rawCD)
	
	code_num,inlA,inlB = zip(*correctionsAB)
	code_num,inlC,inlD = zip(*correctionsCD)
	
	multi_inl = (code_num, inlA, inlB, inlC, inlD)
	
	if save:
		np.savez(fname, zdok0_inl = multi_inl)
		
	return multi_inl
    
 
    def set_inl(self, errors):

        zdok = 0        
        
        #df = np.load(fname)
	#zdok0_inl = df['zdok0_inl']
 
	codes, inl_a, inl_b = zip(*errors)

        adc5g.set_inl_registers(self.roach, zdok, 1, inl_a)
        adc5g.set_inl_registers(self.roach, zdok, 2, inl_b)
        #adc5g.set_inl_registers(self.roach,zdok,3,zdok0_inl[3])
        #adc5g.set_inl_registers(self.roach,zdok,4,zdok0_inl[4])
             
        
    #
    # Sets INL corrections from saved file
    #
    def update_inl(self, fname='inl_default.npz'):

        zdok = 0

	df = np.load(fname)
	zdok0_inl = df['zdok0_inl']

	print
	print "Setting inl for zdok0..."
	print zdok0_inl
	print
	
	self.set_inl(multiple_inl = zdok0_inl, zdok = zdok)
 
    #
    # Sets all INL to 0
    #
    def clear_inl(self):
        
        zdok = 0

        for core in range(1,5):

              adc5g.set_inl_registers(self.roach, zdok, core, np.zeros(17))

        self.roach.progdev(self.bitstream)
        
    
    def save_raw_npz(self, raw, freq, chan):

	t = datetime.datetime.now()
        header = t.strftime("%Y-%m-%d %H:%M") + '_'
	snap_file = '%s_%i_MHZ_raw_%s_snap.npz' %(header, freq, chan)
	np.savez(snap_file, header=header, freq=freq, raw=raw)
	    

    def levels_hist(self, raw):

        raw = np.array(adc5g.get_snapshot(self.roach, 'scope_raw_a0_snap'),
		       np.int32)

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

#	    if (curr_ogp == None):
#		    self.clear_ogp()
#		    time.sleep(0.5)
#	    else:
#		    self.set_ogp(old_ogp, 0, cores=[1,2])

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
