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

#
# Functions used by the fitting routines
#

def syn_func( x, off, a1, a2):
        return off +  a1*np.sin(x) + a2*np.cos(x)

def syn_func_nonlin( x, off, a1, p1):
        return off +  a1*np.sin(x+p1) 

def fitsin(p, sin, cos):
	return p[0] + p[1]*sin + p[2]*cos

def sin_residuals(p, sin, cos, raw):
        res = raw - fitsin(p, sin, cos)
        for i in range(raw.size):
            if raw[i] == -128 or raw[i] == 127:
                res[i] = 0
        return res

#
# ADC Calibration Tools
#
# 
#
#

class ADC5g_Calibration_Tools (object):
    
    def __init__(self, roach, roach_id = '172.30.51.97', clk=1600):

        self.roach_id = roach_id
        self.bitstream = 'adc5g_tim_aleks_test_2015_Oct_14_1208.bof.gz'
        self.clk = clk
	self.roach = roach
	self.freqs = np.arange(6)*50 + 100.0
	self.syn = HP8780A()
	#self.roach = katcp_wrapper.FpgaClient(roach_id)
	#self.roach.progdev(self.bitstream)

    #	
    # Takes a raw snap shot of time domain signal, separates both cores, and finds fitted parameters
    # A,B,offset to linear sine function (y = A*cos(x) + B*sin(x) + offset)
    #
    def fit_snap(self, freq, raw): 

        del_phi = 2*math.pi*freq/self.clk

        x  = [del_phi*i for i in range(len(raw))]

        core1 = raw[0::2]
        core2 = raw[1::2]

        x1 = x[0::2]
        x2 = x[1::2]
        
	p0 = [128.,90.,90.]
        
        params1 = curve_fit(syn_func, x1, core1, p0)

        params2 = curve_fit(syn_func, x2, core2, p0)
    
        return params1,params2
    
    #
    # Takes fitted parameters from a linear sine model and determines
    # values for offset, gain and phase corrections on ADC
    #
    def convert_fit_to_ogp(self, freq, params1, params2):

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

        return (ogp1, ogp2)

    #
    # Returns the current ogp values from an ADC (zdok = 0 or 1) for a list
    # of cores (cores = [1,2,3,4])
    #
    def get_ogp(self, zdok = 0, cores = [1,2,3,4]):
        
        multi_ogp = []

        for core in cores:

            off = adc5g.get_spi_offset(self.roach, zdok, core) 
            gain = adc5g.get_spi_gain(self.roach, zdok, core)
            phase = adc5g.get_spi_phase(self.roach, zdok, core)
            ogp_core = [off, gain, phase]
            multi_ogp.append(ogp_core)

        return multi_ogp

    #
    # Sets ogp values contained in 'multi_ogp' array for each core in 'cores'
    # on specified ADC ('zdok')
    #
    def set_ogp(self, multi_ogp, zdok, cores=[1,2,3,4]):

        adc5g.set_spi_control(self.roach, zdok)
        
        for core in cores:
            
            off, gain, phase = multi_ogp[core-1]

            off_spi = math.floor(.5 + off*255/100.) + 0x80
            adc5g.set_spi_offset(self.roach, zdok, core, float(off))

            gain_spi = math.floor(.5 + gain*255/36.) + 0x80
            adc5g.set_spi_gain(self.roach, zdok, core, float(gain))

            phase_spi = math.floor(.5 + phase*255/28.) + 0x80
            adc5g.set_spi_phase(self.roach, zdok, core, float(phase)*0.65)

	self.roach.progdev(self.bitstream)

    #
    # Generates frequency averaged ogp values derived from sine fit for each core
    # using swept input CW tone
    #
    def do_ogp_cw_sweep(self, zdoks=[0], save=False, fname='ogp_default.npz', save_raws=False, 
			raw_fname = 'raw_snaps.npz', raw_len=8192):

	freqs = self.freqs

	zdok = 0

	ogpA = {}
	ogpB = {}
	ogpC = {}
	ogpD = {}

	rawsAB = np.zeros((len(freqs),raw_len))
	rawsCD = np.zeros((len(freqs),raw_len))

	

	for i in range(len(freqs)): # MHz

		freq = freqs[i]
		self.syn.set_freq(freq*1e6)
		time.sleep(2.0)
	
		rawAB = adc5g.get_snapshot(self.roach, 'scope_raw_a%i_snap' %(zdok))
		rawCD = adc5g.get_snapshot(self.roach, 'scope_raw_c%i_snap' %(zdok))  

		rawsAB[i] = rawAB
		rawsCD[i] = rawCD

		paramsA, paramsB = self.fit_snap(freq, rawAB)
		paramsC, paramsD = self.fit_snap(freq, rawCD)

		ogpA[freq], ogpB[freq] = self.convert_fit_to_ogp(freq, 
								 paramsA, paramsB)
		ogpC[freq], ogpD[freq] = self.convert_fit_to_ogp(freq, 
								 paramsC, paramsD)

	ogpA_m = tuple(map(lambda y: sum(y)/float(len(y)), zip(*ogpA.values())))
        ogpB_m = tuple(map(lambda y: sum(y)/float(len(y)), zip(*ogpB.values())))
        ogpC_m = tuple(map(lambda y: sum(y)/float(len(y)), zip(*ogpC.values())))
        ogpD_m = tuple(map(lambda y: sum(y)/float(len(y)), zip(*ogpD.values())))

        multi_ogp = (ogpA_m,ogpB_m,ogpC_m,ogpD_m)

	if save:
		
		np.savez(fname, zdok0_ogp = multi_ogp)

        if save_raws:
		
		np.savez(raw_fname, freqs = freqs, rawsAB = rawsAB, rawsCD = rawsCD)

	return multi_ogp


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
    def get_resid(self, freq, raw):
    
        params1,params2 = self.fit_snap(freq, raw)

        core1 = raw[0::2]
        core2 = raw[1::2]

	del_phi = 2*np.pi*freq/self.clk

	x = [del_phi*i for i in range(len(raw))]

        x1 = x[0::2]
        x2 = x[1::2]

	y1 =  syn_func(x1, params1[0][0], params1[0][1], params1[0][2])
	y2 =  syn_func(x2, params2[0][0], params2[0][1], params2[0][2])
 
        resid1 = y1 - core1
        resid2 = y2 - core2
        
        return resid1, resid2

    #
    # Determines INL corrections based on residual data
    #
    def convert_residuals_to_inl(self, freq, raw):
         
         resid1, resid2 = self.get_resid(freq,raw)
         
         corrections = np.zeros((17,3), dtype = 'float')
         wts = array([1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,
                 15.,14.,13.,12.,11.,10.,9.,8.,7.,6.,5.,4.,3.,2.,1.])
	 start_data = int(resid1[0])
	 file_limit = len(resid1)
	 data_limit = start_data+file_limit

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
		 av1 = sum(resid1[a:b]*wts[wt_a:wt_b])/wt
		 av2 = sum(resid2[a:b]*wts[wt_a:wt_b])/wt
		 corrections[corr_level][0] = 16*corr_level
		 corrections[corr_level][1] = av1
		 corrections[corr_level][2] = av2
        
	 return corrections
        
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
    
 
#    def set_inl(self, fname = 'inl_default.npz'):

#        zdok = 0        
        
#        df = np.load(fname)
#	zdok0_inl = df['zdok0_inl']
 
#        adc5g.set_inl_registers(self.roach,zdok,1,zdok0_inl[1])
#        adc5g.set_inl_registers(self.roach,zdok,2,zdok0_inl[2])
#        adc5g.set_inl_registers(self.roach,zdok,3,zdok0_inl[3])
#        adc5g.set_inl_registers(self.roach,zdok,4,zdok0_inl[4])
             
        
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

              adc5g.set_inl(self.roach, zdok, core, 0)

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




    









    




	






















