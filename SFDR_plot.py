# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 13:21:58 2016

@author: Sandra Bustamante
"""
import numpy as np
import matplotlib as plt

def dic2arr(dic,chans=[0,1]):
    """
    dic should be a tuple of 2 dictionaries. Assumes keys are the same.
    
    """
    dic_keys=dic[0].keys()
    dic_keys.sort() #es una lista con las frecuencia de menor a mayor
    y=len(dic_keys)
    x=len(chans)
    dic_values=np.zeros((x,y))
    for chan in chans:
        for i in range(x):
            dic_key=dic_keys[i]
            dic_values[chan][i]=dic[chan][dic_key]
        
    return dic_values, dic_keys
    
def plot_sinad_sfdr (label, data_x, data_y, chans=[0,1,2,3],
                     titles=['SFDR','SINAD']):
    """
    x   x values of data (same for all chans)
    y   array with shape (2, chans, data)
    """
    n=len(chans)    
    n2=len(titles)
    pos=np.arange(n2*n)+1

    for t in range(n2):
        pos_val=pos[t::n2]
        for chan in chans:
            plt.subplot(n,n2,pos_val[chan])
            plt.plot(data_x,data_y[t][chan],label=label)
            if t==0:
                plt.ylabel('Chan %i' %chan)
            if chan==0:
                plt.title(titles[t])

#from corr import katcp_wrapper
#roach = katcp_wrapper.FpgaClient('172.30.51.97')
#roach.progdev('adc5g_tim_aleks_test_2015_Oct_14_1208.bof.gz')

from adc5g_cal import ADC5g_Calibration_Tools
adc_cal = ADC5g_Calibration_Tools()

freqarray=[50,800,31]
chans=[0,1]

""" No corrections """

adc_cal.clear_ogp()
adc_cal.clear_inl()

sfdr,sinad=adc_cal.do_sfdr_sinad_cw_sweep(chans=chans, freqarray=freqarray)

sinad_values,freqs= dic2arr(sinad)
sfdr_values,f = dic2arr(sfdr)
data_NC=(sfdr_values,sinad_values)

""" OGP tone corrections """
adc_cal.do_ogp_cw_sweep(chans=chans)

multi_ogp=[]
for chan in chans:
    ogp = adc_cal.get_ogp_chan(chan)
    multi_ogp.append(ogp)

sfdr,sinad = adc_cal.do_sfdr_sinad_cw_sweep(freqarray=freqarray)

sinad_values,freqs= dic2arr(sinad)
sfdr_values,f = dic2arr(sfdr)
data_OGPT=(sfdr_values,sinad_values)

""" OGP noise corrections """

label_5='OGP_noise corrections'
adc_cal.do_ogp_noise_sweep(chans=chans)
multi_ogpNoise=[]
for chan in chans:
    ogp = adc_cal.get_ogp_chan(chan)
    multi_ogpNoise.append(ogp)

sfdr,sinad = adc_cal.do_sfdr_sinad_cw_sweep(freqarray=freqarray)

sinad_values,freqs= dic2arr(sinad)
sfdr_values,f = dic2arr(sfdr)

data_5=(sfdr_values,sinad_values)

plot_sinad_sfdr (label=label_5, data_x=freqs, data_y=data_5, chans=[0,1],
                     titles=['SFDR','SINAD'])


#PLOTS
""" Plot No corrections"""
label_NC='No Corrections'
plot_sinad_sfdr (label=label_NC, data_x=freqs, data_y=data_NC, chans=chans,
                     titles=['SFDR','SINAD'])   

""" Plot INL Corrections """


""" Plot OGP tone corrections"""
label_OGPT='OGP_tone corrections'
plot_sinad_sfdr (label=label_OGPT, data_x=freqs, data_y=data_OGPT, chans=chans,
                     titles=['SFDR','SINAD'])



""" INL corrections """
label_3='INL corrections'
adc_cal.clear_ogp()
adc_cal.do_inl(chans=[0,1],freq=10,set_inl=True)

sfdr,sinad = adc_cal.do_sfdr_sinad_cw_sweep(freqarray=freqarray)

sinad_values,freqs= dic2arr(sinad)
sfdr_values,f = dic2arr(sfdr)

data_3=(sfdr_values,sinad_values)

plot_sinad_sfdr (label=label_3, data_x=freqs, data_y=data_3, chans=[0,1],
                     titles=['SFDR','SINAD'])


""" INL & OGP tone corrections """

label_4='INL & OGP tone corrections'
adc_cal.set_ogp(ogp0,0)
adc_cal.set_ogp(ogp1,1)

sfdr,sinad = adc_cal.do_sfdr_sinad_cw_sweep(freqarray=freqarray)

sinad_values,freqs= dic2arr(sinad)
sfdr_values,f = dic2arr(sfdr)

data_4=(sfdr_values,sinad_values)

plot_sinad_sfdr (label=label_4, data_x=freqs, data_y=data_4, chans=[0,1],
                     titles=['SFDR','SINAD'])



plt.legend(loc=0)

savefig('/home/sandra/wares_spec/Images/SfdrSinad031_OGPnoise.png')
savefig('/home/sandra/wares_spec/Images/SfdrSinad031_OGPnoise.eps')


"""
For channels 2 and 3
"""
adc_cal.clear_ogp()
adc_cal.clear_inl()

freqarray=[50,800,30]

sfdr,sinad=adc_cal.do_sfdr_sinad_cw_sweep(chans=[2,3], freqarray=freqarray)

sinad_values,freqs= dic2arr(sinad)
sfdr_values,f = dic2arr(sfdr)

""" No corrections """
plt.subplot(221)
plt.plot(freqs,sinad_values[0],'-*', label='No corrections')
plt.title('SINAD')
plt.ylabel('Chan 2')

plt.subplot(223)
plt.plot(freqs,sinad_values[1],'-*', label='No corrections')
plt.ylabel('Chan 3')

plt.subplot(222)
plt.plot(freqs,sfdr_values[0],'-*', label='No corrections')
plt.title('SFDR')

plt.subplot(224)
plt.plot(freqs,sfdr_values[1],'-*', label='No corrections')

""" OGP corrections """
adc_cal.do_ogp_cw_sweep(chans=[2,3])
ogp2 = adc_cal.get_ogp_chan(2)
ogp3 = adc_cal.get_ogp_chan(3)

sfdr,sinad = adc_cal.do_sfdr_sinad_cw_sweep(chans=[2,3], freqarray=freqarray)

sinad_values,freqs= dic2arr(sinad)
sfdr_values,f = dic2arr(sfdr)

plt.subplot(221)
plt.plot(freqs,sinad_values[0],'-*', label='OGP corrections')

plt.subplot(223)
plt.plot(freqs,sinad_values[1],'-*', label='OGP corrections')

plt.subplot(222)
plt.plot(freqs,sfdr_values[0],'-*', label='OGP corrections')

plt.subplot(224)
plt.plot(freqs,sfdr_values[1],'-*', label='OGP corrections')   

""" INL corrections """
adc_cal.clear_ogp()
adc_cal.do_inl(chans=[2,3],freq=10,set_inl=True)

sfdr,sinad = adc_cal.do_sfdr_sinad_cw_sweep(chans=[2,3], freqarray=freqarray)

sinad_values,freqs= dic2arr(sinad)
sfdr_values,f = dic2arr(sfdr)

plt.subplot(221)
plt.plot(freqs,sinad_values[0],'-*', label='INL corrections')

plt.subplot(223)
plt.plot(freqs,sinad_values[1],'-*', label='INL corrections')

plt.subplot(222)
plt.plot(freqs,sfdr_values[0],'-*', label='INL corrections')

plt.subplot(224)
plt.plot(freqs,sfdr_values[1],'-*', label='INL corrections')   

""" INL & OGP corrections """
adc_cal.set_ogp(ogp2,2)
adc_cal.set_ogp(ogp3,3)

sfdr,sinad = adc_cal.do_sfdr_sinad_cw_sweep(chans=[2,3], freqarray=freqarray)

sinad_values,freqs= dic2arr(sinad)
sfdr_values,f = dic2arr(sfdr)

plt.subplot(221)
plt.plot(freqs,sinad_values[0],'-*', label='INL & OGP corrections')

plt.subplot(223)
plt.plot(freqs,sinad_values[1],'-*', label='INL & OGP corrections')

plt.subplot(222)
plt.plot(freqs,sfdr_values[0],'-*', label='INL & OGP corrections')

plt.subplot(224)
plt.plot(freqs,sfdr_values[1],'-*', label='INL & OGP corrections')

plt.legend(loc=0)

savefig('/home/sandra/wares_spec/Images/SfdrSinad031_OGP23.png')
savefig('/home/sandra/wares_spec/Images/SfdrSinad031_OGP23.eps')



""" OGP Noise corrections """
plt.subplot(221)
plt.plot(freqs,sinad_values[0],'-*', label='OGP Noise corrections')

plt.subplot(223)
plt.plot(freqs,sinad_values[1],'-*', label='OGP Noise corrections')

plt.subplot(222)
plt.plot(freqs,sfdr_values[0],'-*', label='OGP Noise corrections')

plt.subplot(224)
plt.plot(freqs,sfdr_values[1],'-*', label='OGP Noise corrections')

plt.legend(loc=0)