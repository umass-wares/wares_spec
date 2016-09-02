import numpy as np
from numpy import vstack
from matplotlib import pyplot as plt

Nout = 4000
maxk = Nout/4
sampleint = 1.0

lochan = 1500
hichan = 1700 

BW = 800.e6/(2048.)

data = np.load('allan_spectra_800mhz_4000samps_1int_aug24_1.npz')['spec_matrix']
#start_t = np.load('allan_800mhz_nocal_really.npz')['times']
#read_t = np.load('a;l')['read_t']

alanvar = []

for k in range(1,maxk+1):

    M = Nout/k

    for n in range(M):

        tmpR = data[(n*k):(n+1)*k, lochan:hichan].sum(axis=0)/k
        if n == 0:
            R = tmpR
        else:
            R = vstack((tmpR, R))

    sigmak = 0.0
    avgspec = np.zeros((hichan-lochan,), dtype='float')

    for n in range(M-1):
            
        spec = (R[n+1,:] - R[n,:])/R[n,:]
        avgspec += spec
        sigmak += spec.var()

    avgspec = avgspec*np.sqrt(k)

    if k == 1:
        avgspectra = avgspec
    else:
        avgspectra = vstack((avgspec,avgspectra))

    sigmak = sigmak/(M-1)

    alanvar.append(np.sqrt(sigmak))

alanvar = np.array(alanvar)
time = sampleint+np.arange(maxk)*sampleint

theory = np.sqrt(2.)/(np.sqrt(BW)*np.sqrt(time))
#theory = 1./np.sqrt(time)

#plt.plot(time, alanvar, marker='o')
plt.plot(time, alanvar, label='alanvar')
plt.plot(time, theory, label='theory')
plt.xscale('log')
plt.yscale('log')
#plt.ylim(1.0, 3.0)
plt.xlabel('Integration Time')
plt.ylabel('$\sigma_{Allan}$')
plt.title('Allan Plot, 800 MHz Spectrometer, 1s Int, Chans 1300-1500, Calibrated')
plt.legend(loc='best')
plt.show()
