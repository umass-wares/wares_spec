import numpy as np
from numpy import vstack
from matplotlib import pyplot as plt

Nout = 500
maxk = Nout/4
sampleint = 0.5

lochan = 1200
hichan = 1800 

BW = 800.e6/(2048.)

data = np.load('allan_20db_q.npz')['arr_0']

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

theory = np.sqrt(2)/(np.sqrt(BW)*np.sqrt(time))
#theory = 1./np.sqrt(time)

#plt.plot(time, alanvar, marker='o')
#plt.plot(time, theory)
plt.plot(time, alanvar/theory)
plt.xscale('log')
plt.yscale('log')
plt.show()
