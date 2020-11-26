# -*- coding: utf-8 -*-
"""
Created on Sun Dec  1 18:01:06 2019

@author: Ramadan
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Nov 30 13:53:17 2019

@author: Ramadan
"""

import pandas as pd
from scipy.signal import butter, lfilter,filtfilt,firls,firwin
from scipy import fft, arange
from scipy import signal
import math
import numpy as np
from scipy.optimize import leastsq
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import freqz


def butter_highpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='high', analog=False)
    return b, a

def butter_highpass_filter(data, cutoff, fs, order=5):
    b, a = butter_highpass(cutoff, fs, order=order)
    y = filtfilt(b, a, data)
    return y

def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = filtfilt(b, a, data)
    return y 

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

def plotSpectrum(y,Fs):
 """
 Plots a Single-Sided Amplitude Spectrum of y(t)
 """
 n = len(y) # length of the signal
 k = arange(n)
 T = n/Fs
 frq = k/T # two sides frequency range
 frq = frq[range(n//2)] # one side frequency range
 frq = frq-2.087e9
 Y = fft(y)/n # fft computing and normalization
 Y = Y[range(n//2)]
 plt.figure(figsize=(120,30))

 plt.plot(frq,abs(Y),'r') # plotting the spectrum
 plt.xlabel('Freq (Hz)')
 plt.grid()
 plt.ylabel('|Y(freq)|')
   
def residuals(p, y):  
    A,k,theta = p  
    err = y-A*np.sin(2*np.pi*k+theta)  
    return err  

def eye(bIeye,bQeye,rx,n=60):
    plt.figure()
    for k in range(rx):
        plt.plot(bIeye[k,0:n])
        plt.plot(bQeye[k,0:n])

def constellation(bIeye,bQeye,s=1,n=150):
    plot_lims = [-s,s]
    plt.figure(figsize=(15,14))
    plt.plot(bIeye[0,0:n], bQeye[0,0:n], '.')
    plt.xlim(plot_lims)
    plt.ylim(plot_lims)
    plt.title('QPSK constellation at an SNR of  dB')
    plt.xlabel('real part')
    plt.ylabel('imaginary part')

def lowpass_freq_resp(cutoff,Fs,Ts,lowpass_order = 5):
    lowpass_delay = (lowpass_order // 2)/Fs  # a lowpass of order N delays the signal by N/2 samples (see plot)
    # design the filter
    lowpass =firwin(lowpass_order, cutoff/(Fs/2))
    t_lp = np.arange(len(lowpass))/Fs
    f_lp = np.linspace(-Fs/2, Fs/2, 2048, endpoint=False)
    H = np.fft.fftshift(np.fft.fft(lowpass, 2048))
    
    plt.subplot(121)
    plt.plot(t_lp/Ts, lowpass)
    plt.gca().annotate(r'$\tau_{LP}$', xy=(lowpass_delay/Ts,0.08), xytext=(lowpass_delay/Ts+0.3, 0.08), arrowprops=dict(arrowstyle='->'))
    
    plt.subplot(122)
    plt.plot(f_lp, 20*np.log10(abs(H)))

def order_response():
    plt.figure(1)
    plt.clf()
    for order in [3, 5,6, 9]:
        #b, a = butter_bandpass(lowcut, highcut, fs, order=order)
        b, a = butter_lowpass(lowcut, fs, order=order)
        w, h = freqz(b, a, worN=2000)
        plt.plot((fs * 0.5 / np.pi) * w, abs(h), label="order = %d" % order)

    plt.plot([0, 0.5 * fs], [np.sqrt(0.5), np.sqrt(0.5)],
             '--', label='sqrt(0.5)')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Gain')
    plt.grid(True)
    plt.legend(loc='best')
    
def costaLoop(st,fc,fs):
    N = len(st)
    t = np.linspace(0,N,N)
    phi = np.zeros([N,1])
    s1 = np.zeros([N,1])
    s2 = np.zeros([N,1])
    y1 = np.zeros([N,1])
    y2 = np.zeros([N,1])
    for i in range(0,N-1):
        if i>0:
            phi[i] = phi[i-1] - (0)*np.pi*np.sign(y1[i-1]*y2[i-1])
            
        #s1[i] = st[i] * np.cos(2*np.pi*fc*t[i]/fs  + phi[i]);
        #s2[i] = st[i] * np.sin(2*np.pi*fc*t[i]/fs  + phi[i]);
        s1[i] = st[i]* np.sin(2*np.pi*fc*t[i]/fs+phi[i])
        s2[i] = st[i]* -np.cos(2*np.pi*fc*t[i]/fs+phi[i])
        if i<100:
            for j in range(0,i):
                y1[i] = y1[i]+s1[j]
                y2[i] = y2[i]+s2[j]
        else:
            for j in range(i-100,i):
                y1[i] = y1[i]+s1[j]
                y2[i] = y2[i]+s2[j]
    
     #plt.figure()
     #plt.plot(phi)            
    return y1,y2
        
def downsample(smp, scale=1.0):
    """
    Keyword arguments:
    scale - scale factor for length of sound (2.0 means double length)
    """
    # calculate new length of sample
    n = round(len(smp) / scale)
    # use linear interpolation
    # endpoint keyword means than linspace doesn't go all the way to 1.0
    # If it did, there are some off-by-one errors
    # e.g. scale=2.0, [1,2,3] should go to [1,1.5,2,2.5,3,3]
    # but with endpoint=True, we get [1,1.4,1.8,2.2,2.6,3]
    # Both are OK, but since resampling will often involve
    # exact ratios (i.e. for 44100 to 22050 or vice versa)
    # using endpoint=False gets less noise in the resampled sound
    return np.interp(
        np.linspace(0.0, 1.0, n, endpoint=False), # where to interpret
        np.linspace(0.0, 1.0, len(smp), endpoint=False), # known positions
        smp, # known data points
        )
    
if __name__ == "__main__":

    # Sample rate and desired cutoff frequencies (in Hz).
    fs = 20.0e9 #sampling freq
    fc = 2.087e9 #carrier freq
    #lowcut = 1e9
    #highcut = 3e9
    rs = 1e9
    sps = int(fs/rs) #sample per symbol
    ts = 1/fs
    q = 4
    # Plot the frequency response for a few different orders.
    lowcut = 1e9# rs/2
    highcut = 3e9# 2*fc - rs/2

    #file = "E:/phd/Jacqueline/1Gbps_Data_20Gbps.csv"
    file = "data/signal.csv"
    #file = "data/iq1gbps_20gsps_PD_detected.csv"
    #file = "data/iq1gbps_20gsps_to_MZM.csv"
    
    dataset = pd.read_csv(file,header=None)
    x = dataset.iloc[:,1]
    t = dataset.iloc[:,0]
    
    #plotSpectrum(x,fs)
    
    # Filter a noisy signal.
    p = butter_lowpass_filter(x,highcut,fs,order=5)
    y = butter_bandpass_filter(p, lowcut, highcut, fs,order=7)

    #y = butter_highpass_filter()
    #plotSpectrum(p,fs)
    #plotSpectrum(y,fs)
    
    #demodulation
    
    
    
    
    phi1 = 0
    phi2 = 0
    A = 1
    iconv = A* np.sin(2*np.pi*fc*t+phi1)
    qconv = A* -np.cos(2*np.pi*fc*t+phi2)
    rI =  iconv*y
    rQ =  qconv*y
    
    
    #rI,rQ = costaLoop(y,fc,fs)
    
    print(len(rI))
    print(len(rQ))
    
    plt.figure(figsize=(28,7))
    plt.plot(rI[:1000])
    plt.plot(rQ[:1000])

    fpass = 9e8 #rs/2
    fstop = 3e9 #2*fc - rs/2
    wpass = 1
    wstop = 1
    n = 51
    bl = firls(n,np.divide([0,fpass,fstop,fs/2],(fs/2)),[1,1,0,0],[wpass,wstop])
    
   
    wI = filtfilt(bl,1,rI)
    wQ = filtfilt(bl,1,rQ)
       

    dt = np.vstack((wI,wQ)).T
    np.savetxt("data/data.csv",dt,delimiter=',')
    '''
    cutoff = 3.0e9
    #filtering
    b,a  = butter_lowpass(cutoff,fs,11)
    wI =  lfilter(b,a,wI)
    wQ =  lfilter(b,a,wQ)


    #plt.legend()
    w = wI+wQ
'''
    
 #   plt.figure()
  #  plt.plot(wI[:1000],label='i2')
 #   plt.plot(wQ[:1000],label='q2')
    
#    plt.figure()
#    fig, axs = plt.subplots(6, sharex=True, sharey=True)
#    axs[0].plot(x[:1000])
#    axs[1].plot(y[:1000])
#    axs[2].plot(rI[:1000])
#    axs[3].plot(wI[:1000])
#    axs[4].plot(wQ[:1000])
#    axs[5].plot(w[:1000])
    
    #plotSpectrum(rI,fs)
    #plotSpectrum(wI,fs)
    



    
    #uI = signal.decimate(wI,10)
    #uQ = signal.decimate(wQ,10)
    
    uI = wI#downsample(wI,sps)
    uQ = wQ#downsample(wQ,sps)
    #plt.plot(uI[:10])
    
    fctr = 4
    
    rx = int(sps* fctr)
    rn = int(math.floor(len(uI)/rx)*rx)
    ry = int(math.floor(len(uQ)/rx))

    
    plotSpectrum(rn,fs)
    bI = np.reshape(uI[0:rn],[rx,ry], order="F").copy()
    bQ = np.reshape(uQ[0:rn],[rx,ry], order="F").copy()
    
    
    og = pd.read_csv("data/signal.csv",header=None)
    x = og.iloc[:,1]
    
    #plt.figure(figsize=(28,7))
    #plt.plot(uI[:1000])
    #plt.plot(uQ[:1000])
    #plt.plot(w[:1000])
    #plt.plot(x[:1000])
    #fig, axs = plt.subplots(2, sharex=True, sharey=True)
    #fig.figsize = (28,7)
    #axs[0].plot(wI[:1000])
    #axs[1].plot(wQ[:1000])
    
    di = bI[range(0,rx)]
    dq = bQ[range(0,rx)]
    
    #plt.figure()
    #fig, axs = plt.subplots(2, sharex=True, sharey=True)
    #axs[0].plot(di)
    #axs[1].plot(dq)

    #plt.plot(di)
    #plt.plot(dq)
    #eye(bI,bQ,rx,10)
    #constellation(bI,bQ)

    
    
    
    
