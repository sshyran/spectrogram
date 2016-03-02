#!/usr/bin/ipython -pylab

# Copyright 2010, Bill Cox, all rights reserved.  This program may be freely
# distributed under the terms of the GPL license, version 2.

from pylab import *
from scipy import *

sampleRate = 10000

def computeNoiseBandwidth(response):
    total = 0.0
    peak = -1e50
    for i in range(sampleRate/2):
        if response[i] > peak:
	    peak = response[i]
    peak = peak*peak
    for i in range(sampleRate/2):
        total += response[i]*response[i]
    return total/peak

def genResponse(freq):
    signal=zeros(sampleRate*2)
    for i in range(sampleRate*2):
        signal[i] = sin(2*pi*freq*i/sampleRate)
    
    window=zeros(sampleRate*2)
    for i in range(sampleRate*2):
        window[i] = (1 - cos(pi*i/sampleRate))/2
    
    # Apply the window function to the signal
    windowedSignal = window*signal
    
    # Now do the overlap-add
    fftIn = windowedSignal[:sampleRate] + windowedSignal[sampleRate:]
    
    # Compute the FFT
    fftOut = fft(fftIn)
    return abs(fftOut[:sampleRate/2])

worstValue = -1e50
worstFreq = 0.0
for i in range(100):
    freq = 1000.0 + i/100.0
    response = genResponse(freq)
    testFreq = 1.1*freq
    if response[testFreq] > worstValue:
	worstValue = response[testFreq]
        worstFreq = freq
	print "Worst freq = %f" % worstFreq

response = genResponse(worstFreq)
plot(20*log10(response[900:1100]))
response = genResponse(1000)
B = computeNoiseBandwidth(response)
print "B = %f" % B
