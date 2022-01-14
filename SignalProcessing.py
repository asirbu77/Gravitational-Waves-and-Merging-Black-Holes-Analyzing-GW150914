##
# Taken from: https://www.gw-openscience.org/GW150914data/GW150914_tutorial.html
##

# Standard python numerical analysis imports:
import numpy as np
from scipy import signal
from scipy.interpolate import interp1d
from scipy.signal import butter, filtfilt, iirdesign, zpk2tf, freqz

# the ipython magic below must be commented out in the .py file, since it doesn't work.
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import h5py

# LIGO-specific readligo.py
import readligo as rl


# generate linear time-domain filter coefficients, common to both H1 and L1.
# First, define some functions:

# This function will generate digital filter coefficients for bandstops (notches).
# Understanding it requires some signal processing expertise, which we won't get into here.
def iir_bandstops(fstops, fs, order=4):
    """ellip notch filter
    fstops is a list of entries of the form [frequency (Hz), df, df2]
    where df is the pass width and df2 is the stop width (narrower
    than the pass width). Use caution if passing more than one freq at a time,
    because the filter response might behave in ways you don't expect.
    """
    nyq = 0.5 * fs

    # Zeros zd, poles pd, and gain kd for the digital filter
    zd = np.array([])
    pd = np.array([])
    kd = 1

    # Notches
    for fstopData in fstops:
        fstop = fstopData[0]
        df = fstopData[1]
        df2 = fstopData[2]
        low = (fstop - df) / nyq
        high = (fstop + df) / nyq
        low2 = (fstop - df2) / nyq
        high2 = (fstop + df2) / nyq
        z, p, k = iirdesign([low, high], [low2, high2], gpass=1, gstop=6,
                            ftype='ellip', output='zpk')
        zd = np.append(zd, z)
        pd = np.append(pd, p)

    # Set gain to one at 100 Hz...better not notch there
    bPrelim, aPrelim = zpk2tf(zd, pd, 1)
    outFreq, outg0 = freqz(bPrelim, aPrelim, 100 / nyq)

    # Return the numerator and denominator of the digital filter
    b, a = zpk2tf(zd, pd, k)
    return b, a


def get_filter_coefs(fs):
    # assemble the filter b,a coefficients:
    coefs = []

    # bandpass filter parameters
    lowcut = 43
    highcut = 260
    order = 4

    # bandpass filter coefficients
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    bb, ab = butter(order, [low, high], btype='band')
    coefs.append((bb, ab))

    # Frequencies of notches at known instrumental spectral line frequencies.
    # You can see these lines in the ASD above, so it is straightforward to make this list.
    notchesAbsolute = np.array(
        [14.0, 34.70, 35.30, 35.90, 36.70, 37.30, 40.95, 60.00,
         120.00, 179.99, 304.99, 331.49, 510.02, 1009.99])

    # notch filter coefficients:
    for notchf in notchesAbsolute:
        bn, an = iir_bandstops(np.array([[notchf, 1, 0.1]]), fs, order=4)
        coefs.append((bn, an))

    # Manually do a wider notch filter around 510 Hz etc.
    bn, an = iir_bandstops(np.array([[510, 200, 20]]), fs, order=4)
    coefs.append((bn, an))

    # also notch out the forest of lines around 331.5 Hz
    bn, an = iir_bandstops(np.array([[331.5, 10, 1]]), fs, order=4)
    coefs.append((bn, an))

    return coefs


# and then define the filter function:
def filter_data(data_in, coefs):
    data = data_in.copy()
    for coef in coefs:
        b, a = coef
        # filtfilt applies a linear filter twice, once forward and once backwards.
        # The combined filter has linear phase.
        data = filtfilt(b, a, data)
    return data


# ----------------------------------------------------------------
# Load LIGO data from a single file
# ----------------------------------------------------------------
# First from H1
fn_H1 = 'H-H1_LOSC_4_V2-1126259446-32.hdf5'
strain_H1, time_H1, chan_dict_H1 = rl.loaddata(fn_H1, 'H1')
# and then from L1
fn_L1 = 'L-L1_LOSC_4_V2-1126259446-32.hdf5'
strain_L1, time_L1, chan_dict_L1 = rl.loaddata(fn_L1, 'L1')

# sampling rate:
fs = 4096
# both H1 and L1 will have the same time vector, so:
time = time_H1

# read in the NR template
NRtime, NR_H1 = np.genfromtxt('GW150914_4_NR_waveform.txt').transpose()

tevent = 1126259462.422         # Mon Sep 14 09:50:45 GMT 2015

# get filter coefficients
coefs = get_filter_coefs(fs)

# filter the data:
strain_H1_filt = filter_data(strain_H1, coefs)
strain_L1_filt = filter_data(strain_L1, coefs)
# filter NR template as we do with the data:
NR_H1_filt = filter_data(NR_H1, coefs)

# plot the data prior to filtering:
plt.figure()
plt.plot(time-tevent,strain_H1,'r',label='H1 strain')
plt.plot(time-tevent,strain_L1,'g',label='L1 strain')
plt.xlim([-0.2,0.1])
plt.xlabel('time (s) since '+str(tevent))
plt.ylabel('strain')
plt.legend(loc='lower right')
plt.title('aLIGO strain data near GW150914')
plt.savefig('GW150914_H1_strain_unfiltered.png')

# plot the data after filtering:
# first, shift L1 by 7 ms, and invert. See the GW150914 detection paper for why!
strain_L1_fils = -np.roll(strain_L1_filt,int(0.007*fs))
# We also have to shift the NR template by 2 ms to get it to line up properly with the data
plt.figure()
plt.plot(time-tevent,strain_H1_filt,'r',label='H1 strain')
plt.plot(time-tevent,strain_L1_fils,'g',label='L1 strain')
#plt.plot(NRtime+0.002,NR_H1_filt,'k',label='Matched NR waveform')
plt.xlim([-0.2,0.1])
plt.ylim([-1.5e-21,1.5e-21])
plt.xlabel('Time (s) since '+str(tevent))
plt.ylabel('Strain')
plt.legend(loc='lower left')
plt.title('aLIGO FILTERED strain data near GW150914')
plt.savefig('GW150914_H1_strain_filtered.png')

# Arrange LIGO data into matrix
ligoDat = np.zeros((131072, 3))
ligoDat[:, 0] = time-tevent
ligoDat[:, 1] = strain_H1_filt
ligoDat[:, 2] = strain_L1_fils

# Arrange numerical relativity waveform data into matrix
NRDat = np.zeros((2769, 2))
NRDat[:, 0] = NRtime
NRDat[:, 1] = NR_H1_filt

# Save as csv file
np.savetxt("ligoDat.csv", ligoDat, delimiter=",")
np.savetxt("NRDat.csv", NRDat, delimiter=",")


