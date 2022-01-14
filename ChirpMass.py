import numpy as np
from PyAstronomy import pyaC
import matplotlib.pyplot as plt

tevent = 1126259462.422         # Mon Sep 14 09:50:45 GMT 2015

# Write data from csv's into numpy array
ligoDat = np.genfromtxt('ligoDat.csv', delimiter=',')
NRDat = np.genfromtxt('NRDat.csv', delimiter=',')
# Remove area outside interval of interest
ligoPos1 = np.argmax(ligoDat > -0.09, axis=0)[0]
ligoPos2 = np.argmax(ligoDat > 0.006, axis=0)[0]
NRPos1 = np.argmax(NRDat > -0.09, axis=0)[0]
NRPos2 = np.argmax(NRDat > 0.006, axis=0)[0]
ligoDat = ligoDat[ligoPos1:ligoPos2, :]
NRDat = NRDat[NRPos1:NRPos2, :]

# Find zero crossings
H1zeros = np.array(pyaC.zerocross1d(ligoDat[:, 0], ligoDat[:, 1]))
L1zeros = np.array(pyaC.zerocross1d(ligoDat[:, 0], ligoDat[:, 2]))

# Area of interest plot
# Hanford
plt.plot(ligoDat[:, 0], ligoDat[:, 1], 'r', label="H1 strain")
plt.plot(NRDat[:, 0], NRDat[:, 1], ':', color='k', label='Matched NR waveform')
plt.hlines(0, -0.09, 0.006, color='grey')
plt.scatter(H1zeros[2:5], [0, 0, 0], color='cyan')
plt.legend(loc="lower left")
plt.xlabel("Time (s) since "+str(tevent))
plt.ylabel("Strain")
plt.show()
# Livingston
plt.plot(ligoDat[:, 0], ligoDat[:, 2], 'g', label="L1 strain")
plt.plot(NRDat[:, 0], NRDat[:, 1], ':', color='k', label='Matched NR waveform')
plt.hlines(0, -0.09, 0.006, color='grey')
plt.scatter(L1zeros[0:2], [0, 0], color='magenta')
plt.scatter(L1zeros[4:7], [0, 0, 0], color='yellow')
plt.legend(loc="lower left")
plt.xlabel("Time (s) since "+str(tevent))
plt.ylabel("Strain")
plt.show()

# Average values to better fit numerical relativity template
L1_avg_1 = np.average(L1zeros[0:2])
L1_avg_2 = np.average(L1zeros[4:7])
H1_avg_1 = np.average(H1zeros[2:5])
H1zeros = np.delete(H1zeros, [2, 3, 4])
H1zeros = np.insert(H1zeros, 2, H1_avg_1)
L1zeros = np.delete(L1zeros, [0, 1])
L1zeros = np.delete(L1zeros, [2, 3, 4])
L1zeros = np.insert(L1zeros, 0, L1_avg_1)
L1zeros = np.insert(L1zeros, 3, L1_avg_2)

# Calculate the time intervals between consecutive zero crossings
H1diff = np.diff(H1zeros)
L1diff = np.diff(L1zeros)

# Calculate GW frequencies from time intervals
H1fgw = 1/(2*H1diff)
L1fgw = 1/(2*L1diff)

H1fgw_83 = H1fgw**(-8/3)
L1fgw_83 = L1fgw**(-8/3)

H1fgw_83_fit = np.polyfit()
plt.plot(L1fgw_83)
plt.show()
