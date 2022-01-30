import numpy as np
from PyAstronomy import pyaC
import matplotlib.pyplot as plt
from astropy import constants as const
import math

tevent = 1126259462.422  # Mon Sep 14 09:50:45 GMT 2015

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
plt.xlabel("Time (s) since " + str(tevent))
plt.ylabel("Strain")
plt.show()
# Livingston
plt.plot(ligoDat[:, 0], ligoDat[:, 2], 'g', label="L1 strain")
plt.plot(NRDat[:, 0], NRDat[:, 1], ':', color='k', label='Matched NR waveform')
plt.hlines(0, -0.09, 0.006, color='grey')
plt.scatter(L1zeros[0:2], [0, 0], color='purple')
plt.scatter(L1zeros[4:7], [0, 0, 0], color='red')
plt.legend(loc="lower left")
plt.xlabel("Time (s) since " + str(tevent))
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
H1fgw = 1 / (2 * H1diff)
L1fgw = 1 / (2 * L1diff)

# Find frequency at merger at second last (max amplitude) calculated frequency
H1fgw_merge = H1fgw[-2]
L1fgw_merge = L1fgw[-2]
print("H1 frequency at merger: ", round(H1fgw_merge, 2), "Hz")
print("L1 frequency at merger: ", round(L1fgw_merge, 2), "Hz")

# Calculate frequencies to the power of -8/3
H1fgw_83 = H1fgw ** (-8 / 3)
L1fgw_83 = L1fgw ** (-8 / 3)

# Calculate times for each freq. as midpoints of zero crossings
H1fgw_times = np.zeros(12)
L1fgw_times = np.zeros(12)
for i in range(12):
    H1fgw_times[i] = np.average([H1zeros[i], H1zeros[i + 1]])
    L1fgw_times[i] = np.average([L1zeros[i], L1zeros[i + 1]])

# Calculate least squares linear fit
H1fgw_83_fit, H1_fit_residuals, _, _, _ = np.polyfit(H1fgw_times, H1fgw_83, 1, full=True)
L1fgw_83_fit, L1_fit_residuals, _, _, _ = np.polyfit(L1fgw_times, L1fgw_83, 1, full=True)
print('H1 Fit Residuals:', H1_fit_residuals[0])
print('L1 Fit Residuals:', L1_fit_residuals[0])

# Frequency plots
# Livingston
plt.plot(L1fgw_times, L1fgw_83_fit[0] * L1fgw_times + L1fgw_83_fit[1], color='red', label='Least Squares Fit')
plt.scatter(L1fgw_times, L1fgw_83, color='red', label='Livingston Data')
plt.figtext(0.21,0.55,r'$\frac{df^{-8/3}}{dt}$', fontsize=14)
plt.figtext(0.28,0.55,r'$=-3.4\times 10^{-4}$', fontsize=10)
plt.legend()
plt.xlabel("Time (s)")
plt.ylabel(r'(Frequency/Hz)$^{-8/3}$')
plt.show()
# Handford
plt.plot(H1fgw_times, H1fgw_83_fit[0] * H1fgw_times + H1fgw_83_fit[1], label='Least Squares Fit')
plt.scatter(H1fgw_times, H1fgw_83, label='Handford Data')
plt.figtext(0.21,0.72,r'$\frac{df^{-8/3}}{dt}$', fontsize=14)
plt.figtext(0.28,0.72,r'$=-4.5\times 10^{-4}$', fontsize=10)
plt.legend()
plt.xlabel("Time (s)")
plt.ylabel(r'(Frequency/Hz)$^{-8/3}$')
plt.show()

# Calculate frequency-time slopes (df^(-8/3)/dt)
H1_f_slope = H1fgw_83_fit[0]
L1_f_slope = L1fgw_83_fit[0]
print("H1 frequency-time slope:", H1_f_slope)
print("L1 frequency-time slope:", L1_f_slope)

# Calculate total masses (in units of solar masses)
comp_ratio = 1.7
M_tot_H1 = const.c ** (3) / (2 * math.sqrt(2) * np.pi * H1fgw_merge * const.G * comp_ratio ** (3 / 2))
M_tot_H1_M0 = (M_tot_H1 / const.M_sun).value
M_tot_L1 = const.c ** (3) / (2 * math.sqrt(2) * np.pi * L1fgw_merge * const.G * comp_ratio ** (3 / 2))
M_tot_L1_M0 = (M_tot_L1 / const.M_sun).value
M_tot = np.average([M_tot_H1_M0, M_tot_L1_M0])
print("Estimated Total Masses: H1:", round(M_tot_H1_M0, 3), "solar masses, L1:", round(M_tot_L1_M0, 3), "solar masses")
print("Average Estimated Total Mass:", round(M_tot, 3), "solar masses")

# Calculate chirp masses (in units of solar masses)
M_ch_H1 = (const.c ** 3 / const.G) * (-H1_f_slope * 5 * (8 * np.pi) ** (-8 / 3)) ** (3 / 5)
M_ch_H1_M0 = (M_ch_H1 / const.M_sun).value
M_ch_L1 = (const.c ** 3 / const.G) * (-L1_f_slope * 5 * (8 * np.pi) ** (-8 / 3)) ** (3 / 5)
M_ch_L1_M0 = (M_ch_L1 / const.M_sun).value
M_ch = np.average([M_ch_H1_M0, M_ch_L1_M0])
print("Estimated Chirp Masses: H1:", round(M_ch_H1_M0, 3), "solar masses, L1:", round(M_ch_L1_M0, 3), "solar masses")
print("Average Estimated Chirp Mass:", round(M_ch, 3), "solar masses")

# Find m1 and m2 from total and chirp mass
# Solve quadratic of m2=M_tot-m1 subbed into M_ch equation:
m = np.roots([1, -M_tot, (M_ch ** 5 * M_tot) ** (1 / 3)])
print("Estimated Black Hole Masses:", round(m[0], 2), "and", round(m[1], 2), "solar masses")
