from astropy import constants as const
import numpy as np
import matplotlib.pyplot as plt

# Define quantities
m = np.array([43.35, 23.02]) * const.M_sun.value  # Black hole masses
Mc = ((m[0] * m[1]) ** (3 / 5)) / ((m[0] + m[1]) ** (1 / 5))  # Chirp mass
M = m[0] + m[1]  # Total mass
G = const.G.value
c = const.c.value
tc = 0.5  # Time of coalescence (determines number of periods to plot)
t = np.linspace(0, tc, num=1000, endpoint=False)  # Time array
tau = tc-t

# Define arrays
omega = (5/(256*tau))**(3/8) * (G*Mc/(c**3))**(-5/8)  # Orbital frequency
v = (G*M*omega)**(1/3)  # Orbital velocity
r = (G*M)**(1/3) * omega**(-2/3)  # Orbital separation

# Find data at merger
rs1 = 2*G*m[0]/c**2
rs2 = 2*G*m[1]/c**2
rMerge = 1.7 * (rs1 + rs2)
tMerge = (rMerge*G**(-3/4)*M**(-1/3)*(Mc / c**3)**(-5/12))**4 * 5 / 256
omegaMerge = (5/(256*tMerge))**(3/8) * (G*Mc/(c**3))**(-5/8)
vMerge = (G*M*omegaMerge)**(1/3)

# Plotting
# Orbital frequency
plt.plot(t, omega, label="Orbital frequency")
plt.xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5], [-0.5, -0.4, -0.3, -0.2, -0.1, 0.0])
plt.axhline(omegaMerge, linestyle='--', color='red', label='Orbital frequency at merger')
plt.xlabel("Time before merger (s)")
plt.ylabel(r'Orbital frequency $(\frac{rad}{s})$')
plt.ylim(0, 700)
plt.legend()
plt.show()
# Orbital velocity
plt.plot(t, v, label="Orbital velocity")
plt.xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5], [-0.5, -0.4, -0.3, -0.2, -0.1, 0.0])
plt.axhline(vMerge, linestyle='--', color='red', label='Orbital velocity at merger')
plt.xlabel("Time before merger (s)")
plt.ylabel(r'Orbital velocity $\frac{m}{s}$')
plt.legend()
plt.show()
# Orbital separation
plt.plot(t, r, label='Orbital separation')
plt.xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5], [-0.5, -0.4, -0.3, -0.2, -0.1, 0.0])
plt.axhline(rMerge, linestyle='--', color='red', label='Separation at merger')
plt.legend()
plt.xlabel("Time before merger (s)")
plt.ylabel("Separation distance (m)")
plt.show()

