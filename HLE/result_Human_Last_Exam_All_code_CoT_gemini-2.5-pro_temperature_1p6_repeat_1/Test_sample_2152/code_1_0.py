import numpy as np

# Physical constants
e = 1.60217663e-19  # Elementary charge in Coulombs
hbar = 1.05457182e-34  # Reduced Planck constant in J.s
m_e = 9.1093837e-31  # Electron mass in kg
c = 299792458  # Speed of light in m/s
mu_B = 9.27401008e-24 # Bohr magneton in J/T

# Given parameters
# The problem gives relations e_m = e/16 and mu = 25*mu_B.
# This implies the quantities e_m and mu are related to the elementary charge e
# and fundamental constants. We will substitute these relations into our ratio formula.
# e_m = e/16
# mu = 25 * mu_B = 25 * e*hbar/(2*m_e)
# The particle is charged, we assume it's an electron, so its mass m = m_e.

v = c / 100.0  # Particle speed
theta = np.pi / 30.0  # Scattering angle

# Calculate momentum p and momentum transfer q
p = m_e * v
q = 2 * p * np.sin(theta / 2.0)

# Ratio formula R = (3 * m^2 * e_m^2) / (mu^2 * q^4)
# We substitute the expressions for e_m and mu in terms of fundamental constants.
# e_m^2 = (e/16)^2
# mu^2 = (25 * e*hbar / (2*m_e))^2
# m^2 = m_e^2

# So, the ratio R becomes:
# R = (3 * m_e^2 * (e/16)^2) / ((25 * e*hbar / (2*m_e))^2 * q^4)
# R = (3 * m_e^2 * e^2 / 256) / (625 * e^2 * hbar^2 / (4*m_e^2) * q^4)
# The e^2 terms cancel out.
# R = (3 * m_e^2 / 256) * (4 * m_e^2 / (625 * hbar^2)) / q^4
# R = (12 * m_e^4) / (160000 * hbar^2 * q^4)
# R = (3 * m_e^4) / (40000 * hbar^2 * q^4)

# Substitute q = 2 * m_e * v * sin(theta/2)
# q^4 = 16 * m_e^4 * v^4 * sin(theta/2)^4
# R = (3 * m_e^4) / (40000 * hbar^2 * (16 * m_e^4 * v^4 * np.sin(theta/2)**4))
# The m_e^4 terms cancel out.
# R = 3 / (640000 * hbar^2 * v^4 * np.sin(theta/2)**4)

numerator = 3
denominator = 640000 * hbar**2 * v**4 * (np.sin(theta/2.0))**4
ratio = numerator / denominator

print(f"The calculation for the ratio proceeds as follows:")
print(f"Numerator in the simplified ratio formula: {numerator}")
print(f"Denominator value: {denominator:.4e}")
print(f"Final ratio of the differential cross-sections (Monopole/Dipole) is: {ratio:.4e}")

# The result is a large number, but this is the consequence of the derived formulas.
# Let's show the equation with the final numbers plugged in.
print("\nThe final equation with all values plugged in is:")
print(f"Ratio = (3 * ({m_e:.4e})^2 * ({e/16:.4e})^2) / (({25*mu_B:.4e})^2 * ({q:.4e})^4)")
final_ratio_from_raw = (3 * m_e**2 * (e/16)**2) / ((25*mu_B)**2 * q**4)
print(f"Calculated Ratio = {final_ratio_from_raw:.4e}")
