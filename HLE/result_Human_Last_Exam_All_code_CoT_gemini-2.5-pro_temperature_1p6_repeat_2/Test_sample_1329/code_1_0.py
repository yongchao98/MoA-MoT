import math

# This script calculates the maximum gyroradius for a beta particle (electron) with 1 MeV
# of kinetic energy in a 166 mT magnetic field. This helps to verify the physical
# plausibility of the parameters given in the problem's answer choices.

# 1. Define physical constants and problem parameters
q = 1.60217663e-19      # Elementary charge in Coulombs
m_e_rest_MeV = 0.510998 # Electron rest mass in MeV/c^2
c = 299792458           # Speed of light in m/s
T_max_MeV = 1.0         # Maximum kinetic energy in MeV from the problem
B_T = 0.166             # Magnetic field in Tesla (166 mT)

# 2. Perform relativistic calculations for momentum
# The total energy (E) is the sum of rest energy and kinetic energy (T)
E_total_MeV = T_max_MeV + m_e_rest_MeV

# From the relativistic energy-momentum relation E^2 = (pc)^2 + (m_e*c^2)^2,
# we can find the momentum (p).
# pc = sqrt(E_total^2 - (m_e*c^2)^2)
p_times_c_MeV = math.sqrt(E_total_MeV**2 - m_e_rest_MeV**2)

# Convert momentum from units of MeV/c to SI units (kg*m/s)
# p = (p_times_c_MeV * 1e6 * q) / c
p_si = (p_times_c_MeV * 1e6 * q) / c

# 3. Calculate the gyroradius
# The radius of the spiral path (gyroradius) is given by r = p_perp / (q * B).
# We calculate the maximum possible radius, where all momentum is perpendicular to the B-field.
gyroradius_m = p_si / (q * B_T)

# 4. Print the results in a clear format
print("Calculation of maximum gyroradius for a 1 MeV electron in a 166 mT magnetic field.")
print("-" * 75)
print(f"Given Kinetic Energy (T): {T_max_MeV} MeV")
print(f"Given Magnetic Field (B): {B_T} T")
print("-" * 75)
print(f"Relativistic Momentum (p): {p_si:.4e} kg*m/s")
print(f"Electron Charge (q):       {q:.4e} C")
print("-" * 75)
print("Final Equation for Gyroradius (r): r = p / (q * B)")
print(f"r = {p_si:.4e} / ({q:.4e} * {B_T})")
print(f"r = {gyroradius_m:.4f} meters")
print(f"r = {gyroradius_m * 100:.2f} cm")
print("-" * 75)
print("A maximum gyroradius of a few centimeters is a manageable scale for a tabletop experiment,")
print("confirming that the proposed magnetic field strength is physically reasonable for this application.")
