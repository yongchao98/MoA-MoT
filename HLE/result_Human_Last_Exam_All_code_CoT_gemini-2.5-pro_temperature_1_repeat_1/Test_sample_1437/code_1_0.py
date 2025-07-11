import math

# The Lindhard polarization function at zero external frequency and momentum transfer,
# Pi(q=0, w=0), is equal to the negative of the density of states at the Fermi
# energy, -g(E_F). This code calculates its numerical value in the natural units
# of the 3D homogeneous electron gas.

# 1. Define the constants in the natural unit system where h_bar = 1, k_F = 1, and E_F = 1.
h_bar = 1.0
k_F = 1.0

# 2. In this system of units, the electron mass 'm' is derived from the Fermi energy relation:
# E_F = (h_bar^2 * k_F^2) / (2*m)
# 1 = (1^2 * 1^2) / (2*m) which gives m = 0.5.
m = 0.5

# 3. The formula for the density of states at the Fermi energy is:
# g(E_F) = (m * k_F) / (pi^2 * h_bar^2)
# The Lindhard function is Pi(0,0) = -g(E_F).
Pi_00 = - (m * k_F) / (math.pi**2 * h_bar**2)

# 4. Print the final equation and the resulting numerical value.
print("The Lindhard function at q=0, w=0 is Pi(0,0) = -g(E_F).")
print("The density of states at the Fermi energy is g(E_F) = (m * k_F) / (pi^2 * h_bar^2).")
print("\nIn natural units for the electron gas, the parameters are:")
print(f"m = {m}")
print(f"k_F = {k_F}")
print(f"h_bar = {h_bar}")
print("\nSubstituting these values into the equation for Pi(0,0):")
print(f"Pi(0,0) = - ( {m} * {k_F} ) / ( pi^2 * {h_bar}^2 )")

print("\nThe final numerical value is:")
print(Pi_00)