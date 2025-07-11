import math

# Description of the problem parameters
# Number of cavities (sites)
N = 7
# Number of photons (particles)
M = 4

# Step 1: The on-site energy term gives a constant offset.
# The total number of photons M is 4.
# This energy is M * omega.
omega_coefficient = M

# Step 2: Calculate the ground state energy from the hopping term.
# In the U->infinity limit, this system is equivalent to 4 free spinless fermions on a 7-site ring.
# The single-particle energies are E_q = -2*J*cos(2*pi*q/N).
# To find the ground state energy, we sum the M=4 lowest single-particle energies.

# The energy levels for N=7 are, in increasing order:
# E_0 (for q=0)
# E_1 = E_6 (for q=1 and q=6, degenerate)
# E_2 = E_5 (for q=2 and q=5, degenerate)
# ...and so on.

# The M=4 fermions will occupy the states with the lowest energies.
# These correspond to q=0, q=1, q=6, and one of q=2 or q=5. We'll pick q=2.
# The choice does not change the total energy as E_2 = E_5.

# The occupied states' energies are E_0, E_1, E_6, E_2.
# Their contribution to the total energy is the sum of these energies.
# E_hopping = E_0 + E_1 + E_6 + E_2
# E_0 = -2*J*cos(0) = -2J
# E_1 = -2*J*cos(2*pi/7)
# E_6 = -2*J*cos(12*pi/7) = -2*J*cos(2*pi/7)
# E_2 = -2*J*cos(4*pi/7)
# The total coefficient for J is therefore -2*(cos(0) + cos(2pi/7) + cos(12pi/7) + cos(4pi/7))

# Calculate the values of the cosine terms
cos_0 = math.cos(0)
cos_2pi_over_7 = math.cos(2 * math.pi / 7)
cos_4pi_over_7 = math.cos(4 * math.pi / 7)

# The coefficient for J is the sum of the coefficients from the four lowest energy levels.
# Since E_q's coefficient is -2*cos(2*pi*q/7), we have:
j_coefficient = -2 * (cos_0 + 2 * cos_2pi_over_7 + cos_4pi_over_7)

# Step 3: Print the final equation for the total ground state energy E_g.
# E_g = (energy from on-site term) + (energy from hopping term)
# E_g = omega_coefficient * omega + j_coefficient * J

print("The ground state energy E_g of the system is given by the equation:")
print(f"E_g = {omega_coefficient} * \u03C9 + ({j_coefficient:.4f}) * J")
print("\nIn a more compact form:")
print(f"E_g = {omega_coefficient}\u03C9 {j_coefficient:+.4f}J")
