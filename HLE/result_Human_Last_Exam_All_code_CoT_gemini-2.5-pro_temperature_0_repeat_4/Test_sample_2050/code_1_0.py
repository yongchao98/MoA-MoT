import math

# The problem reduces to finding the ground state energy of 4 non-interacting fermions on a 7-site ring.
# The single-particle energies are given by E_m = -2*J*cos(2*pi*m/N) for N=7.
# We need to fill the 4 lowest energy levels.
# The quantum numbers for the levels are m = 0, 1, 2, 3, 4, 5, 6.
# The energies are ordered as E_0 < E_1=E_6 < E_2=E_5 < E_3=E_4.
# The 4 fermions occupy the states with m=0, m=1, m=6, and one of m=2 or m=5.
# The energy contribution from the hopping term is E' = E_0 + E_1 + E_6 + E_2.
# E' = -2*J*cos(0) - 2*J*cos(2*pi/7) - 2*J*cos(2*pi*6/7) - 2*J*cos(2*pi*2/7)
# Since cos(2*pi*6/7) = cos(2*pi/7), this simplifies to:
# E' = -2*J * (1 + 2*cos(2*pi/7) + cos(4*pi/7))

# Number of photons (and sites for the constant energy term)
num_photons = 4

# Calculate the coefficient of J
# cos(2*pi/7)
cos_2pi_over_7 = math.cos(2 * math.pi / 7)
# cos(4*pi/7)
cos_4pi_over_7 = math.cos(4 * math.pi / 7)

# Coefficient C such that E' = C * J
j_coefficient = -2 * (1 + 2 * cos_2pi_over_7 + cos_4pi_over_7)

# The total ground state energy is E = 4*omega + E'
# We print the numbers in the final equation for the energy.
print("The ground state energy is given by the equation:")
print(f"E_g = {num_photons} * omega + ({j_coefficient}) * J")
print("\nThe numbers in the final equation are:")
print(f"Coefficient of omega: {num_photons}")
print(f"Coefficient of J: {j_coefficient}")
