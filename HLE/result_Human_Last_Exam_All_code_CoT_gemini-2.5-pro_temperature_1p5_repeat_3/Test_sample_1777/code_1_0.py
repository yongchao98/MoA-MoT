import math

# Given parameters for the 2D semiconductor
band_gap_Eg = 3.0  # eV
resonance_peak_E1s = 1.0  # eV
quantum_number_n = 3

# Step 1: Calculate the binding energy of the 1s ground state exciton (E_b_1s)
# The resonance peak energy is related to the band gap and binding energy by: E_peak = E_gap - E_binding
E_b_1s = band_gap_Eg - resonance_peak_E1s

# Step 2: Calculate the effective 2D Rydberg constant (R*)
# The binding energy series for a 2D system is E_b(n) = R* / (n - 1/2)^2
# For the ground state (n=1), this gives E_b(1s) = R* / (1 - 0.5)^2 = R* / 0.25 = 4 * R*
# Therefore, R* = E_b(1s) / 4
R_star = E_b_1s / 4.0

# Step 3: Calculate the binding energy for the n=3 state (E_b_3)
# This is the "Rydberg energy for n=3" requested by the user.
E_b_n = R_star / (quantum_number_n - 0.5)**2

# Step 4: Print the result, showing the final equation with all numbers.
print(f"The Rydberg energy (binding energy) for n = {quantum_number_n} is calculated as follows:")
print(f"1. Find the 1s binding energy: E_b(1s) = E_gap - E_1s_peak = {band_gap_Eg} - {resonance_peak_E1s} = {E_b_1s} eV")
print(f"2. Find the effective Rydberg constant R*: R* = E_b(1s) / 4 = {E_b_1s} / 4 = {R_star} eV")
print(f"3. Calculate the binding energy for n={quantum_number_n}: E_b({quantum_number_n}) = R* / ({quantum_number_n} - 0.5)^2")
print(f"   E_b({quantum_number_n}) = {R_star} / ({quantum_number_n - 0.5})^2 = {E_b_n} eV")