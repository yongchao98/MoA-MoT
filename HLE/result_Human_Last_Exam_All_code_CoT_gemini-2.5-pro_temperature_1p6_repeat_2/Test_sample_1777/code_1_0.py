# Define the given physical constants
band_gap_Eg = 3.0  # Band gap in eV
resonance_peak_E1s = 1.0  # 1s exciton resonance peak in eV
principal_quantum_number_n = 3

print("This script calculates the binding energy of the n=3 exciton state.")
print("---")

# Step 1: Calculate the exciton Rydberg energy (Ry).
# This is the binding energy of the ground state (n=1) and is the difference
# between the band gap and the ground state energy level.
print("Step 1: Calculate the exciton Rydberg energy (Ry).")
print("Formula: Ry = Eg - E_1s")
rydberg_energy_Ry = band_gap_Eg - resonance_peak_E1s
print(f"Ry = {band_gap_Eg} eV - {resonance_peak_E1s} eV = {rydberg_energy_Ry} eV")
print("---")


# Step 2: Calculate the binding energy for the n=3 state.
# For a hydrogen-like system, the binding energy E_b(n) scales as Ry / n^2.
# The phrase 'Rydberg energy for n=3' is interpreted as this binding energy.
print("Step 2: Calculate the binding energy for n=3 (E_b(3)).")
print("Formula: E_b(n) = Ry / n^2")
n_squared = principal_quantum_number_n**2
binding_energy_n3 = rydberg_energy_Ry / n_squared

# Output each number in the final equation as requested
print(f"E_b({principal_quantum_number_n}) = {rydberg_energy_Ry} eV / {principal_quantum_number_n}^2")
print(f"E_b({principal_quantum_number_n}) = {rydberg_energy_Ry} eV / {n_squared}")
print(f"E_b({principal_quantum_number_n}) = {binding_energy_n3:.4f} eV")
print("---")

# Final result
print(f"The final calculated Rydberg energy for n=3 is: {binding_energy_n3:.4f} eV")