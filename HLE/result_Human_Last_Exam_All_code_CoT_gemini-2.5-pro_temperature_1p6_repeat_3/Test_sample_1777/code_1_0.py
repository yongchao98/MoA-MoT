# Step 1: Define the given physical constants.
band_gap_Eg = 3.0  # Band gap in eV
peak_energy_1s = 1.0  # Resonance peak energy for the 1s exciton in eV
quantum_number_n = 3 # The principal quantum number for the target state

# Step 2: Calculate the binding energy of the 1s ground state exciton (E_b_1s).
# The energy of an exciton peak is the band gap minus its binding energy.
# E_peak = E_g - E_b  =>  E_b = E_g - E_peak
print("Step 1: Calculate the ground state (1s) exciton binding energy (E_b_1s).")
print(f"The formula is: E_b_1s = Band Gap - 1s Peak Energy")
binding_energy_1s = band_gap_Eg - peak_energy_1s
print(f"E_b_1s = {band_gap_Eg} eV - {peak_energy_1s} eV = {binding_energy_1s} eV")
print("-" * 30)

# Step 3: Calculate the binding energy for the n=3 state using the 2D exciton model.
# For a 2D system, the binding energies of the exciton states follow the series:
# E_b_n = E_b_1s / (2*n - 1)^2
print(f"Step 2: Calculate the binding energy for n = {quantum_number_n} using the 2D exciton model.")
print("The formula is: E_b_n = E_b_1s / (2*n - 1)^2")

# Calculate the denominator first for clarity
denominator = (2 * quantum_number_n - 1)**2
binding_energy_n = binding_energy_1s / denominator

# Print the calculation with all numbers included
print("\nFinal Calculation:")
print(f"E_b_{quantum_number_n} = {binding_energy_1s} eV / (2 * {quantum_number_n} - 1)^2")
print(f"E_b_{quantum_number_n} = {binding_energy_1s} eV / ({int(denominator**0.5)})^2")
print(f"E_b_{quantum_number_n} = {binding_energy_1s} eV / {denominator}")
print(f"E_b_{quantum_number_n} = {binding_energy_n:.3f} eV")

# The final answer is the binding energy for n=3.
print(f"\nThe Rydberg energy for n = {quantum_number_n} is {binding_energy_n:.3f} eV.")